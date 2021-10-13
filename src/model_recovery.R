require(pmwg)
require(rtdists)
library(dplyr)
library(MCMCpack)
library(stringi)
library(mcce)

print(sessionInfo())

# For debugging:
# Sys.setenv(DCE_EST_EXP="SymbolicVDCE", VDCE_DISPLAY="Absent", NCPUS=3, DCE_REC_MODEL="IEX", DCE_REC_MED="median_alpha_exp2_abs.RDS", DCE_MODEL_FILE="Task2_Absent_1069903.rcgbcm_CorrectedTry1.RData", DCE_REC_DATA="SymbolicVDCE_IEX_Absent_0P6QGThVC9il_untagged_data.RDS")
# Sys.setenv(DCE_EST_EXP="NumericVDCE", NCPUS=3, DCE_REC_MODEL="IEX", DCE_REC_MED="median_alpha_exp1.RDS", DCE_MODEL_FILE="NumericVDCE_1878182.rcgbcm_Estimation5Model.RData")
# Sys.setenv(DCE_EST_EXP="NumericVDCE", NCPUS=3, DCE_REC_MODEL="CB", DCE_REC_MED="median_alpha_exp1.RDS", DCE_MODEL_FILE="NumericVDCE_1878182.rcgbcm_Estimation5Model.RData", DCE_REC_DATA="5ModelRecovery/NumericVDCE_CB_1878514.rcgbcm_5ModelRecovery_data.RDS")
# Sys.setenv(DCE_EST_EXP="PrefDCE", NCPUS=3, DCE_REC_MODEL="CB", DCE_REC_MED="median_alpha_pref.RDS", DCE_MODEL_FILE="PrefDCE_1896523.rcgbcm_Estimation5Model.RData")
# Get environment variables to normal vars
known_vars <- c("DCE_EST_EXP", "VDCE_DISPLAY", "NCPUS", "PBS_JOBID", "VDCE_TAG",
                "DCE_REC_MODEL", "DCE_REC_MED", "DCE_MODEL_FILE", "DCE_REC_DATA")

envars <- Sys.getenv(known_vars)
print(envars)
envars <- as.list(envars)
experiment <- envars$DCE_EST_EXP
displaytype <- envars$VDCE_DISPLAY
cores <- ifelse(envars$NCPUS == "", 1, as.numeric(envars$NCPUS))
jobid <- ifelse(
  envars$PBS_JOBID == "",
  stri_rand_strings(1, 12),
  envars$PBS_JOBID
)
tag <- ifelse(envars$VDCE_TAG == "", "untagged", envars$VDCE_TAG)
test_model <- envars$DCE_REC_MODEL
median_file <- envars$DCE_REC_MED
model_file <- envars$DCE_MODEL_FILE
recovery_data <- envars$DCE_REC_DATA

get_estimation_data <- function(filename) {
  load(here::here("data", "output", filename), ex <- new.env())
  ex$sampled$data %>%
    tibble()
}

medians <- readRDS(here::here("data", "output", median_file))

if (test_model == "") {
  stop("DCE_REC_MODEL environment variable must be set for recovery")
}

sample_func <- select_ll(test_model, sample=TRUE)
if (is.null(sample_func)) {
  stop("DCE_REC_MODEL env variable must be exposed in mcce package rll_funcs vector")
}

if (! (experiment %in% c("NumericVDCE", "SymbolicVDCE", "PrefDCE"))) {
  stop("System Environment Variable DCE_EST_EXP not defined or unknown value")
}

# Experiment specific details/checks
if (experiment == "SymbolicVDCE") {
  if (! (displaytype %in% c("Absent", "Greyed"))) {
    stop("System Environment variable VDCE_DISPLAY should be defined")
  }
  filename <- paste(experiment, test_model, displaytype, jobid, tag, sep = "_")
} else {
  filename <- paste(experiment, test_model, jobid, tag, sep = "_")
}

# Get output filename and input data
outfile <- here::here("data", "output", paste0(filename, ".RData"))
datafile <- here::here("data", "output", paste0(filename, "_data.RDS"))

model_data <- get_estimation_data(model_file)
subjects <- unique(model_data$subject)

if (recovery_data == "") {
  print("Creating test dataset")
  test_data <- lapply(subjects, FUN = function(subjectid) {
    s_idx <- match(subjectid, subjects)
    subject_data <- model_data %>% filter(subject == subjectid)
    pars <- medians[s_idx, ]

    rmodel_wrapper(pars, subject_data, sample_func)
  })

  test_data <- test_data %>%
    bind_rows()

  saveRDS(test_data, file = datafile)
} else {
  test_data <- readRDS(file = here::here("data", "output", recovery_data))
}

# Model specification - should be identical to model estimation file
# < 0.3 participants were penalised, max trial length was 4.5 seconds
min_rt <- 0
max_rt <- 4.5
p_contam <- 0.02

acc_rej_drift <- c("v_acc_p", "v_acc_r", "v_rej_p", "v_rej_r")
stim_levels <- c("H", "L", "D")

parameters <- c(
  # alpha (dirichlet mixture pars) for each likelihood function exposed in mcce
  apply(expand.grid("alpha", names_ll()), 1, paste, collapse = "_"),
  # A - start point variability (sampled from U(0, A) where U is uniform dist)
  "A",
  # b_acc - threshold to accept based on evidence accumulaton in channel
  "b_acc",
  # b_rej - threshold to reject based on evidence accumulation in channel
  "b_rej",
  # t0 - residual time, bounded above by min response time for participant k
  "t0",
  # Drift rates rto accept/reject for different stimulus levels/attributes
  apply(expand.grid(acc_rej_drift, stim_levels), 1, paste, collapse = "_")
)

# Mixture counts should always come first
mix_counts <- 1:sum(startsWith(parameters, "alpha"))

priors <- list(
  theta_mu_mean = rep(0, length(parameters)),
  theta_mu_var = diag(rep(1, length(parameters)))
)
# Set alpha values to be mu 1, sigma 2
priors$theta_mu_mean[mix_counts] <- 1
diag(priors$theta_mu_var)[mix_counts] <- 2

# Create the Particle Metropolis within Gibbs sampler object ------------------

sampler <- pmwgs(
  data = test_data,
  pars = parameters,
  ll_func = dirichlet_mix_ll,
  prior = priors
)

sampler <- init(sampler)

burned <- run_stage(sampler, stage = "burn", iter = 5000, particles = 500, n_cores = cores)

save.image(outfile)

adapted <- run_stage(burned, stage = "adapt", iter = 10000, particles = 500, n_cores = cores, n_unique = 40)

save.image(outfile)

sampled <- run_stage(adapted, stage = "sample", iter = 10000, particles = 100, n_cores = cores, pdist_update_n = NA)

save.image(outfile)
