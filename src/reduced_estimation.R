require(pmwg)
require(rtdists)
library(dplyr)
library(MCMCpack)
library(stringi)
library(mcce)
library(purrr)

print(sessionInfo())

# For debugging:
# Sys.setenv(DCE_EST_EXP="SymbolicVDCE", VDCE_DISPLAY="Absent", NCPUS=3)
# Sys.setenv(DCE_EST_EXP="NumericVDCE", NCPUS=3, DCE_EXP_DATA="Task1_preprocessed.RDS", DCE_MIN_RT=0, DCE_MAX_RT=4.5, DCE_CONTAM=0.02)
# Sys.setenv(DCE_EST_EXP="PrefDCE", NCPUS=3, DCE_EXP_DATA="Pref_preprocessed.RDS", DCE_MIN_RT=0.35, DCE_MAX_RT=10, DCE_CONTAM=0.02)
# Get environment variables to normal vars
known_vars <- c("DCE_EST_EXP", "VDCE_DISPLAY", "NCPUS", "PBS_JOBID", "VDCE_TAG",
                "DCE_EXP_DATA", "DCE_MIN_RT", "DCE_MAX_RT", "DCE_CONTAM",
                "RANDOM_SEED")

envars <- Sys.getenv(known_vars)
envars <- as.list(envars)
print(envars %>% data.frame %>% t)
experiment <- envars$DCE_EST_EXP
# < 0.3 participants were penalised, max trial length was 4.5 seconds
min_rt <- as.numeric(envars$DCE_MIN_RT)
max_rt <- as.numeric(envars$DCE_MAX_RT)
p_contam <- as.numeric(envars$DCE_CONTAM)
displaytype <- envars$VDCE_DISPLAY
cores <- ifelse(envars$NCPUS == "", 1, as.numeric(envars$NCPUS))
jobid <- ifelse(
  envars$PBS_JOBID == "",
  stri_rand_strings(1, 12),
  envars$PBS_JOBID
)
tag <- ifelse(envars$VDCE_TAG == "", "untagged", envars$VDCE_TAG)
experimental_data <- envars$DCE_EXP_DATA
if (envars$RANDOM_SEED != "") {
  set.seed(as.numeric(envars$RANDOM_SEED))
}

#Tests
if (! (experiment %in% c("NumericVDCE", "SymbolicVDCE", "PrefDCE"))) {
  stop("System Environment Variable DCE_EST_EXP not defined or unknown value")
}

# Experiment specific details/checks
if (experiment == "SymbolicVDCE") {
  if (! (displaytype %in% c("Absent", "Greyed"))) {
    stop("System Environment variable VDCE_DISPLAY should be defined")
  }
  filename <- paste(experiment, displaytype, jobid, tag, sep = "_")
} else {
  filename <- paste(experiment, jobid, tag, sep = "_")
}

# Numeric env variable tests
if (is.na(min_rt) | is.na(max_rt) | is.na(p_contam)) {
  stop("All of DCE_MIN_RT, DCE_MAX_RT, DCE_CONTAM must be provided and numeric")
}
if (!between(p_contam, 0, 1)) {
  stop("DCE_CONTAM should be a value between 0 and 1")
}

# Get output filename and input data
outfile <- here::here("data", "output", paste0(filename, ".RData"))
model_data <- readRDS(here::here("data", "output", experimental_data))



acc_rej_drift <- c("v_acc_p", "v_acc_r")
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
  apply(expand.grid(acc_rej_drift, stim_levels), 1, paste, collapse = "_"),
  "beta0_p", "beta0_r", "beta1_p", "beta1_r"
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

dirichlet_func <- partial(dirichlet_reduced_mix, contaminant_prob = p_contam, alpha_indices = mix_counts, min_rt = min_rt, max_rt = max_rt)

sampler <- pmwgs(
  data = model_data,
  pars = parameters,
  ll_func = dirichlet_func,
  prior = priors
)

sampler <- init(sampler)

sampler <- run_stage(sampler, stage = "burn", iter = 200, particles = 500, epsilon = 1, n_cores = cores)
save.image(outfile)
sampler <- run_stage(sampler, stage = "burn", iter = 300, particles = 500, epsilon = 0.8, n_cores = cores)
save.image(outfile)
sampler <- run_stage(sampler, stage = "burn", iter = 500, particles = 500, epsilon = 0.6, n_cores = cores)
save.image(outfile)
sampler <- run_stage(sampler, stage = "burn", iter = 4000, particles = 500, n_cores = cores)
save.image(outfile)
