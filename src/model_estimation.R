require(pmwg)
require(rtdists)
library(dplyr)
library(MCMCpack)
library(stringi)
library(mcce)
library(purrr)

print(sessionInfo())

# For debugging:
# Sys.setenv(MCCE_EST_EXP="SymbolicVDCE", VDCE_DISPLAY="Absent", NCPUS=3)
# Sys.setenv(MCCE_EST_EXP="NumericVDCE", NCPUS=3, MCCE_EXP_DATA="Task1_preprocessed.RDS", MCCE_MIN_RT=0, MCCE_MAX_RT=4.5, MCCE_CONTAM=0.02)
# Sys.setenv(MCCE_EST_EXP="PrefDCE", NCPUS=3, MCCE_EXP_DATA="Pref_preprocessed.RDS", MCCE_MIN_RT=0.35, MCCE_MAX_RT=10, MCCE_CONTAM=0.02)
# For testing
# Sys.setenv(MCCE_EST_EXP="PrefDCE", NCPUS=1, MCCE_EXP_DATA="Pref_preprocessed.RDS", MCCE_MIN_RT=0.35, MCCE_MAX_RT=10, MCCE_CONTAM=0.02, MCCE_TAG="test_new", RANDOM_SEED=101)
# For continuation
# Sys.setenv(MCCE_EST_EXP="PrefDCE", NCPUS=3, MCCE_ORIG_JOB_DATA="PrefDCE_S6I7q4nycmRv_short_burn_cont.RData", MCCE_TAG="TestContinue", MCCE_STAGES="sample")
# Get environment variables to normal vars
known_vars <- c("MCCE_EST_EXP", "VDCE_DISPLAY", "NCPUS", "PBS_JOBID", "MCCE_TAG",
                "MCCE_EXP_DATA", "MCCE_MIN_RT", "MCCE_MAX_RT", "MCCE_CONTAM",
                "RANDOM_SEED", "MCCE_MODEL", "MCCE_METHOD",
                "MCCE_ORIG_JOB_DATA", "MCCE_STAGES")

envars <- Sys.getenv(known_vars)
envars <- as.list(envars)
print(envars %>% data.frame %>% t)
experiment <- envars$MCCE_EST_EXP
# < 0.3 participants were penalised, max trial length was 4.5 seconds
min_rt <- as.numeric(envars$MCCE_MIN_RT)
max_rt <- as.numeric(envars$MCCE_MAX_RT)
p_contam <- as.numeric(envars$MCCE_CONTAM)
displaytype <- envars$VDCE_DISPLAY
cores <- ifelse(envars$NCPUS == "", 1, as.numeric(envars$NCPUS))
jobid <- ifelse(
  envars$PBS_JOBID == "",
  stri_rand_strings(1, 12),
  envars$PBS_JOBID
)
tag <- ifelse(envars$MCCE_TAG == "", "untagged", envars$MCCE_TAG)
early_data <- envars$MCCE_ORIG_JOB_DATA
stages_to_run <- envars$MCCE_STAGES


experimental_data <- envars$MCCE_EXP_DATA
if (envars$RANDOM_SEED != "") {
  set.seed(as.numeric(envars$RANDOM_SEED))
}
model <- envars$MCCE_MODEL
method <- envars$MCCE_METHOD

#Tests
if (! (experiment %in% c("NumericVDCE", "SymbolicVDCE", "PrefDCE"))) {
  stop("System Environment Variable MCCE_EST_EXP not defined or unknown value")
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
if (is.na(min_rt) || is.na(max_rt) || is.na(p_contam)) {
  stop("All of MCCE_MIN_RT, MCCE_MAX_RT, MCCE_CONTAM must be provided and numeric")
}
if (!between(p_contam, 0, 1)) {
  stop("MCCE_CONTAM should be a value between 0 and 1")
}

# Model tests
if (! (model %in% c("std", "reduced", "reduced_more"))) {
  stop("System Environment Variable MCCE_MODEL not defined or unknown value")
}

# Method tests
if (! (method %in% c("test", "model", "continue"))) {
  stop("System Environment Variable MCCE_METHOD not defined or unknown value")
}

# Get output filename and input data
outfile <- here::here("data", "output", paste0(filename, ".RData"))
model_data <- readRDS(here::here("data", "output", experimental_data))
if (method == "test") {
  test_data_files <- c(std = "PrefDCE_test.RData",
                       reduced = "PrefDCE_testred.RData",
                       reduced_more = "PrefDCE_testredmore.RData")
  testfile <- here::here("data", "output", test_data_files[model])
} else if (method == "continue") {
  if (early_data == "") {
    stop("MCCE_ORIG_JOB_DATA environment variable must be set to continue")
  }
  if (stages_to_run == "") {
    stop("MCCE_STAGES envvar must be a comma separated list of PMwG stages (burn,adapt,sample)")
  }
  # Pull in vars etc from earlier run
  load(file = here::here("data", "output", early_data), envir = e <- new.env())

  for (varname in ls(e)) {
    if (!varname %in% c("outfile", "cores", "stages_to_run")) {
      print(paste('Assigning', varname))
      assign(varname, e[[varname]], envir = .GlobalEnv)
    }
  }
  outfile <- new_outfile
  cores <- new_cores
  stages_to_run <- strsplit(stages_to_run, ",")[[1]]
}



#Get the parameters from a separate file
source(paste0(model, "_pars.R"))

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
dirichlet_func <- partial(dirichlet_mix_ll, contaminant_prob = p_contam, alpha_indices = mix_counts, min_rt = min_rt, max_rt = max_rt, tforms = model)


source(paste0(method, "_sampling.R"))
