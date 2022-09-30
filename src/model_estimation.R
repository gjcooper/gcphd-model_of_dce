require(pmwg)
require(rtdists)
library(dplyr)
library(MCMCpack)
library(stringi)
library(mcce)
library(purrr)

print(sessionInfo())

# Example Environment Variables for illustration and for debugging:
# For testing
# Sys.setenv(MCCE_EST_EXP="PrefDCE",
#            NCPUS=1,
#            MCCE_EXP_DATA="Pref_preprocessed.RDS",
#            MCCE_MIN_RT=0.35,
#            MCCE_MAX_RT=10,
#            MCCE_CONTAM=0.02,
#            MCCE_TAG="test_new",
#            RANDOM_SEED=101,
#            MCCE_MODEL="reduced",
#            MCCE_METHOD="model")
# For continuation
# Sys.setenv(MCCE_EST_EXP="PrefDCE",
#            NCPUS=3,
#            MCCE_ORIG_JOB_DATA="PrefDCE_S6I7q4nycmRv_short_burn_cont.RData",
#            MCCE_TAG="TestContinue",
#            MCCE_STAGES="sample")
# For recovery
# Sys.setenv(MCCE_EST_EXP="PrefDCE",
#            NCPUS=3,
#            MCCE_REC_MODEL="CB",
#            MCCE_REC_MED="median_alpha_pref.RDS",
#            MCCE_MODEL_FILE="PrefDCE_2506730.rcgbcm_Estimation5Model.RData",
#            MCCE_MIN_RT=0.35,
#            MCCE_MAX_RT=10,
#            MCCE_CONTAM=0.02,
#            RANDOM_SEED=3,
#            MCCE_MODEL="std",
#            MCCE_METHOD="recovery")
# Get environment variables to normal vars
# < 0.3 participants were penalised, max trial length was 4.5 seconds

# Experiment specific details/checks
vars <- check_env()
filename <- get_base_filename(vars)

# Get output filename and input data
outfile <- here::here("data", "output", paste0(filename, ".RData"))
if (vars$method %in% c("test", "model")) {
  model_data <- readRDS(here::here("data", "output", vars$experimental_data))
} else if (vars$method == "test") {
  test_data_files <- c(std = "PrefDCE_test.RData",
                       reduced = "PrefDCE_testred.RData",
                       reduced_more = "PrefDCE_testredmore.RData")
  testfile <- here::here("data", "output", test_data_files[vars$model])
} else if (vars$method == "continue") {
  # Pull in vars etc from earlier run
  load(file = here::here("data", "output", vars$early_data), envir = e <- new.env())
  for (varname in ls(e)) {
    if (!varname %in% c("outfile", "cores", "stages_to_run")) {
      print(paste("Assigning", varname))
      assign(varname, e[[varname]], envir = .GlobalEnv)
    }
  }
  vars$stages_to_run <- strsplit(stages_to_run, ",")[[1]]
} else if (vars$method == "recovery") {
  datafile <- here::here("data", "output", paste0(filename, "_data.RDS"))
  model_data <- get_estimation_data(vars$model_file)
  if (vars$recovery_data == "") {
    print("Creating test dataset")
    test_data <- generate_data(model_data, vars)
    saveRDS(test_data, file = datafile)
  } else {
    test_data <- readRDS(file = here::here("data", "output", vars$recovery_data))
  }
}

#Get the parameters from a separate file
source(paste0(vars$model, "_pars.R"))

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
dirichlet_func <- partial(dirichlet_mix_ll,
                          contaminant_prob = vars$p_contam,
                          alpha_indices = mix_counts,
                          min_rt = vars$min_rt,
                          max_rt = vars$max_rt,
                          tforms = vars$model)

source(paste0(vars$method, "_sampling.R"))
