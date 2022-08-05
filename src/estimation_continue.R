require(pmwg)
require(rtdists)
library(dplyr)
library(MCMCpack)
library(stringi)
library(mcce)

print(sessionInfo())

# For debugging:
# Sys.setenv(DCE_EST_EXP="PrefDCE", NCPUS=3, DCE_REC_MODEL="CB", DCE_ORIG_JOB_DATA="PrefDCE_2506730.rcgbcm_Estimation5Model.RData")
# Get environment variables to normal vars
known_vars <- c("DCE_EST_EXP", "VDCE_DISPLAY", "NCPUS", "PBS_JOBID", "VDCE_TAG",
                "DCE_EXP_DATA", "DCE_MIN_RT", "DCE_MAX_RT", "DCE_CONTAM",
                "RANDOM_SEED", "DCE_ORIG_JOB_DATA")

envars <- Sys.getenv(known_vars)
envars <- as.list(envars)
print(envars %>% data.frame %>% t)
experiment <- envars$DCE_EST_EXP
# < 0.3 participants were penalised, max trial length was 4.5 seconds
min_rt <- as.numeric(envars$DCE_MIN_RT)
max_rt <- as.numeric(envars$DCE_MAX_RT)
p_contam <- as.numeric(envars$DCE_CONTAM)
displaytype <- envars$VDCE_DISPLAY
new_cores <- ifelse(envars$NCPUS == "", 1, as.numeric(envars$NCPUS))
jobid <- ifelse(
  envars$PBS_JOBID == "",
  stri_rand_strings(1, 12),
  envars$PBS_JOBID
)
tag <- ifelse(envars$VDCE_TAG == "", "untagged", envars$VDCE_TAG)
early_data <- envars$DCE_ORIG_JOB_DATA

if (envars$RANDOM_SEED != "") {
  set.seed(as.numeric(envars$RANDOM_SEED))
}

if (early_data == "") {
  stop("DCE_ORIG_JOB_DATA environment variable must be set to continue")
}

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
new_outfile <- here::here("data", "output", paste0(filename, ".RData"))

# Pull in vars etc from earlier run
load(file = here::here("data", "output", early_data))

outfile <- new_outfile
cores <- new_cores

if (!("burn" %in% sampler$samples$stage)) {
	sampler <- run_stage(sampler, stage = "burn", iter = 5000, particles = 500, n_cores = cores)
}

save.image(outfile)

if (!("adapt" %in% sampler$samples$stage)) {
	sampler <- run_stage(sampler, stage = "adapt", iter = 10000, particles = 500, n_cores = cores, n_unique = 40)
}

save.image(outfile)

if (!("sample" %in% sampler$samples$stage)) {
	sampler <- run_stage(sampler, stage = "sample", iter = 10000, particles = 100, n_cores = cores, pdist_update_n = NA)
}

save.image(outfile)
