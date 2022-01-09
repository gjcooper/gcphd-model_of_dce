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
                "DCE_REC_MODEL", "DCE_ORIG_JOB_DATA")

envars <- Sys.getenv(known_vars)
print(envars)
envars <- as.list(envars)
experiment <- envars$DCE_EST_EXP
displaytype <- envars$VDCE_DISPLAY
new_cores <- ifelse(envars$NCPUS == "", 1, as.numeric(envars$NCPUS))
jobid <- ifelse(
  envars$PBS_JOBID == "",
  stri_rand_strings(1, 12),
  envars$PBS_JOBID
)
tag <- ifelse(envars$VDCE_TAG == "", "untagged", envars$VDCE_TAG)
test_model <- envars$DCE_REC_MODEL
early_data <- envars$DCE_ORIG_JOB_DATA


if (test_model == "") {
  stop("DCE_REC_MODEL environment variable must be set for recovery")
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
  filename <- paste(experiment, test_model, displaytype, jobid, tag, sep = "_")
} else {
  filename <- paste(experiment, test_model, jobid, tag, sep = "_")
}

# Get output filename and input data
new_outfile <- here::here("data", "output", paste0(filename, ".RData"))

# Pull in vars etc from earlier run
load(file = here::here("data", "output", early_data))

outfile <- new_outfile
cores <- new_cores

if (!exists("burned")) {
	burned <- run_stage(sampler, stage = "burn", iter = 5000, particles = 500, n_cores = cores)
}

save.image(outfile)

if (!exists("adapted")) {
	adapted <- run_stage(burned, stage = "adapt", iter = 10000, particles = 500, n_cores = cores, n_unique = 40)
}

save.image(outfile)

if (!exists("sampled")) {
	sampled <- run_stage(adapted, stage = "sample", iter = 10000, particles = 100, n_cores = cores, pdist_update_n = NA)
}

save.image(outfile)
