require(pmwg)
require(rtdists)
library(dplyr)
library(MCMCpack)
devtools::load_all()

# For debugging:
# Sys.setenv(DCE_EST_EXP="NumericVDCE")
# Sys.setenv(DCE_EST_EXP="SymbolicVDCE", VDCE_DISPLAY="Absent")
# Get experiment to model
experiment <- Sys.getenv("DCE_EST_EXP")
if (! (experiment %in% c("NumericVDCE", "SymbolicVDCE"))) {
  stop("System Environment Variable DCE_EST_EXP not defined or unknown value")
}

if (experiment == "NumericVDCE") {
  fileprefix <- "NumericVDCE_"
  infile <- "Task1_preprocessed.RDS"
} else if (experiment == "SymbolicVDCE") {
  #Get display to analyse
  displaytype <- Sys.getenv("VDCE_DISPLAY")
  # If VDCE_DISPLAY not defined should error out
  if (! (displaytype %in% c("Absent", "Greyed"))) {
    stop("System Environment variable VDCE_DISPLAY should be defined")
  }
  fileprefix <- paste0("SymbolicVDCE_", displaytype, "_")
  infile <- paste0("Task2_preprocessed_", displaytype, ".RDS")
}

cores <- Sys.getenv("NCPUS")
if (cores == "") {
  cores <- 1
}

# Get output filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  tag <- "_untagged"
} else {
  tag <- paste0("_", args[1])
}

jobid <- Sys.getenv("PBS_JOBID")
if (jobid == "") {
  filename <- tempfile(pattern = fileprefix, tmpdir = ".", fileext = tag)
} else {
  filename <- paste0(fileprefix, jobid, tag)
}

outfile <- here::here("data", "output", paste0(filename, ".RData"))

cleaned_data <- readRDS(here::here("data", "output", infile))

# Only accept trials
# Create simplifed data for modelling with rtdists
# add drift parameter names
mod_data <- cleaned_data %>%
  filter(acceptAND) %>%
  transmute(
    rt = RT / 1000,
    subject = subject_id,
    accept = as.numeric(Accept) + 1,
    price = Price,
    rating = Rating) %>%
  mutate(
    v_acc_p = paste0("v_acc_p_", price),
    v_rej_p = paste0("v_rej_p_", price),
    v_acc_r = paste0("v_acc_r_", rating),
    v_rej_r = paste0("v_rej_r_", rating)
  )

# < 0.3 participants were penalised, max trial length was 4.5 seconds
min_rt <- 0
max_rt <- 4.5
p_contam <- 0.02

acc_rej_drift <- c("v_acc_p", "v_acc_r", "v_rej_p", "v_rej_r")
stim_levels <- c("H", "L", "D")

parameters <- c(
  # alpha (dirichlet mixture pars) for each likelihood function exposed in mcce
  apply(expand.grid("alpha", names(ll_funcs)), 1, paste, collapse = "_"),
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

# Likelihood functions from separate file

priors <- list(
  theta_mu_mean = rep(0, length(parameters)),
  theta_mu_var = diag(rep(1, length(parameters)))
)
# Set alpha values to be mu 1, sigma 2
priors$theta_mu_mean[mix_counts] <- 1
diag(priors$theta_mu_var)[mix_counts] <- 2

# Create the Particle Metropolis within Gibbs sampler object ------------------

sampler <- pmwgs(
  data = mod_data,
  pars = parameters,
  ll_func = dirichlet_mix_ll,
  prior = priors
)

sampler <- init(sampler)

burned <- run_stage(sampler, stage = "burn", iter = 2000, particles = 500, n_cores = cores)

save.image(outfile)

adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 500, n_cores = cores)

save.image(outfile)

sampled <- run_stage(adapted, stage = "sample", iter = 5000, particles = 100, n_cores = cores)

save.image(outfile)
