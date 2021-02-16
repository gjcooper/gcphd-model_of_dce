require(pmwg)
require(rtdists)
library(dplyr)
library(MCMCpack)
devtools::load_all()

outdir <- here::here("data", "output")

get_estimation_data <- function(filename) {
  load(file.path(outdir, filename), ex <- new.env())
  ex$sampled$data %>%
    tibble()
}

medians <- readRDS(file.path(outdir, "median_alpha_exp1.RDS"))

test_model <- Sys.getenv("DCE_REC_MODEL")
if (test_model == "") {
  stop("DCE_REC_MODEL environment variable must be set for recovery")
}

sample_func <- ll_funcs[[match(test_model, names(ll_funcs))]]$sample
if (is.null(sample_func)) {
  stop("DCE_REC_MODEL env variable must be exposed in mcce package rll_funcs vector")
}

# Get output filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  tag <- paste0("_untagged_recovery_", test_model)
} else {
  tag <- paste0("_", args[1], "_recovery_", test_model)
}

jobid <- Sys.getenv("PBS_JOBID")
if (jobid == "") {
  filename <- tempfile(pattern = "Task1_", tmpdir = ".", fileext = tag)
} else {
  filename <- paste0("Task1_", jobid, tag)
}

outfile <- file.path(outdir, paste0(filename, ".RData"))
datafile <- file.path(outdir, paste0(filename, "_data.RDS"))

model_data <- get_estimation_data("Task1_1069902.rcgbcm_CorrectedTry1.RData")
subjects <- unique(model_data$subject)

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

# Model specification - should be identical to model estimation file
# < 0.3 participants were penalised, max trial length was 4.5 seconds
min_rt <- 0
max_rt <- 4.5
p_contam <- 0.02

acc_rej_drift <- c("v_acc_p", "v_acc_r", "v_rej_p", "v_rej_r")
stim_levels <- c("H", "L", "D")

parameters <- c(
  # Parallel mixture counts
  "alpha_IST", "alpha_IEX",
  # Coactive mixture probabilities
  "alpha_CB",
  # A - start point variability (sampled from U(0, A) where U is uniform dist)
  "A",
  # b_acc - threshold to accept based on evidence accumulaton in channel
  "b_acc",
  # b_rej - threshold to reject based on evidence accumulation in channel
  "b_rej",
  # t0 - residual time, bounded above by min response time for participant k
  "t0"
)
parameters <- c(parameters, apply(
  expand.grid(acc_rej_drift, stim_levels), 1, paste,
  collapse = "_"
))
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

burned <- run_stage(sampler, stage = "burn", iter = 2000, particles = 500, n_cores = 26)

save.image(outfile)

adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 500, n_cores = 26)

save.image(outfile)

sampled <- run_stage(adapted, stage = "sample", iter = 5000, particles = 100, n_cores = 26)

save.image(outfile)
