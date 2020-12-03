require(pmwg)
require(rtdists)
library(MCMCpack)

# Get output filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  jobid <- Sys.getenv()["PBS_JOBID"]
  if (is.na(jobid)) {
    args[1] <- tempfile(pattern = "Task1_", tmpdir = ".", fileext = ".RData")
  } else {
    args[1] <- paste0("Task1_", jobid, ".RData")
  }
}
outfile <- here::here("data", "output", args[1])

task1_data <- readRDS(here::here("data", "output", "Task1_preprocessed.RDS"))

# Only accept trials
# Create simplifed data for modelling with rtdists
# add drift parameter names
mod_data <- task1_data %>%
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
  # Parallel mixture counts
  "alpha_IST", "alpha_IEX",
  # Coactive mixture probabilities
  "alpha_CYST", "alpha_CYEX", "alpha_CNST", "alpha_CNEX", "alpha_CB",
  # A - start point variability (sampled from U(0, A) where U is uniform dist)
  "A",
  # b_acc - threshold to accept based on evidence accumulaton in channel
  "b_acc",
  # b_rej - threshold to reject based on evidence accumulation in channel
  "b_neg",
  # t0 - residual time, bounded above by min response time for participant k
  "t0"
)
parameters <- c(parameters, apply(
  expand.grid(acc_rej_drift, stim_levels), 1, paste,
  collapse = "_"
))
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

burned <- run_stage(sampler, stage = "burn", iter = 2000, particles = 500, n_cores = 36)

save.image(outfile)

adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 500, n_cores = 36)

save.image(outfile)

sampled <- run_stage(adapted, stage = "sample", iter = 5000, particles = 100, n_cores = 36)

save.image(outfile)
