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

test_model <- "IST"

sample_func <- rll_funcs[[match(test_model, ll_names)]]

# Get output filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  tag <- paste0("_untagged_recovery_", test_model)
} else {
  tag <- paste0("_", args[1], "_recovery_", test_model)
}

jobid <- Sys.getenv()["PBS_JOBID"]
if (is.na(jobid)) {
  filename <- tempfile(pattern = "Task1_", tmpdir = ".", fileext = tag)
} else {
  filename <- paste0("Task1_", jobid, tag)
}

outfile <- here::here("data", "output", paste0(filename, ".RData"))
datafile <- here::here("data", "output", paste0(filename, "_data.RDS"))

subjects <- unique(model_data$subject)

print("Creating test dataset")
test_data <- lapply(subjects, FUN = function(subjectid) {
  s_idx <- match(subjectid, subjects)
  subject_data <- model_data %>% filter(subject == subjectid)
  pars <- medians[s_idx, ]

  sample_func(pars, subject_data)
})

test_data <- test_data %>%
  bind_rows()

saveRDS(test_data, file = datafile)

# Model specification - should be identical to model estimation file
# Sum and difference of evidence rates for positive and negative accumulators
sum_diff <- c("v_pos", "v_neg")
stim_levels <- c("HH", "HL", "HD", "LH", "LL", "LD", "DH", "DL", "DD")

parameters <- c(
  # Parallel mixture counts
  "alpha_IST", "alpha_IEX",
  # Coactive mixture probabilities
  "alpha_CYST", "alpha_CYEX", "alpha_CNST", "alpha_CNEX", "alpha_CB",
  # A - start point variability (sampled from U(0, A) where U is uniform dist)
  "A",
  # b_pos - threshold for positive evidence accumulation
  "b_pos",
  # b_neg - threshold for negative evidence accumulation
  "b_neg",
  # t0 - residual time, bounded above by min response time for participant k
  "t0"
)
# beta and delta - 9 versions each for each of beta and delta,
# corresponding to the 9 cells of the experimental design
parameters <- c(parameters, apply(
  expand.grid(sum_diff, stim_levels), 1, paste,
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

start_mu <- c(0, 0, 0, 0, 0, 0, 0, .4, .2, .2, -2, rep(c(1.3, .3), 9))
start_sig <- MCMCpack::riwish(sampler$n_pars * 2, diag(sampler$n_pars))

min_rt <- 0
max_rt <- 4.5
p_contam <- 0.02

sampler <- init(sampler, start_mu = start_mu, start_sig = start_sig)

burned <- run_stage(sampler, stage = "burn", iter = 2000, particles = 500, n_cores = 26)

save.image(outfile)

adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 500, n_cores = 26)

save.image(outfile)

sampled <- run_stage(adapted, stage = "sample", iter = 5000, particles = 100, n_cores = 26)

save.image(outfile)
