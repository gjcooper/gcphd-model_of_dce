require(pmwg)
require(rtdists)
library(MCMCpack)
require(tidyverse)
source(here::here("src", "read_expyriment.R"))
source(here::here("src", "loglike.R"))

#Get display to analyse
displaytype <- Sys.getenv("VDCE_DISPLAY")
# If VDCE_DISPLAY not defined should error out
if (! (displaytype %in% c("Absent", "Greyed"))) {
  stop("System Environment variable VDCE_DISPLAY should be defined")
}
#Get output filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  jobid <- Sys.getenv()["PBS_JOBID"]
  if (is.na(jobid)) {
    args[1] <- tempfile(pattern = paste0("Task2_", displaytype, "_"), tmpdir = ".", fileext = ".RData")
  } else {
    args[1] <- paste0("Task2_", displaytype, "_", jobid, ".RData")
  }
}
outfile <- here::here("data/output", args[1])

cfix <- function(x) {
  substr(x, 3, nchar(x) - 1)
}

short_codes <- c(H = "High", L = "Low", D = "OutOfBounds")

# Read in data, turn Correct column into logical
# Get double target, single target and double distractor trials
# Clean based on minimum % correct in all of four trial categories,
# Clean the column names and drop rows with no response (NA values)
# Clean values in the renamed columns too
task2_data <- read.expyriment.data(here::here("data/input/Task2"), "S*") %>%
  mutate(Correct = as.logical(Correct)) %>%
  filter(Display == displaytype) %>%
  mutate(
    price_match = PriceSalience %in% c("High", "Low"),
    rating_match = RatingSalience %in% c("High", "Low")
  ) %>%
  mutate(trial_cat = case_when(
    price_match & rating_match ~ "both",
    price_match & !rating_match ~ "psing",
    rating_match & !price_match ~ "rsing",
    !(price_match | rating_match) ~ "neither"
  )) %>%
  group_by(subject_id, trial_cat) %>%
  mutate(pc_correct = mean(Correct)) %>%
  group_by(subject_id) %>%
  filter(min(pc_correct) >= 0.8) %>%
  rename(
    PriceRatingOrder = `b'PriceRatingOrder' `,
    GreyedItemDisplay = `b'GreyedItemDisplay' `,
    AcceptRejectFocus = `b'AcceptRejectFocus' `
  ) %>%
  mutate(
    PriceRatingOrder = fct_relabel(PriceRatingOrder, cfix),
    GreyedItemDisplay = fct_relabel(GreyedItemDisplay, cfix),
    AcceptRejectFocus = fct_relabel(AcceptRejectFocus, cfix)
  ) %>%
  mutate(
    PriceSalience = fct_recode(PriceSalience, !!!short_codes),
    RatingSalience = fct_recode(RatingSalience, !!!short_codes)
  ) %>%
  drop_na()

# Two char labels for each cell of design,
# first char is price, second is quality, H=High, L=Low, D=Distractor
stim_levels <- c("HH", "HL", "HD", "LH", "LL", "LD", "DH", "DL", "DD")
accept_trials <- task2_data %>% filter(AcceptRejectFocus == "Accept")
mod_data <- data.frame(
  rt = task2_data$RT / 1000,
  subject = task2_data$subject_id,
  response = as.numeric(task2_data$Correct) + 1,
  cell = factor(
    paste0(
      as.character(task2_data$PriceSalience),
      as.character(task2_data$RatingSalience)
    ),
    labels = stim_levels
  )
)
mod_data$v_pos <- paste0("v_pos_", mod_data$cell)
mod_data$v_neg <- paste0("v_neg_", mod_data$cell)

#< 0.3 participants were penalised, max trial length was 4.5 seconds
min_rt <- 0
max_rt <- 4.5
p_contam <- 0.02

# Sum and difference of evidence rates for positive and negative accumulators
sum_diff <- c("v_pos", "v_neg")

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

start_mu <- c(0, 0, 0, 0, 0, 0, 0, .4, .2, .2, -2, rep(c(1.3, .3), 9))
start_sig <- MCMCpack::riwish(sampler$n_pars * 2, diag(sampler$n_pars))

sampler <- init(sampler, start_mu = start_mu, start_sig = start_sig)

burned <- run_stage(sampler, stage = "burn", iter = 2000, particles = 500, n_cores = 34)

save.image(outfile)

adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 500, n_cores = 34)

save.image(outfile)

sampled <- run_stage(adapted, stage = "sample", iter = 5000, particles = 100, n_cores = 34)

save.image(outfile)
