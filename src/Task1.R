require(pmwg)
require(rtdists)
library(MCMCpack)
require(tidyverse)
source("read_expyriment.R")
source("loglike.R")

#Get output filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  jobid <- Sys.getenv()["PBS_JOBID"]
  if (is.na(jobid)) {
    args[1] <- tempfile(pattern = "Task1_", tmpdir = ".", fileext = ".RData")
  } else {
    args[1] <- paste0("Task1_", jobid, ".RData")
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
task1_data <- read.expyriment.data(here::here("data/input/Task1"), "S*") %>%
  mutate(Correct = as.logical(Correct)) %>%
  mutate(price_match = PriceSalience %in% c("High", "Low"),
         rating_match = RatingSalience %in% c("High", "Low")) %>%
  mutate(trial_cat = case_when(price_match & rating_match ~ "both",
                               price_match & !rating_match ~ "psing",
                               rating_match & !price_match ~ "rsing",
                               !(price_match | rating_match) ~ "neither")) %>%
  group_by(subject_id, trial_cat) %>%
  mutate(pc_correct = mean(Correct)) %>%
  group_by(subject_id) %>%
  filter(min(pc_correct) >= 0.8) %>%
  rename(PriceRatingOrder = `b'PriceRatingOrder' `,
         ResponseCounterbalancing = `b'ResponseCounterbalancing' `,
         AcceptRejectFocus = `b'AcceptRejectFocus' `) %>%
  mutate(PriceRatingOrder = fct_relabel(PriceRatingOrder, cfix),
         ResponseCounterbalancing = fct_relabel(ResponseCounterbalancing, cfix),
         AcceptRejectFocus = fct_relabel(AcceptRejectFocus, cfix)) %>%
  mutate(PriceSalience = fct_recode(PriceSalience, !!!short_codes),
         RatingSalience = fct_recode(RatingSalience, !!!short_codes)) %>%
  drop_na()

# Two char labels for each cell of design,
# first char is price, second is quality, H=High, L=Low, D=Distractor
stim_levels <- c("HH", "HL", "HD", "LH", "LL", "LD", "DH", "DL", "DD")
# Clean no responses
accept_trials <- task1_data %>% filter(AcceptRejectFocus == "Accept")
mod_data <- data.frame(
  rt = task1_data$RT / 1000,
  subject = task1_data$subject_id,
  response = as.numeric(task1_data$Correct) + 1,
  cell = factor(
    paste0(
      as.character(task1_data$PriceSalience),
      as.character(task1_data$RatingSalience)
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
#Mixture counts should always come first
mix_counts <- 1:sum(startsWith(parameters, "alpha"))

# Specify likelihoods of parameters for individual models ---------------------

# Independant parallel self terminating
ll_IST <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"]  #nolint
  t0 <- x["t0"]

  yes <- 2 *
    dlba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    (1 - plba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1)) *
    (1 - plba_norm(ydat$rt, A, x["b_neg"], t0, x[ydat$v_neg], 1)**2)
  no <- 2 *
    dlba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    plba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    (1 - plba_norm(ndat$rt, A, x["b_pos"], t0, x[ndat$v_pos], 1))**2
  c(yes, no)
}

ll_IEX <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"]  #nolint
  t0 <- x["t0"]
  yes <- 2 *
    dlba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    plba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    (1 - plba_norm(ydat$rt, A, x["b_neg"], t0, x[ydat$v_neg], 1))**2
  no <- 2 *
    dlba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    (1 - plba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1)) *
    (1 - plba_norm(ndat$rt, A, x["b_pos"], t0, x[ndat$v_pos], 1)**2)
  c(yes, no)
}

ll_CYST <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"]  #nolint
  t0 <- x["t0"]

  yes <- dlba_norm(ydat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ydat$v_pos], sqrt(2)) *  #nolint
    (1 - plba_norm(ydat$rt, A, x["b_neg"], t0, x[ydat$v_neg], 1)**2)
  no <- 2 *
    dlba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    plba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    (1 - plba_norm(ndat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ndat$v_pos], sqrt(2)))  #nolint
  c(yes, no)
}

ll_CYEX <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"]  #nolint
  t0 <- x["t0"]

  yes <- dlba_norm(ydat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ydat$v_pos], sqrt(2)) *  #nolint
    (1 - plba_norm(ydat$rt, A, x["b_neg"], t0, x[ydat$v_neg], 1))**2
  no <- 2 *
    dlba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1) *
    (1 - plba_norm(ndat$rt, A, x["b_neg"], t0, x[ndat$v_neg], 1)) *
    (1 - plba_norm(ndat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ndat$v_pos], sqrt(2)))  #nolint
  c(yes, no)
}

ll_CNST <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"]  #nolint
  t0 <- x["t0"]

  yes <- 2 *
    dlba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    (1 - plba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1)) *
    (1 - plba_norm(ydat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ydat$v_neg], sqrt(2)))  #nolint
  no <- dlba_norm(ndat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ndat$v_neg], sqrt(2)) *  #nolint
    (1 - plba_norm(ndat$rt, A, x["b_pos"], t0, x[ndat$v_pos], 1))**2
  c(yes, no)
}

ll_CNEX <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"]  #nolint
  t0 <- x["t0"]

  yes <- 2 *
    dlba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    plba_norm(ydat$rt, A, x["b_pos"], t0, x[ydat$v_pos], 1) *
    (1 - plba_norm(ydat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ydat$v_neg], sqrt(2)))  #nolint
  no <- dlba_norm(ndat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ndat$v_neg], sqrt(2)) *  #nolint
    (1 - plba_norm(ndat$rt, A, x["b_pos"], t0, x[ndat$v_pos], 1)**2)
  c(yes, no)
}

ll_CB <- function(x, data) { # nolint
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"]  #nolint
  t0 <- x["t0"]

  yes <- dlba_norm(ydat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ydat$v_pos], sqrt(2)) *  #nolint
    (1 - plba_norm(ydat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ydat$v_neg], sqrt(2)))  #nolint
  no <- dlba_norm(ndat$rt, 2 * A, 2 * x["b_neg"], t0, 2 * x[ndat$v_neg], sqrt(2)) *  #nolint
    (1 - plba_norm(ndat$rt, 2 * A, 2 * x["b_pos"], t0, 2 * x[ndat$v_pos], sqrt(2)))  #nolint
  c(yes, no)
}

ll_funcs <- c(ll_IST, ll_IEX, ll_CYST, ll_CYEX, ll_CNST, ll_CNEX, ll_CB)

# Specify the log likelihood function -----------------------------------------
# For one accumulator, following terminology in Cox and Criss (2019)
# Start point variability is uniform 0 -> α
# Drift sampled from normal distribution with mean v and std dev s
# Threshold is ω
# Residual time is R (time to detect and start processing, and respond)
# For rtdists, we have A == α, b == ω, t0 == R, mean_v == v, sd_v == s
dirichlet_mix_ll <- function(x, data) {
  x <- exp(x)

  #Enforce alphas to be greater than 0.01 and less than 100
  if (any(x[mix_counts] < 0.01) || any(x[mix_counts] > 100)) {
    return(-1e10)
  }

  # Enforces b cannot be less than A, b parameter is thus threshold - A
  x["b_pos"] <- x["b_pos"] + x["A"]
  x["b_neg"] <- x["b_neg"] + x["A"]

  # all decision rules
  rdev <- rdirichlet(1, x[mix_counts])
  func_idx <- sample(mix_counts, 1, prob = rdev)
  ll_func <- ll_funcs[[func_idx]]
  trial_ll <- ll_func(x, data)
  new_like <- (1 - p_contam) * trial_ll +
    p_contam * (dunif(data$rt, min_rt, max_rt) / 2)
  sum(log(pmax(new_like, 1e-10)))
}

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

burned <- run_stage(sampler, stage = "burn", iter = 2000, particles = 500, n_cores = 36)

save.image(outfile)

adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 500)

save.image(outfile)

sampled <- run_stage(adapted, stage = "sample", iter = 5000, particles = 100)

save.image(outfile)
