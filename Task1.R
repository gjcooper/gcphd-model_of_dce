require(psamplers)
require(rtdists)
library(MCMCpack)
source("read_expyriment.R")

#Get output filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  jobid <- Sys.getenv()["PBS_JOBID"]
  if (is.na(jobid)) {
    args[1] <- tempfile(pattern = "cce_burn_", tmpdir = ".", fileext = ".RData")
  } else {
    args[1] <- paste0("cce_burn_", jobid, ".RData")
  }
}
outfile <- paste0("data/output/", args[1])

task1_data <- read.expyriment.data("data/input/Task1/", "S*")
# Two char labels for each cell of design,
# first char is price, second is quality, H=High, L=Low, D=Distractor
stim_levels <- c("HH", "HL", "HD", "LH", "LL", "LD", "DH", "DL", "DD")
# Clean no responses
task1_data <- task1_data[complete.cases(task1_data), ]
task1_data <- task1_data[task1_data$`b'AcceptRejectFocus' ` == "b'Accept'", ]
levels(task1_data$PriceSalience) <- c("H", "L", "D")
levels(task1_data$RatingSalience) <- c("H", "L", "D")
mod_data <- data.frame(
  rt = task1_data$RT / 1000,
  subject = task1_data$subject_id,
  response = as.numeric(task1_data$Correct),
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

decision_rules <- c("IST", "IEX", "CYST", "CYEX", "CNST", "CNEX", "CB")

parameters <- c(
  # Mixture counts for each decision rule
  apply(expand.grid("alpha", decision_rules), 1, paste, collapse="_"),
  # A - start point variability (sampled from U(0, A) where U is uniform dist)
  "A",
  # b_pos - threshold for positive evidence accumulation
  "b_pos",
  # b_neg - threshold for negative evidence accumulation
  "b_neg",
  # t0 - residual time, bounded above by min response time for participant k
  "t0",
  # Positive (accept) and negative (reject) drift rate by cell of design
  apply(expand.grid(c("v_pos", "v_neg"), stim_levels), 1, paste, collapse = "_")
)
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


# Specify the model and surround it with log like wrapper ---------------------
# Based on a given string (the model ID) then select a specific model to run
# for all rows of data, for all subjects, for all iterations. (No model
# selection within the process)
selection_ll <- function(decision_rule) {
  ll_func <- ll_funcs[[match(decision_rule, decision_rules)]]
  wrapped_ll <- function(x, data) {
    x["b_pos"] <- x["b_pos"] + x["A"]
    x["b_neg"] <- x["b_neg"] + x["A"]

    # all decision rules
    trial_ll <- ll_func(x, data)
    new_like <- (1 - p_contam) * trial_ll +
      p_contam * (dunif(data$rt, min_rt, max_rt) / 2)
    sum(log(pmax(new_like, 1e-10)))
  }
}

priors <- list(
  theta_mu = rep(0, length(parameters)),
  theta_sig = diag(rep(1, length(parameters)))
)
# Set alpha values to be mu 1, sigma 2
priors$theta_mu[mix_counts] <- 1
diag(priors$theta_sig)[mix_counts] <- 2

# Create the Particle Metropolis within Gibbs sampler object ------------------

sampler <- pmwgs(
  data = mod_data,
  pars = parameters,
  ll_func = dirichlet_mix_ll,
  prior = priors
)

start_mu <- c(0, 0, 0, 0, 0, 0, 0, .4, .2, .2, -2, rep(c(1.3, .3), 9))
start_sig <- MCMCpack::riwish(sampler$n_pars * 2, diag(sampler$n_pars))

sampler <- init(sampler, theta_mu = start_mu, theta_sig = start_sig)

burned <- run_stage(sampler, stage = "burn", iter = 200, particles = 200)

adapted <- run_stage(burned, stage = "adapt", iter = 200, particles = 200)

sampled <- run_stage(adapted, stage = "sample", iter = 100, particles = 100)

save.image(outfile)
