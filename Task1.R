require(psamplers)
require(rtdists)
source("read_expyriment.R")

task1_data <- read.expyriment.data("data/input/Task1/", "S*")
# Two char labels for each cell of design,
# first char is price, second is quality, H=High, L=Low, D=Distractor
stim_levels <- c("HH", "HL", "HD", "LH", "LL", "LD", "DH", "DL", "DD")
# Clean no responses
task1_data <- task1_data[complete.cases(task1_data), ]
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
mod_data$beta <- paste0("beta_", mod_data$cell)
mod_data$delta <- paste0("delta_", mod_data$cell)

# Sum and difference of evidence rates for positive and negative accumulators
sum_diff <- c("beta", "delta")

parameters <- c(
  # Parallel mixture probabilities
  "lambda_IST", "lambda_IEX",
  # Coactive mixture probabilities
  "lambda_CYST", "lambda_CYEX", "lambda_CNST", "lambda_CNEX",
  # "lambda_CB", -> set to 1 - sum(other lambdas)
  # ω - Sum of response thresholds
  "omega",
  # ν - how much threshold is subject to start point variability
  "nu",
  # w - proportion of total threshold given to negative accunmulator
  "w",
  # R - residual time, bounded above by min response time for participant k
  "R"
)
# beta and delta - 9 versions each for each of beta and delta,
# corresponding to the 9 cells of the experimental design
parameters <- c(parameters, apply(
  expand.grid(sum_diff, stim_levels), 1, paste,
  collapse = "_"
))

pars <- rnorm(length(parameters))
names(pars) <- parameters

# Helper methods --------------------------------------------------------------
alpha <- function(pars) pars["nu"] * pars["omega"]

omega <- function(pars, osign) {
  if (osign == "+") {
    return(
      pars["nu"] * pars["omega"] +
        (1 - pars["nu"]) * (1 - pars["w"]) * pars["omega"]
    )
  }
  # otherwise
  pars["nu"] * pars["omega"] + (1 - pars["nu"]) * pars["w"] * pars["omega"]
}

yes <- function(data) data[data$response == 2, ]
no <- function(data) data[data$response == 1, ]

drift <- function(pars, data, osign) {
  # drift is derived from β (sum) and δ (difference) of evidence rates for pos
  # and neg accumulators
  # β is sum(v+, v-) and δ is diff(v+, v-)
  # So v+ is (β + δ)/2 and v- is (β - δ)/2
  if (osign == "+") {
    return((pars[data$beta] + pars[data$delta]) / 2)
  }
  (pars[data$beta] - pars[data$delta]) / 2
}


# individual race pdf and cdfs ------------------------------------------------
pos_f <- function(x, data, drifts, scale = 1) {
  dlba_norm(
    data,
    A = scale * alpha(x),
    b = scale * omega(x, "+"),
    t0 = x["R"],
    mean_v = scale * drifts,
    sd_v = sqrt(scale)
  )
}

neg_f <- function(x, data, drifts, scale = 1) {
  dlba_norm(
    data,
    A = scale * alpha(x),
    b = scale * omega(x, "-"),
    t0 = x["R"],
    mean_v = scale * drifts,
    sd_v = sqrt(scale)
  )
}

pos_F <- function(x, data, drifts, scale = 1) {
  plba_norm(
    data,
    A = scale * alpha(x),
    b = scale * omega(x, "+"),
    t0 = x["R"],
    mean_v = scale * drifts,
    sd_v = sqrt(scale)
  )
}

neg_F <- function(x, data, drifts, scale = 1) {
  plba_norm(
    data,
    A = scale * alpha(x),
    b = scale * omega(x, "-"),
    t0 = x["R"],
    mean_v = scale * drifts,
    sd_v = sqrt(scale)
  )
}

# A table of calls to individual (d/p)lba_norms -------------------------------

#    func ║ Correct │ Incorrect
#   ══════╬═════════╪═══════════
#    f+   ║ 4       │
#    f-   ║         │ 4
#    F+   ║ 4       │ 4
#    F-   ║ 4       │ 4
#    f+²  ║ 3       │
#    f-²  ║         │ 3
#    F+²  ║         │ 3
#    F-²  ║ 3       │
#
#
# Enough repeats here to be worthwhile calling each once and passing resulting
# arrays around to individual model likelihood functions.

# Specify likelihoods of parameters for individual models ---------------------

# Independant parallel self terminating
ll_IST <- function(dfunc) {  # nolint
  yes <- 2 * dfunc[["f+"]] * (1 - dfunc[["F+"]]) * (1 - dfunc[["F-"]]**2)
  no <- 2 * dfunc[["f-"]] * dfunc[["F-"]] * (1 - dfunc[["F+"]])**2
  c(yes, no)
}

ll_IEX <- function(dfunc) {  # nolint
  yes <- 2 * dfunc[["f+"]] * dfunc[["F+"]] * (1 - dfunc[["F-"]])**2
  no <- 2 * dfunc[["f-"]] * (1 - dfunc[["F-"]]) * (1 - dfunc[["F+"]]**2)
  c(yes, no)
}

ll_CYST <- function(dfunc) {  # nolint
  yes <- dunc[["f+2"]] * (1 - dfunc[["F-"]]**2)
  no <- 2 * dfunc[["f-"]] * dfunc[["F-"]] *  (1 - dfunc[["F+2"]])
  c(yes, no)
}

ll_CYEX <- function(dfunc) {  # nolint
  yes <- dfunc[["f+2"]] * (1 - dfunc[["F-"]])**2
  no <- 2 * dfunc[["f-"]] * (1 - dfunc[["F-"]]) *  (1 - dfunc[["F+2"]])
  c(yes, no)
}

ll_CNST <- function(dfunc) {  # nolint
  yes <- 2 * dfunc[["f+"]] * (1 - dfunc[["F+"]]) *  (1 - dfunc[["F-2"]])
  no <- dfunc[["f-2"]] * (1 - dfunc[["F+2"]])**2
  c(yes, no)
}

ll_CNEX <- function(dfunc) {  # nolint
  yes <- 2 * dfunc[["f+"]] * dfunc[["F+"]] * (1 - dfunc[["F-2"]])
  no <- dfunc[["f-2"]] * (1 - dfunc[["F+"]]**2)
  c(yes, no)
}

ll_CB <- function(dfunc) {  # nolint
  yes <- dfunc[["f+2"]] * (1 - dfunc[["F-2"]])
  no <- dfunc[["f-2"]] * (1 - dfunc[["F+2"]])
  c(yes, no)
}

# Specify the log likelihood function -----------------------------------------
# For one accumulator, following terminology in Cox and Criss (2019)
# Start point variability is uniform 0 -> α
# Drift sampled from normal distribution with mean v and std dev s
# Threshold is ω
# Residual time is R (time to detect and start processing, and respond)
# For rtdists, we have A == α, b == ω, t0 == R, mean_v == v, sd_v == s
cox_ll <- function(x, data) {
  x <- exp(x)
  if (any(data$RT < x["R"])) {
    return(-1e10)
  }

  # MEAT OF THE mixture LL here

  bad <- (out < 1e-10) | (!is.finite(out))
  out[bad] <- 1e-10
  out <- sum(log(out))
  out
}
