require(psamplers)
require(rtdists)
source("read_expyriment.R")

task1_data <- read.expyriment.data("data/input/Task1/", "S*")

#Sum and difference of evidence rates for positive and negative accumulators
sum_diff <- c("beta", "delta")
stim_levels <- c("high", "low", "distract")
channels <- c("price", "quality")

parameters <- c(
  # Parallel mixture probabilities
  "lambda_IST", "lambda_IEX",
  # Coactive mixture probabilities
  "lambda_CYST", "lambda_CYEX", "lambda_CNST", "lambda_CNEX", "lambda_CB",
  # ω - Sum of response thresholds
  "omega",
  # ν - how much threshold is subject to start point variability
  "nu",
  # w - proportion of total threshold given to negative accunmulator
  "w",
  # R - residual time, bounded above by the minimum observed response time for participant k
  "R"
)
# beta and delta - 12 versions each for each of the 12 cells of the experimental design
parameters <- c(parameters, apply(expand.grid(sum_diff, stim_levels, channels), 1, paste, collapse="_"))


# Specify the log likelihood function -----------------------------------------
# NOT YET IMPLEMENTED
cox_ll <- function(x, data, sample = FALSE) {
  x <- exp(x)
  if (any(data$rt < x["t0"])) {
    return(-1e10)
  }
  # This is faster than "paste".
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$condition]
  
  out <- rtdists::dLBA(rt = data$rt, # nolint
                       response = data$correct,
                       A = x["A"],
                       b = bs,
                       t0 = x["t0"],
                       mean_v = x[c("v1", "v2")],
                       sd_v = c(1, 1),
                       distribution = "norm",
                       silent = TRUE)
  bad <- (out < 1e-10) | (!is.finite(out))
  out[bad] <- 1e-10
  out <- sum(log(out))
  out
}