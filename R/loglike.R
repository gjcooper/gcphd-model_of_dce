#' Independant Self Terminating
#'
#' @param rt A vector of response times
#' @param A Start point variability 
#' @param b_pos positive evidence threshold
#' @param b_neg negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_IST <- function(rt, A, b_pos, b_neg, t0, v_pos, v_neg, positive) { # nolint
  if (positive) {
    ll <- 2 *
      dlba_norm(rt, A, b_pos, t0, v_pos, 1) *
      (1 - plba_norm(rt, A, b_pos, t0, v_pos, 1)) *
      (1 - plba_norm(rt, A, b_neg, t0, v_neg, 1)**2)
  } else {
    ll <- 2 *
      dlba_norm(rt, A, b_neg, t0, v_neg, 1) *
      plba_norm(rt, A, b_neg, t0, v_neg, 1) *
      (1 - plba_norm(rt, A, b_pos, t0, v_pos, 1))**2
  }
  ll
}
rll_IST <- function(x, data) {
  data$response <- NA
  data$rt <- NA
  x <- exp(x)
  x["b_pos"] <- x["b_pos"] + x["A"]
  x["b_neg"] <- x["b_neg"] + x["A"]

  for (row in seq_len(nrow(data))) {
    pos <- rlba_norm(2, x[["A"]], x[["b_pos"]], x[["t0"]], x[[data$v_pos[row]]], 1)
    neg <- rlba_norm(2, x[["A"]], x[["b_neg"]], x[["t0"]], x[[data$v_neg[row]]], 1)
    minpos <- min(pos[, "rt"])
    maxneg <- max(neg[, "rt"])
    if (minpos < maxneg) {
      data$rt[row] <- minpos
      data$response[row] <- 2
    } else {
      data$rt[row] <- maxneg
      data$response[row] <- 1
    }
  }
  data
}


#' Independant Exhaustive model
#'
#' @param rt A vector of response times
#' @param A Start point variability 
#' @param b_pos positive evidence threshold
#' @param b_neg negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_IEX <- function(rt, A, b_pos, b_neg, t0, v_pos, v_neg, positive) { # nolint
  if (positive) {
    ll <- 2 *
      dlba_norm(rt, A, b_pos, t0, v_pos, 1) *
      plba_norm(rt, A, b_pos, t0, v_pos, 1) *
      (1 - plba_norm(rt, A, b_neg, t0, v_neg, 1))**2
  } else {
    ll <- 2 *
      dlba_norm(rt, A, b_neg, t0, v_neg, 1) *
      (1 - plba_norm(rt, A, b_neg, t0, v_neg, 1)) *
      (1 - plba_norm(rt, A, b_pos, t0, v_pos, 1)**2)
  }
  ll
}
rll_IEX <- function(x, data) {
  stop("Not implemented yet")
}


#' Coactive Yes Self terminating No
#'
#' @param rt A vector of response times
#' @param A Start point variability 
#' @param b_pos positive evidence threshold
#' @param b_neg negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CYST <- function(rt, A, b_pos, b_neg, t0, v_pos, v_neg, positive) { # nolint
  if (positive) {
    ll <- dlba_norm(rt, 2 * A, 2 * b_pos, t0, 2 * v_pos, sqrt(2)) * # nolint
      (1 - plba_norm(rt, A, b_neg, t0, v_neg, 1)**2)
  } else {
    ll <- 2 *
      dlba_norm(rt, A, b_neg, t0, v_neg, 1) *
      plba_norm(rt, A, b_neg, t0, v_neg, 1) *
      (1 - plba_norm(rt, 2 * A, 2 * b_pos, t0, 2 * v_pos, sqrt(2))) # nolint
  }
  ll
}
rll_CYST <- function(x, data) {
  stop("Not implemented yet")
}


#' Coactive Yes Exhaustive No
#'
#' @param rt A vector of response times
#' @param A Start point variability 
#' @param b_pos positive evidence threshold
#' @param b_neg negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CYEX <- function(rt, A, b_pos, b_neg, t0, v_pos, v_neg, positive) { # nolint
  if (positive) {
    ll <- dlba_norm(rt, 2 * A, 2 * b_pos, t0, 2 * v_pos, sqrt(2)) * # nolint
      (1 - plba_norm(rt, A, b_neg, t0, v_neg, 1))**2
  } else {
    ll <- 2 *
      dlba_norm(rt, A, b_neg, t0, v_neg, 1) *
      (1 - plba_norm(rt, A, b_neg, t0, v_neg, 1)) *
      (1 - plba_norm(rt, 2 * A, 2 * b_pos, t0, 2 * v_pos, sqrt(2))) # nolint
  }
  ll
}
rll_CYEX <- function(x, data) {
  stop("Not implemented yet")
}


#' Coactive No/Self terminating Yes
#'
#' @param rt A vector of response times
#' @param A Start point variability 
#' @param b_pos positive evidence threshold
#' @param b_neg negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CNST <- function(rt, A, b_pos, b_neg, t0, v_pos, v_neg, positive) { # nolint
  if (positive) {
    ll <- 2 *
      dlba_norm(rt, A, b_pos, t0, v_pos, 1) *
      (1 - plba_norm(rt, A, b_pos, t0, v_pos, 1)) *
      (1 - plba_norm(rt, 2 * A, 2 * b_neg, t0, 2 * v_neg, sqrt(2))) # nolint
  } else {
    ll <- dlba_norm(rt, 2 * A, 2 * b_neg, t0, 2 * v_neg, sqrt(2)) * # nolint
      (1 - plba_norm(rt, A, b_pos, t0, v_pos, 1))**2
  }
  ll
}
rll_CNST <- function(x, data) {
  stop("Not implemented yet")
}


#' Coactive No / Exhaustive Yes
#'
#' @param rt A vector of response times
#' @param A Start point variability 
#' @param b_pos positive evidence threshold
#' @param b_neg negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CNEX <- function(rt, A, b_pos, b_neg, t0, v_pos, v_neg, positive) { # nolint
  if (positive) {
    ll <- 2 *
      dlba_norm(rt, A, b_pos, t0, v_pos, 1) *
      plba_norm(rt, A, b_pos, t0, v_pos, 1) *
      (1 - plba_norm(rt, 2 * A, 2 * b_neg, t0, 2 * v_neg, sqrt(2))) # nolint
  } else {
    ll <- dlba_norm(rt, 2 * A, 2 * b_neg, t0, 2 * v_neg, sqrt(2)) * # nolint
      (1 - plba_norm(rt, A, b_pos, t0, v_pos, 1)**2)
  }
  ll
}
rll_CNEX <- function(x, data) {
  stop("Not implemented yet")
}


#' Coactive Both
#'
#' @param rt A vector of response times
#' @param A Start point variability 
#' @param b_pos positive evidence threshold
#' @param b_neg negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CB <- function(rt, A, b_pos, b_neg, t0, v_pos, v_neg, positive) { # nolint
  if (positive) {
    ll <- dlba_norm(rt, 2 * A, 2 * b_pos, t0, 2 * v_pos, sqrt(2)) * # nolint
      (1 - plba_norm(rt, 2 * A, 2 * b_neg, t0, 2 * v_neg, sqrt(2))) # nolint
  } else {
    ll <- dlba_norm(rt, 2 * A, 2 * b_neg, t0, 2 * v_neg, sqrt(2)) * # nolint
      (1 - plba_norm(rt, 2 * A, 2 * b_pos, t0, 2 * v_pos, sqrt(2))) # nolint
  }
  ll
}
rll_CB <- function(x, data) {
  stop("Not implemented yet")
}

ll_funcs <- c(ll_IST, ll_IEX, ll_CYST, ll_CYEX, ll_CNST, ll_CNEX, ll_CB)
rll_funcs <- c(rll_IST, rll_IEX, rll_CYST, rll_CYEX, rll_CNST, rll_CNEX, rll_CB)
ll_names <- c("IST", "IEX", "CYST", "CYEX", "CNST", "CNEX", "CB")


#' Wrapper for individual model log likelihood function
#'
#' This function performs some common steps such as rearrangeing the data
#' pulling out parameter items into variable and combining the results of
#' applying to model to different subsets of the data (accept and reject
#' responses)
#'
#' @inheritParams dirichlet_mix_ll
#' @param model The model to be wrapped and returned
#'
#' @return The log of the likelihood for the data under parameter values x
model_wrapper <- function(x, data, model) {
  ydat <- data[data$response == 2, ]
  ndat <- data[data$response == 1, ]
  A <- x["A"] # nolint
  t0 <- x["t0"]
  b_pos <- x["b_pos"]
  b_neg <- x["b_neg"]

  yes <- switch(
    nrow(ydat) != 0,
    model(ydat$rt, A, b_pos, b_neg, t0, x[ydat$v_pos], x[ydat$v_neg], 1)
  )
  no <- switch(
    nrow(ndat) != 0,
    model(ndat$rt, A, b_pos, b_neg, t0, x[ndat$v_pos], x[ndat$v_neg], 1)
  )
  c(yes, no)
}


#' Top level log-likelihood function that implements the dirichlet sampling
#'
#' This function selects one of the possible architectures using a random
#' dirichlet sample weighted by the α values from the parameter vector. It
#' then delegates the likelihood for the data to an architecture specific
#' log-likelihood function.
#'
#' @section The parameter vector:
#'
#' The vector x should contain the following elements:
#' A number of \eqn{\alpha}
#' values 
#' 
#' \itemize{
#'   \item A number of \strong{\eqn{\alpha}} parameter values matching the number and
#'     order of the ll_funcs vector defined in this same file.
#'   \item \strong{A} - the start point variability
#'   \item \strong{b^a} and \strong{b^r}, the thresholds to either accept or
#'     reject the item.
#'   \item \strong{t0} - the residual time, bounded above by the minimum
#'     response time for the participant
#'   \item 12 drift rates. For each attribute there are three stimulus levels.
#'     For each of these 6 attribute levels there are two drift rates, one drift
#'     rate to accept (\strong{v^a}) and one to reject (\strong{v^r})
#' }
#'
#' @param x A named vector containing parameter values to test
#' @param data The data for a single subject for which the likelihood should be
#'   calculated
#'
#' @return The log of the likelihood for the data under parameter values x
dirichlet_mix_ll <- function(x, data) {
  x <- exp(x)

  # Enforce alphas to be greater than 0.01 and less than 100
  if (any(x[mix_counts] < 0.01) || any(x[mix_counts] > 100)) {
    return(-1e10)
  }

  # Enforces b cannot be less than A, b parameter is thus threshold - A
  x["b_pos"] <- x["b_pos"] + x["A"]
  x["b_neg"] <- x["b_neg"] + x["A"]

  # all decision rules
  rdev <- MCMCpack::rdirichlet(1, x[mix_counts])
  func_idx <- sample(mix_counts, 1, prob = rdev)
  ll_func <- ll_funcs[[func_idx]]
  trial_ll <- model_wrapper(x, data, ll_func)
  new_like <- (1 - p_contam) * trial_ll +
    p_contam * (stats::dunif(data$rt, min_rt, max_rt) / 2)
  sum(log(pmax(new_like, 1e-10)))
}
