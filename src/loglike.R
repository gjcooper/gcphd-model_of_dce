library(MCMCpack)

# Simulation for sampling

# Specify likelihoods of parameters for individual models ---------------------

# Independant parallel self terminating
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
#' dirichlet sample weighted by the \alpha values from the parameter vector. It
#' then delegates the likelihood for the data to an architecture specific
#' log-likelihood function.
#'
#' @section The parameter vector:
#'
#' The vector x should contain the following elements:
#' A number of /alpha
#' values 
#' 
#' \itemize{
#'   \item A number of \strong{\alpha} parameter values matching the number and
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
#' @return The log of the lilelihood for the data under paramweter values x
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
  rdev <- rdirichlet(1, x[mix_counts])
  func_idx <- sample(mix_counts, 1, prob = rdev)
  ll_func <- ll_funcs[[func_idx]]
  trial_ll <- model_wrapper(x, data, ll_func)
  new_like <- (1 - p_contam) * trial_ll +
    p_contam * (dunif(data$rt, min_rt, max_rt) / 2)
  sum(log(pmax(new_like, 1e-10)))
}

test.likelihood <- function(x, num.values = 9, fake.data, dimensions, ll_func, server = FALSE, ...) {
  pars <- names(x)
  x.tmp <- x
  sequence <- seq(from = -.22, to = .22, length.out = num.values)
  op <- par(mfrow = dimensions)
  par_likelihoods <- lapply(setNames(pars, pars), function(p) {
    testvalues <- sequence + true.x[[p]]
    tmp <- unlist(lapply(testvalues, function(i) { # for each test vlaue, it will apply the function where i = testvalue.
      x.tmp[[p]] <- i
      ll_func(x = x.tmp, data = fake.data)
    }))

    plot(
      x = testvalues,
      y = tmp,
      type = "b",
      main = p,
      xlab = "log par values",
      ylab = "loglikelihood",
      ...
    )
    abline(v = true.x[[p]], col = "red")
    return(tmp)
  })
  par(op)
  return(par_likelihoods)
}
