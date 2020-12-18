# Channel A is Price
# Channel B is Rating
# pdf to accept Price is: dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)
# cdf to accept Price is: plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)
# pdf to reject Price is: dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)
# cdf to rejept Price is: plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)
# pdf to accept Rating is: dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)
# cdf to accept Rating is: plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)
# pdf to reject Rating is: dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)
# cdf to rejept Rating is: plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)

#' Independant Self Terminating
#'
#' @param rt A vector of response times
#' @param A Start point variability
#' @param b_acc positive evidence threshold
#' @param b_rej negative evidence threshold
#' @param t0 non decision time parameter
#' @param drifts An array of drift rates, 4 x length(rt) with the columns being
#'   vectors of drift rates for accept price, reject price, accept rating,
#'   reject rating respectively
#' @param accept Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_IST <- function(rt, A, b_acc, b_rej, t0, drifts, accept) { # nolint
  if (accept) {
    ll <- (dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
            (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
            dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
            (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1))) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
            plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- (dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
            plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) +
            dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
            plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1))
  }
  ll
}

#' Independant Self Terminating random samples
#'
#' @param data the data for one subject
#' @inheritParams ll_IST
#'
#' @return A new data object with the same shape and new randomly drawn
#'   responses and RT's.
rll_IST <- function(data, A, b_acc, b_rej, t0, drifts) {
  for (row in seq_len(nrow(data))) {
    acc_price <- rlba_norm(1, A, b_acc, t0, drifts$AccPrice[[row]], 1)
    acc_rating <- rlba_norm(1, A, b_acc, t0, drifts$AccRating[[row]], 1)
    rej_price <- rlba_norm(1, A, b_rej, t0, drifts$RejPrice[[row]], 1)
    rej_rating <- rlba_norm(1, A, b_rej, t0, drifts$RejRating[[row]], 1)
    minacc <- min(acc_price[, "rt"], acc_rating[, "rt"])
    maxrej <- max(rej_price[, "rt"], rej_rating[, "rt"])
    if (minacc < maxrej) {
      data$rt[row] <- minacc
      data$response[row] <- 2
    } else {
      data$rt[row] <- maxrej
      data$response[row] <- 1
    }
  }
  data
}


#' Independant Exhaustive model
#'
#' @inheritParams ll_IST
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_IEX <- function(rt, A, b_acc, b_rej, t0, drifts, accept) { # nolint
  if (accept) {
    ll <- (dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
      plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) +
      dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
        plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
      (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
      (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- (dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
      (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
      dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
        (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1))) *
      (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
        plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1))
  }
  ll
}

#' Independant Exhaustive random samples
#'
#' @param data the data for one subject
#' @inheritParams ll_IST
#'
#' @return A new data object with the same shape and new randomly drawn
#'   responses and RT's.
rll_IEX <- function(data, A, b_acc, b_rej, t0, drifts) {
  data$response <- NA
  data$rt <- NA
  x <- exp(x)
  x["b_acc"] <- x["b_acc"] + x["A"]
  x["b_rej"] <- x["b_rej"] + x["A"]

  for (row in seq_len(nrow(data))) {
    acc_price <- rlba_norm(1, A, b_acc, t0, drifts$AccPrice[[row]], 1)
    acc_rating <- rlba_norm(1, A, b_acc, t0, drifts$AccRating[[row]], 1)
    rej_price <- rlba_norm(1, A, b_rej, t0, drifts$RejPrice[[row]], 1)
    rej_rating <- rlba_norm(1, A, b_rej, t0, drifts$RejRating[[row]], 1)
    minacc <- max(acc_price[, "rt"], acc_rating[, "rt"])
    maxrej <- min(rej_price[, "rt"], rej_rating[, "rt"])
    if (maxacc < minrej) {
      data$rt[row] <- minacc
      data$response[row] <- 2
    } else {
      data$rt[row] <- maxrej
      data$response[row] <- 1
    }
  }
  data
}


#' Coactive Yes Self terminating No
#'
#' @param rt A vector of response times
#' @param A Start point variability
#' @param b_acc positive evidence threshold
#' @param b_rej negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CYST <- function(rt, A, b_acc, b_rej, t0, v_pos, v_neg, accept) { # nolint
  if (accept) {
    ll <- dlba_norm(rt, 2 * A, 2 * b_acc, t0, 2 * v_pos, sqrt(2)) * # nolint
      (1 - plba_norm(rt, A, b_rej, t0, v_neg, 1)**2)
  } else {
    ll <- 2 *
      dlba_norm(rt, A, b_rej, t0, v_neg, 1) *
      plba_norm(rt, A, b_rej, t0, v_neg, 1) *
      (1 - plba_norm(rt, 2 * A, 2 * b_acc, t0, 2 * v_pos, sqrt(2))) # nolint
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
#' @param b_acc positive evidence threshold
#' @param b_rej negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CYEX <- function(rt, A, b_acc, b_rej, t0, v_pos, v_neg, accept) { # nolint
  if (accept) {
    ll <- dlba_norm(rt, 2 * A, 2 * b_acc, t0, 2 * v_pos, sqrt(2)) * # nolint
      (1 - plba_norm(rt, A, b_rej, t0, v_neg, 1))**2
  } else {
    ll <- 2 *
      dlba_norm(rt, A, b_rej, t0, v_neg, 1) *
      (1 - plba_norm(rt, A, b_rej, t0, v_neg, 1)) *
      (1 - plba_norm(rt, 2 * A, 2 * b_acc, t0, 2 * v_pos, sqrt(2))) # nolint
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
#' @param b_acc positive evidence threshold
#' @param b_rej negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CNST <- function(rt, A, b_acc, b_rej, t0, v_pos, v_neg, accept) { # nolint
  if (accept) {
    ll <- 2 *
      dlba_norm(rt, A, b_acc, t0, v_pos, 1) *
      (1 - plba_norm(rt, A, b_acc, t0, v_pos, 1)) *
      (1 - plba_norm(rt, 2 * A, 2 * b_rej, t0, 2 * v_neg, sqrt(2))) # nolint
  } else {
    ll <- dlba_norm(rt, 2 * A, 2 * b_rej, t0, 2 * v_neg, sqrt(2)) * # nolint
      (1 - plba_norm(rt, A, b_acc, t0, v_pos, 1))**2
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
#' @param b_acc positive evidence threshold
#' @param b_rej negative evidence threshold
#' @param t0 non decision time parameter
#' @param v_pos positive evidence drift rate (vector of)
#' @param v_neg vector of negative evidence drift rates
#' @param positive Whether we are looking at accept or reject trials
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CNEX <- function(rt, A, b_acc, b_rej, t0, v_pos, v_neg, accept) { # nolint
  if (accept) {
    ll <- 2 *
      dlba_norm(rt, A, b_acc, t0, v_pos, 1) *
      plba_norm(rt, A, b_acc, t0, v_pos, 1) *
      (1 - plba_norm(rt, 2 * A, 2 * b_rej, t0, 2 * v_neg, sqrt(2))) # nolint
  } else {
    ll <- dlba_norm(rt, 2 * A, 2 * b_rej, t0, 2 * v_neg, sqrt(2)) * # nolint
      (1 - plba_norm(rt, A, b_acc, t0, v_pos, 1)**2)
  }
  ll
}
rll_CNEX <- function(x, data) {
  stop("Not implemented yet")
}


#' Coactive Both
#'
#' @inheritParams ll_IST
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_CB <- function(rt, A, b_acc, b_rej, t0, drifts, accept) { # nolint
  acc_co_drifts <- drifts$AccPrice + drifts$AccRating
  rej_co_drifts <- drifts$RejPrice + drifts$RejRating
  if (accept) {
    ll <- dlba_norm(rt, 2 * A, 2 * b_acc, t0, acc_co_drifts, sqrt(2)) * # nolint
      (1 - plba_norm(rt, 2 * A, 2 * b_rej, t0, rej_co_drifts, sqrt(2))) # nolint
  } else {
    ll <- dlba_norm(rt, 2 * A, 2 * b_rej, t0, rej_co_drifts, sqrt(2)) * # nolint
      (1 - plba_norm(rt, 2 * A, 2 * b_acc, t0, acc_co_drifts, sqrt(2))) # nolint
  }
  ll
}

#' Coactive Both random samples
#'
#' @param data the data for one subject
#' @inheritParams ll_IST
#'
#' @return A new data object with the same shape and new randomly drawn
#'   responses and RT's.
rll_CB <- function(data, A, b_acc, b_rej, t0, drifts) {
  data$response <- NA
  data$rt <- NA
  x <- exp(x)
  x["b_acc"] <- x["b_acc"] + x["A"]
  x["b_rej"] <- x["b_rej"] + x["A"]
  acc_co_drifts <- drifts$AccPrice + drifts$AccRating
  rej_co_drifts <- drifts$RejPrice + drifts$RejRating

  for (row in seq_len(nrow(data))) {
    acc_coactive <- rlba_norm(1, 2*A, 2*b_acc, t0, acc_co_drifts[[row]], sqrt(2))
    rej_coactive <- rlba_norm(1, 2*A, 2*b_rej, t0, rej_co_drifts[[row]], sqrt(2))
    if (acc_coactive < rej_coactive) {
      data$rt[row] <- acc_coactive
      data$response[row] <- 2
    } else {
      data$rt[row] <- rej_coactive
      data$response[row] <- 1
    }
  }
  data
}

ll_funcs <- c(ll_IST, ll_IEX, ll_CB)
rll_funcs <- c(rll_IST, rll_IEX, rll_CB)
ll_names <- c("IST", "IEX", "CB")


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
  adat <- data[data$accept == 2, ]
  rdat <- data[data$accept == 1, ]


  A <- x["A"] # nolint
  t0 <- x["t0"]
  # Enforces b cannot be less than A, b parameter is thus threshold - A
  b_acc <- x["b_acc"] + A
  b_rej <- x["b_rej"] + A
  acc_drifts <- tibble(
    AccPrice = x[adat$v_acc_p],
    RejPrice = x[adat$v_rej_p],
    AccRating = x[adat$v_acc_r],
    RejRating = x[adat$v_rej_r]
  )
  rej_drifts <- tibble(
    AccPrice = x[rdat$v_acc_p],
    RejPrice = x[rdat$v_rej_p],
    AccRating = x[rdat$v_acc_r],
    RejRating = x[rdat$v_rej_r]
  )

  accept <- switch(
    nrow(adat) != 0,
    model(adat$rt, A, b_acc, b_rej, t0, acc_drifts, TRUE)
  )
  reject <- switch(
    nrow(rdat) != 0,
    model(rdat$rt, A, b_acc, b_rej, t0, rej_drifts, FALSE)
  )
  c(accept, reject)
}

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
rmodel_wrapper <- function(x, data, model) {
  data$response <- NA
  data$rt <- NA
  x <- exp(x)

  A <- x["A"] # nolint
  t0 <- x["t0"]
  # Enforces b cannot be less than A, b parameter is thus threshold - A
  b_acc <- x["b_acc"] + A
  b_rej <- x["b_rej"] + A
  drifts <- tibble(
    AccPrice = x[data$v_acc_p],
    RejPrice = x[data$v_rej_p],
    AccRating = x[data$v_acc_r],
    RejRating = x[data$v_rej_r]
  )
  # Return newly generated data
  model(data, A, b_acc, b_rej, t0, drifts)
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

  # all decision rules
  rdev <- MCMCpack::rdirichlet(1, x[mix_counts])
  func_idx <- sample(mix_counts, 1, prob = rdev)
  ll_func <- ll_funcs[[func_idx]]
  trial_ll <- model_wrapper(x, data, ll_func)
  new_like <- (1 - p_contam) * trial_ll +
    p_contam * (stats::dunif(data$rt, min_rt, max_rt) / 2)
  sum(log(pmax(new_like, 1e-10)))
}


#' Top level log-likelihood function that implements the single model sampling
#'
#' This function runs one of the possible architectures repeatedly, horrible,
#' uses a global variable for the selected model.
#'
#' @section The parameter vector:
#'
#' The vector x should contain the following elements:
#' A number of \eqn{\alpha}
#' values
#'
#' \itemize{
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
single_model_ll <- function(x, data) {
  x <- exp(x)

  # all decision rules
  func_idx <- match(architecture, ll_names)
  ll_func <- ll_funcs[[func_idx]]
  trial_ll <- model_wrapper(x, data, ll_func)
  new_like <- (1 - p_contam) * trial_ll +
    p_contam * (stats::dunif(data$rt, min_rt, max_rt) / 2)
  sum(log(pmax(new_like, 1e-10)))
}
