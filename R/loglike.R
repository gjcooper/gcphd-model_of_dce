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
      data$accept[row] <- 2
    } else {
      data$rt[row] <- maxrej
      data$accept[row] <- 1
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
#' @inheritParams rll_IST
#'
#' @return A new data object with the same shape and new randomly drawn
#'   responses and RT's.
rll_IEX <- function(data, A, b_acc, b_rej, t0, drifts) {
  for (row in seq_len(nrow(data))) {
    acc_price <- rlba_norm(1, A, b_acc, t0, drifts$AccPrice[[row]], 1)
    acc_rating <- rlba_norm(1, A, b_acc, t0, drifts$AccRating[[row]], 1)
    rej_price <- rlba_norm(1, A, b_rej, t0, drifts$RejPrice[[row]], 1)
    rej_rating <- rlba_norm(1, A, b_rej, t0, drifts$RejRating[[row]], 1)
    maxacc <- max(acc_price[, "rt"], acc_rating[, "rt"])
    minrej <- min(rej_price[, "rt"], rej_rating[, "rt"])
    if (maxacc < minrej) {
      data$rt[row] <- maxacc
      data$accept[row] <- 2
    } else {
      data$rt[row] <- minrej
      data$accept[row] <- 1
    }
  }
  data
}


#' First Past the Post
#'
#' @inheritParams ll_IST
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_FPP <- function(rt, A, b_acc, b_rej, t0, drifts, accept) { # nolint
  if (accept) {
    ll <- dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) ) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) ) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) ) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1))
  }
  ll
}

#' First Past the Post random samples
#'
#' @inheritParams rll_IST
#'
#' @return A new data object with the same shape and new randomly drawn
#'   responses and RT's.
rll_FPP <- function(data, A, b_acc, b_rej, t0, drifts) {
  for (row in seq_len(nrow(data))) {
    acc_price <- rlba_norm(1, A, b_acc, t0, drifts$AccPrice[[row]], 1)
    acc_rating <- rlba_norm(1, A, b_acc, t0, drifts$AccRating[[row]], 1)
    rej_price <- rlba_norm(1, A, b_rej, t0, drifts$RejPrice[[row]], 1)
    rej_rating <- rlba_norm(1, A, b_rej, t0, drifts$RejRating[[row]], 1)
    minacc <- min(acc_price[, "rt"], acc_rating[, "rt"])
    minrej <- min(rej_price[, "rt"], rej_rating[, "rt"])
    if (minacc < minrej) {
      data$rt[row] <- minacc
      data$accept[row] <- 2
    } else {
      data$rt[row] <- minrej
      data$accept[row] <- 1
    }
  }
  data
}


#' Max Winner model
#'
#' A model where both accumulators for Price and Rating must terminate before
#' either the option is accepted or rejected.
#'
#' @inheritParams ll_IST
#'
#' @return The log likelihood of the rts for the accept or reject trials given
#'   the provided parameter values
ll_MW <- function(rt, A, b_acc, b_rej, t0, drifts, accept) { # nolint
  if (accept) {
    ll <- dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)  *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1))
  }
  ll
}

#' Max Winner random samples
#'
#' @inheritParams rll_IST
#'
#' @return A new data object with the same shape and new randomly drawn
#'   responses and RT's.
rll_MW <- function(data, A, b_acc, b_rej, t0, drifts) {
  for (row in seq_len(nrow(data))) {
    acc_price <- rlba_norm(1, A, b_acc, t0, drifts$AccPrice[[row]], 1)
    acc_rating <- rlba_norm(1, A, b_acc, t0, drifts$AccRating[[row]], 1)
    rej_price <- rlba_norm(1, A, b_rej, t0, drifts$RejPrice[[row]], 1)
    rej_rating <- rlba_norm(1, A, b_rej, t0, drifts$RejRating[[row]], 1)
    maxacc <- max(acc_price[, "rt"], acc_rating[, "rt"])
    maxrej <- max(rej_price[, "rt"], rej_rating[, "rt"])
    if (maxacc < maxrej) {
      data$rt[row] <- maxacc
      data$accept[row] <- 2
    } else {
      data$rt[row] <- maxrej
      data$accept[row] <- 1
    }
  }
  data
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
#' @inheritParams rll_IST
#'
#' @return A new data object with the same shape and new randomly drawn
#'   responses and RT's.
rll_CB <- function(data, A, b_acc, b_rej, t0, drifts) {
  acc_co_drifts <- drifts$AccPrice + drifts$AccRating
  rej_co_drifts <- drifts$RejPrice + drifts$RejRating
  for (row in seq_len(nrow(data))) {
    acc_coactive <- rlba_norm(1, 2*A, 2*b_acc, t0, acc_co_drifts[[row]], sqrt(2))
    rej_coactive <- rlba_norm(1, 2*A, 2*b_rej, t0, rej_co_drifts[[row]], sqrt(2))
    acc_coactive_time <- acc_coactive[,"rt"]
    rej_coactive_time <- rej_coactive[,"rt"]
    if (acc_coactive_time < rej_coactive_time) {
      data$rt[row] <- acc_coactive_time
      data$accept[row] <- 2 # Accept the option
    } else {
      data$rt[row] <- rej_coactive_time
      data$accept[row] <- 1 # Reject the option
    }
  }
  data
}


ll_funcs <- list(
  IST = list(likelihood = ll_IST, sample = rll_IST),
  IEX = list(likelihood = ll_IEX, sample = rll_IEX),
  FPP = list(likelihood = ll_FPP, sample = rll_FPP),
  MW = list(likelihood = ll_MW, sample = rll_MW),
  CB = list(likelihood = ll_CB, sample = rll_CB)
)


#' Select log likelihood function
#'
#' @param name The name of the log likelihpood function to select
#' @param sample If TRUE then return the sample function instead of the density
#'   function
#'
#' @return the requested function
#' @export
select_ll <- function(name, sample = FALSE) {
  if (sample) {
    return(ll_funcs[[name]]$sample)
  }
  ll_funcs[[name]]$likelihood
}

#' Return names of all implemented models (log likelihood functions)
#'
#' @return names of the ll_funcs package object
#' @export
names_ll <- function() {
  names(ll_funcs)
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
#' @export
model_wrapper <- function(x, data, model) {
  adat <- data[data$accept == 2, ]
  rdat <- data[data$accept == 1, ]


  A <- x["A"] # nolint
  t0 <- x["t0"]
  # Enforces b cannot be less than A, b parameter is thus threshold - A
  b_acc <- x["b_acc"] + A
  b_rej <- x["b_rej"] + A
  acc_drifts <- tibble::tibble(
    AccPrice = x[adat$v_acc_p],
    RejPrice = x[adat$v_rej_p],
    AccRating = x[adat$v_acc_r],
    RejRating = x[adat$v_rej_r]
  )
  rej_drifts <- tibble::tibble(
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
#' @export
rmodel_wrapper <- function(x, data, model) {
  data$accept <- NA
  data$rt <- NA
  x <- exp(x)

  A <- x["A"] # nolint
  t0 <- x["t0"]
  # Enforces b cannot be less than A, b parameter is thus threshold - A
  b_acc <- x["b_acc"] + A
  b_rej <- x["b_rej"] + A
  drifts <- tibble::tibble(
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
#' @export
dirichlet_mix_ll <- function(x, data) {
  x <- exp(x)

  # Enforce alphas to be greater than 0.01 and less than 100
  if (any(x[mix_counts] < 0.01) || any(x[mix_counts] > 100)) {
    return(-1e10)
  }

  # all decision rules
  rdev <- MCMCpack::rdirichlet(1, x[mix_counts])
  func_idx <- sample(mix_counts, 1, prob = rdev)
  ll_func <- ll_funcs[[func_idx]]$likelihood
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
#' @export
single_model_ll <- function(x, data) {
  x <- exp(x)

  # all decision rules
  func_idx <- match(architecture, names_ll())
  ll_func <- ll_funcs[[func_idx]]$likelihood
  trial_ll <- model_wrapper(x, data, ll_func)
  new_like <- (1 - p_contam) * trial_ll +
    p_contam * (stats::dunif(data$rt, min_rt, max_rt) / 2)
  sum(log(pmax(new_like, 1e-10)))
}

#' Simple wrapper for integrate testing
#'
#' This function does the parameter transformations necessary for the model to
#' run.
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
#' @inheritParams ll_IST
#' @param model The loglikelihood function currently being tested
#'
#' @return The log of the likelihood for the rt's for the parameters from the
#'   model
#' @export
simple_model_wrapper <- function(rt, A, b_acc, b_rej, t0, drifts, accept, model) {
  A <- exp(A)
  b_acc <- exp(b_acc) + A
  b_rej <- exp(b_rej) + A
  t0 <- exp(t0)
  drifts <- exp(drifts)

  model(rt, A, b_acc, b_rej, t0, drifts, accept)
}
