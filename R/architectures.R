
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
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
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
    ll <- (dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- (dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1))
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
