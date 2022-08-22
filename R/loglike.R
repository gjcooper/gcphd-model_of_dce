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
rmodel_wrapper <- function(x, data, model, contaminant_prob = 0.02, min_rt = 0, max_rt = 1) {
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
  gen_df <- model(data, A, b_acc, b_rej, t0, drifts)
  gen_df$generator <- "model"
  # Generate contaminant responses
  for (row_idx in sample(nrow(gen_df), contaminant_prob * nrow(gen_df))) {
    gen_df[row_idx, "rt"] <- stats::runif(1, min = min_rt, max = max_rt)
    gen_df[row_idx, "accept"] <- sample.int(2, size = 1)
    gen_df[row_idx, "generator"] <- "contaminant"
  }
  gen_df
}


#' Reduced Model Wrapper for individual model log likelihood function
#'
#' This function performs some common steps such as rearrangeing the data
#' pulling out parameter items into variable and combining the results of
#' applying to model to different subsets of the data (accept and reject
#' responses). One main difference is it calculates reject drift rates using
#' accept drift rates, intercept and beta pars
#'
#' @inheritParams dirichlet_mix_ll
#' @param model The model to be wrapped and returned
#'
#' @return The log of the likelihood for the data under parameter values x
#' @export
reduced_model_wrapper <- function(x, data, model) {
  adat <- data[data$accept == 2, ]
  rdat <- data[data$accept == 1, ]

  A <- x["A"] # nolint
  t0 <- x["t0"]
  # Enforces b cannot be less than A, b parameter is thus threshold - A
  b_acc <- x["b_acc"] + A
  b_rej <- x["b_rej"] + A

  # Get accept/reject drift par namess
  v_acc_p <- apply(expand.grid("v_acc_p", stim_levels), 1, paste, collapse = "_")
  v_rej_p <- apply(expand.grid("v_rej_p", stim_levels), 1, paste, collapse = "_")
  v_acc_r <- apply(expand.grid("v_acc_r", stim_levels), 1, paste, collapse = "_")
  v_rej_r <- apply(expand.grid("v_rej_r", stim_levels), 1, paste, collapse = "_")

  # Calculate reject drift values
  x[v_rej_p] <- x["beta0_p"] - x["beta1_p"] * x[v_acc_p]
  x[v_rej_r] <- x["beta0_r"] - x["beta1_r"] * x[v_acc_r]

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

#' Reduced More Model Wrapper for individual model log likelihood function
#'
#' This function performs some common steps such as rearrangeing the data
#' pulling out parameter items into variable and combining the results of
#' applying to model to different subsets of the data (accept and reject
#' responses). One main difference is it calculates reject drift rates using
#' accept drift rates, intercept and beta pars, shared b0
#'
#' @inheritParams dirichlet_mix_ll
#' @param model The model to be wrapped and returned
#'
#' @return The log of the likelihood for the data under parameter values x
#' @export
reduced_more_model_wrapper <- function(x, data, model) {
  adat <- data[data$accept == 2, ]
  rdat <- data[data$accept == 1, ]

  A <- x["A"] # nolint
  t0 <- x["t0"]
  # Enforces b cannot be less than A, b parameter is thus threshold - A
  b_acc <- x["b_acc"] + A
  b_rej <- x["b_rej"] + A

  # Get accept/reject drift par namess
  v_acc_p <- apply(expand.grid("v_acc_p", stim_levels), 1, paste, collapse = "_")
  v_rej_p <- apply(expand.grid("v_rej_p", stim_levels), 1, paste, collapse = "_")
  v_acc_r <- apply(expand.grid("v_acc_r", stim_levels), 1, paste, collapse = "_")
  v_rej_r <- apply(expand.grid("v_rej_r", stim_levels), 1, paste, collapse = "_")

  # Calculate reject drift values
  x[v_rej_p] <- x["beta0"] - x["beta1_p"] * x[v_acc_p]
  x[v_rej_r] <- x["beta0"] - x["beta1_r"] * x[v_acc_r]

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
#'   \item A number of \strong{\eqn{\alpha}} parameter values matching the
#'     number and order of the ll_funcs vector defined in this same file.
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
#' @param contaminant_prob The probability used for contaminant process in the
#'   modelling. A contaminant process is just a uniform random response in the
#'   allowable time window.
#' @param alpha_indices A vector containing the indicies of the alpha parameters
#'   - that is the parameters that correspond to the dirichlet shape parameters.
#' @param min_rt The smallest possible response time in the data
#' @param max_rt The largest possible response time in the data
#'
#' @return The log of the likelihood for the data under parameter values x
#' @export
dirichlet_mix_ll <- function(x, data, contaminant_prob = 0.02, alpha_indices = c(1, 2), min_rt = 0, max_rt = 1) {
  x <- exp(x)

  # Enforce alphas to be greater than 0.01 and less than 100
  if (any(x[alpha_indices] < 0.01) || any(x[alpha_indices] > 100)) {
    return(-1e10)
  }

  # all decision rules
  rdev <- MCMCpack::rdirichlet(1, x[alpha_indices])
  func_idx <- sample(alpha_indices, 1, prob = rdev)
  ll_func <- ll_funcs[[func_idx]]$likelihood
  trial_ll <- model_wrapper(x, data, ll_func)
  new_like <- (1 - contaminant_prob) * trial_ll +
    contaminant_prob * (stats::dunif(data$rt, min_rt, max_rt) / 2)
  loglike <- sum(log(pmax(new_like, 1e-10)))
  attr(loglike, "extra_info") <- names_ll()[func_idx]
  loglike
}


#' Top level log-likelihood function that implements reduced dirichlet sampling
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
#'   \item A number of \strong{\eqn{\alpha}} parameter values matching the
#'     number and order of the ll_funcs vector defined in this same file.
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
#' @param contaminant_prob The probability used for contaminant process in the
#'   modelling. A contaminant process is just a uniform random response in the
#'   allowable time window.
#' @param alpha_indices A vector containing the indicies of the alpha parameters
#'   - that is the parameters that correspond to the dirichlet shape parameters.
#' @param min_rt The smallest possible response time in the data
#' @param max_rt The largest possible response time in the data
#'
#' @return The log of the likelihood for the data under parameter values x
#' @export
dirichlet_reduced_mix <- function(x, data, contaminant_prob = 0.02, alpha_indices = c(1, 2), min_rt = 0, max_rt = 1, shared=FALSE) {
  force_pos_mask <- !startsWith(names(x), "v")
  x[force_pos_mask] <- exp(x[force_pos_mask])
  if (shared) {
    wrapper_func = reduced_more_model_wrapper
  } else {
    wrapper_func = reduced_model_wrapper
  }

  # Enforce alphas to be greater than 0.01 and less than 100
  if (any(x[alpha_indices] < 0.01) || any(x[alpha_indices] > 100)) {
    return(-1e10)
  }

  # all decision rules
  rdev <- MCMCpack::rdirichlet(1, x[alpha_indices])
  func_idx <- sample(alpha_indices, 1, prob = rdev)
  ll_func <- ll_funcs[[func_idx]]$likelihood
  trial_ll <- reduced_model_wrapper(x, data, ll_func)
  new_like <- (1 - contaminant_prob) * trial_ll +
    contaminant_prob * (stats::dunif(data$rt, min_rt, max_rt) / 2)
  loglike <- sum(log(pmax(new_like, 1e-10)))
  attr(loglike, "extra_info") <- names_ll()[func_idx]
  loglike
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
#' @inheritParams dirichlet_mix_ll
#' @param architecture The name of the architecture model to implement
#'
#' @return The log of the likelihood for the data under parameter values x
#' @export
single_model_ll <- function(x, data, contaminant_prob = 0.02, architecture = "IST", min_rt = 0, max_rt = 1) {
  x <- exp(x)

  # all decision rules
  func_idx <- match(architecture, names_ll())
  ll_func <- ll_funcs[[func_idx]]$likelihood
  trial_ll <- model_wrapper(x, data, ll_func)
  new_like <- (1 - contaminant_prob) * trial_ll +
    contaminant_prob * (stats::dunif(data$rt, min_rt, max_rt) / 2)
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
