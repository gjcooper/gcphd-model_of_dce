# Useful names for drift rates used internally in functions
acc_rej_drift <- c("v_acc_p", "v_acc_r", "v_rej_p", "v_rej_r")
stim_levels <- c("H", "L", "D")
drift_names <- apply(expand.grid(acc_rej_drift, stim_levels), 1, paste, collapse = "_")
v_acc_p <- drift_names[grepl("acc_p", drift_names)]
v_rej_p <- drift_names[grepl("rej_p", drift_names)]
v_acc_r <- drift_names[grepl("acc_r", drift_names)]
v_rej_r <- drift_names[grepl("rej_r", drift_names)]


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
  drifts <- tibble::tibble(
    AccPrice = x[data$v_acc_p],
    RejPrice = x[data$v_rej_p],
    AccRating = x[data$v_acc_r],
    RejRating = x[data$v_rej_r]
  )

  accepts <- data$accept == 2
  acc_args <- c(
    list(rt = data[accepts,]$rt, drifts = drifts[accepts, ], accept = TRUE),
    as.list(x[c('A', 'b_acc', 'b_rej', 't0')]))

  rejects <- data$accept == 1
  rej_args <- c(
    list(rt = data[rejects,]$rt, drifts = drifts[rejects, ], accept = FALSE),
    as.list(x[c('A', 'b_acc', 'b_rej', 't0')]))

  accept <- switch(
    sum(accepts) != 0,
    do.call(model, acc_args)
  )
  reject <- switch(
    sum(rejects) != 0,
    do.call(model, rej_args)
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
  x <- transform_pars(x)
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
#' @inheritParams transform_pars
#'
#' @return The log of the likelihood for the data under parameter values x
#' @export
dirichlet_mix_ll <- function(x, data, contaminant_prob = 0.02, alpha_indices = c(1, 2), min_rt = 0, max_rt = 1, tforms = "std") {
  x <- transform_pars(x, data, tforms)

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
  x <- transform_pars(x, data)

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

#' Transforms parameters for specific models/
#'
#' This function performs common transforms for the LBA based architecture of
#' choice modelling
#'
#' @section Transformation schemes:
#'
#' Transformations scheme currently allowed are provided below.
#' \itemize{
#'   \item std - This is the default, and takes the exponent of all
#'     parameter values, enforces \strong{b_acc} and \strong{b_rej} being
#'     greater than \strong{A} by adding \strong{A}
#'   \item reduced - This works similarly to std, except it does not take the
#'     exponent of drift rates (parameters starting with v). It also calculates
#'     drift rate parameters for reject accumulators from the accept drift rates
#'     combined with a beta_0 and beta_1 parameters for each attribute.
#'   \item reduced_more - This works similarly to reduced, however the beta_0
#'     parameters is shared between attributes.
#' }
#'
#' @param pars The named list of parameters to transform
#' @param tforms A character vector selecting the transformation scheme
#'   scheme to be used.
#'
#' @return A named list of transformed parameter values
#' @export
transform_pars <- function(pars, data, tforms = "std") {
  if (tforms == "std") {
    newpars <- exp(pars)
  } else if (startsWith(tforms, "reduced")) {
    force_pos_mask <- !startsWith(names(pars), "v")
    newpars <- pars
    newpars[force_pos_mask] <- exp(newpars[force_pos_mask])
    if (tforms == "reduced") {
      # Calculate reject drift values
      newpars[v_rej_p] <- newpars["beta0_p"] - newpars["beta1_p"] * newpars[v_acc_p]
      newpars[v_rej_r] <- newpars["beta0_r"] - newpars["beta1_r"] * newpars[v_acc_r]
    } else if (tforms == "reduced_more") {
      # Calculate reject drift values
      newpars[v_rej_p] <- newpars["beta0"] - newpars["beta1_p"] * newpars[v_acc_p]
      newpars[v_rej_r] <- newpars["beta0"] - newpars["beta1_r"] * newpars[v_acc_r]
    }
  }
  # Always need to adjust for threshold > A
  newpars["b_acc"] <- newpars["b_acc"] + newpars["A"]
  newpars["b_rej"] <- newpars["b_rej"] + newpars["A"]
  newpars
}
