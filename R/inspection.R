consolidate_filter <- function(sampler, filter) {
  samples <- sampler$samples

  if (test_character(filter)) {
    assert_subset(filter, c("init", "burn", "adapt", "sample"))
    filter <- which(samples$stage %in% filter)
  } else {
    assert_subset(filter, 1:samples$idx)
  }
  filter
}


#' Extract group level parameters from the samples
#'
#' This function taks a pmwgs sampler object and extracts the group level
#' samples for the specified parameter estimates.
#'
#' @param sampler The pmwgs sampler object
#' @param par_names The names of the parameters to extract - defaults to all
#'   parameters estimated.
#' @param filter The sampling stages to extract - defaults to all stages run.
#'   Can also be a vector of sample indices.
#'
#' @return A tibble with the parameter samples and a sampleid column
#'
#' @import dplyr
#' @export
extract_tmu <- function(sampler,
                        par_names = sampler$par_names,
                        filter = unique(sampler$samples$stage)) {
  filter <- consolidate_filter(sampler, filter)
  stages <- sampler$samples$stage[filter]

  sampler %>%
    pmwg::as_mcmc(filter = filter) %>%
    as_tibble() %>%
    select(all_of(par_names)) %>%
    mutate(sampleid = row_number()) %>%
    mutate(stageid = stages)
}

#' Extract subject level parameters from the samples
#'
#' This function taks a pmwgs sampler object and extracts the subject level
#' samples for the specified parameter estimates.
#'
#' @inheritParams extract_tmu
#'
#' @return A tibble with the parameter samples and a sampleid + subjectid column
#'
#' @import dplyr
#' @export
extract_alpha <- function(sampler,
                          par_names = sampler$par_names,
                          filter = unique(sampler$samples$stage)) {
  filter <- consolidate_filter(sampler, filter)
  stages <- sampler$samples$stage[filter]

  pmwg::as_mcmc(sampler, selection = "alpha", filter = filter) %>%
    lapply(FUN = function(x) {
      x %>%
        as_tibble() %>%
        select(all_of(par_names)) %>%
        mutate(sampleid = row_number()) %>%
        mutate(stageid = stages)
    }) %>%
    bind_rows(.id = "subjectid")
}

#' Extract group level covariances from the samples
#'
#' This function taks a pmwgs sampler object and extracts the group level
#' covariance samples for the specified parameter estimates.
#'
#' @inheritParams extract_tmu
#'
#' @return A 3D array with (M, M, N) where M is parameters and N is samples
#'
#' @import checkmate
#' @export
extract_cov <- function(sampler,
                        par_names = sampler$par_names,
                        filter = unique(sampler$samples$stage)) {
  samples <- sampler$samples
  filter <- consolidate_filter(sampler, filter)

  assert_subset(par_names, sampler$par_names)

  samples$theta_sig[par_names, par_names, filter]
}


#' Extract parameters from the samples
#'
#' This function taks a pmwgs sampler object and extracts both the group level
#' and the individual subject level samples for the specified parameter
#' estimates. The resulting tibble will bein long format
#'
#' @param sampler The pmwgs sampler object
#' @param par_names The names of the parameters to extract - defaults to all
#'   parameters estimated.
#' @inheritParams pmwg::as_mcmc
#'
#' @return A tibble with the parameter samples and a subjectid column
#'
#' @import dplyr
#' @export
extract_parameters <- function(sampler,
                               par_names = sampler$par_names,
                               filter = unique(sampler$samples$stage)) {
  tmu <- extract_tmu(sampler, par_names, filter) %>%
    mutate(subjectid = "theta_mu")

  extract_alpha(sampler, par_names, filter) %>%
    bind_rows(tmu) %>%
    tidyr::pivot_longer(cols = -c(subjectid, sampleid, stageid),
                        names_to = "parameter")
}


#' Return summary values for parameters (from extract_parameters)
#'
#' @param pars The result of extract_parameters, tibble containing samples for
#'   the selected parameter estimates.
#' @param tform Any transformation necessary to the data prior to summarising
#' @param sfunc A summary function (one value), default is median
#'
#' @return A tibble containing the medians of the samples for each subject
#'
#' @import dplyr
#' @importFrom rlang .data
#' @export
get_summary <- function(pars, tform = base::identity, sfunc = stats::median) {
  pars %>%
    mutate(value = tform(.data$value)) %>%
    group_by(.data$subjectid, .data$parameter) %>%
    summarise(value = sfunc(value))
}


#' Reshape median architecture estimates for plotting
#'
#' @param medians
#'
#' @return A tibble augmented ready for plotting
#' @export
arch_medians <- function(medians) {
  medians %>%
    group_by(subjectid) %>%
    mutate(rel_val = value / sum(value)) %>%
    mutate(parameter = str_remove(parameter, "alpha_")) %>%
    mutate(subjectid = case_when(
      subjectid == "theta_mu" ~ "Group",
      TRUE ~ str_pad(subjectid, 2, pad = "0")
    ))
}


#' Takes a sample tibble and rearranges it to a long df with just drift rates.
#'
#' @param sample_df A dataframe with columns for each parameter and sampleid
#'
#' @return A long tibble containing the drift rate type and values
#'
#' @import dplyr
#' @importFrom readr parse_factor
#' @importFrom rlang .data
#' @export
get_drifts <- function(sample_df) {
  sample_df %>%
  select(c(.data$sampleid, starts_with("v_"))) %>%
  pivot_longer(
    -.data$sampleid,
    names_to = c("drift", "response", "attribute", "salience"),
    names_transform = list(
      response = ~ parse_factor(.x, levels = c("acc", "rej")),
      attribute = ~ parse_factor(.x, levels = c("p", "r")),
      salience = ~ parse_factor(.x, levels = c("H", "L", "D"))
    ),
    names_sep = "_"
  ) %>%
  select(-.data$drift) %>%
  rename(drift = .data$value) %>%
  mutate(
    response = fct_recode(.data$response, Accept = "acc", Reject = "rej"),
    attribute = fct_recode(.data$attribute, Price = "p", Rating = "r")
  )
}

#' Reorder architectures and subjects of medians for plotting
#'
#' From a arch_medians tibble, reorder the subjectID's and architectures
#' in order for them to be presented in a plot nicely.
#'
#' @param model_medians
#'
#' @return ordered model medians
#'
#' @export
order_arch_medians <- function(model_medians, arch_order = NULL) {
  if (is.null(arch_order)) {
    arch_order <- model_medians %>%
      group_by(parameter) %>%
      summarise(mean_val = mean(rel_val)) %>%
      arrange(mean_val) %>%
      pull(parameter)
  }

  most_common_arch <- arch_order[length(arch_order)]

  subject_order <- model_medians %>%
    filter(parameter == most_common_arch) %>%
    arrange(desc(rel_val)) %>%
    pull(subjectid)

  model_medians %>%
    mutate(subjectid = factor(subjectid, subject_order)) %>%
    filter(subjectid != "Group") %>%
    mutate(parameter = factor(parameter, arch_order))
}
