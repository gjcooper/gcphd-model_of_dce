#' Extract alpha parameters (architectures) from the samples
#'
#' This function taks a pmwgs sampler object and extracts both the group level
#' and the individual subject level samples for the alpha parameter estimates.
#'
#' @param sampler The pmwgs sampler object
#'
#' @return A tibble with the alpha parameter samples and a subjectid column
#'
#' @import dplyr
#' @export
extract_alphas <- function(sampler) {
  tmu_alphas <- sampler %>%
    pmwg::as_mcmc(filter = "sample") %>%
    data.frame() %>%
    tibble() %>%
    select(starts_with("alpha")) %>%
    mutate(subjectid = "theta_mu")

  pmwg::as_mcmc(sampler, selection = "alpha", filter = "sample") %>%
    lapply(FUN = function(x) {
      x %>%
        data.frame() %>%
        tibble() %>%
        select(starts_with("alpha"))
    }) %>%
    bind_rows(.id = "subjectid") %>%
    tibble() %>%
    bind_rows(tmu_alphas)
}


#' Return median values for alpha parameters (from extract_alphas)
#'
#' @param alphas The result of extract_alphas, tibble containing samples for
#'   the architecture parameter estimates.
#'
#' @return A tibble containing the medians of the samples for each subject as
#'   well as the relative weight for each architecture.
#'
#' @import dplyr
#' @export
get_medians <- function(alphas) {
  alphas %>%
    group_by(subjectid) %>%
    summarise(across(.fns = stats::median)) %>%
    mutate_at(vars(contains("alpha")), exp) %>%
    tidyr::pivot_longer(-subjectid,
                        names_to = c("drop", "Parameter"),
                        names_sep = "_",
                        names_transform = list(Parameter = as.factor)) %>%
    select(-drop) %>%
    group_by(subjectid) %>%
    mutate(rel_val = value / sum(value))
}
