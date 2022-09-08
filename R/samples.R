#' Extract group level summary stats from a sampler
#'
#' This function extracts group level summary statistics from a pmwgs object
#'
#' @param sampler A pmwgs sampler containing mcce samples
#' @param func A summary stat function
#'
#' @return A tibble
#'
#' @import dplyr
#' @export
tmu <- function(sampler, func = median) {
  sampler %>%
    pmwg::as_mcmc(filter = "sample") %>%
    as_tibble() %>%
    setNames(sampler$par_names) %>%
    summarise_all(func)
}

#' Extract subject level summary stats from a sampler
#'
#' This function extracts subject level summary statistics from a pmwgs object
#'
#' @param sampler A pmwgs sampler containing mcce samples
#' @param func A summary stat function
#'
#' @return A tibble
#'
#' @import dplyr
#' @export
alph <- function(sampler, func = median) {
  sampler %>%
    pmwg::as_mcmc(selection = "alpha", filter = "sample") %>%
    sapply(function(x) {
             x %>%
               as_tibble() %>%
               setNames(sampler$par_names) %>%
               summarise_all(func)},
           simplify=FALSE) %>%
    bind_rows(.id = "SubjectID")
}
