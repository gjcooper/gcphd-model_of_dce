#' Create an output basename from envars list
#'
#' From a list extracted from environment variable for a modelling job
#' extract the necessary element for a basename for output files.
#'
#' @param vars The list of object extracted from the environment
#'
#' @return a character vector with 1 element, the basename
#'
#' @export
get_base_filename <- function(vars) {
  exp <- vars$experiment
  id <- vars$jobid
  tag <- vars$tag

  if (vars$experiment == "SymbolicVDCE") {
    tag <- paste(vars$displaytype, tag, sep = "_")
  }

  if (vars$method == "recovery") {
    tag <- paste(vars$test_model, tag, sep = "_")
  }

  paste(exp, id, tag, sep = "_")
}

#' Extract the data from a sampler object stored in an RData file
#'
#' The sampler object is assumed to be called "sampler", no checks are done
#'
#' @param filename The RData file containing the sampler object
#'
#' @return The data tibble that was used as input to the modelling
#'
#' @export
get_estimation_data <- function(filename) {
  load(here::here("data", "output", filename), ex <- new.env())
  ex$sampler$data %>%
    tibble()
}
