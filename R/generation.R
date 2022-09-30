#' Generate data for one architecture for all subjects
#'
#' This function takes some data, median parameter estimates and job variables
#' and outputs a tibble with simulated data for all subjects
#'
#' @param model_data The original data for the task being simulated.
#' @param vars The list of job variables.
#' @param medians The median estimates of parameter values for each participant.
#'
#' @return A tibble with the simulated data values matching the original data.
#'
#' @import dplyr
#' @export
generate_data <- function(model_data, vars) {
  medians <- readRDS(here::here("data", "output", vars$median_file))
  subjects <- unique(model_data$subject)
  sample_func <- select_ll(vars$test_model, sample = TRUE)
  pb <- txtProgressBar(min = 0, max = length(subjects), initial = 0, style = 3)
  test_data <- lapply(subjects, FUN = function(subjectid) {
    s_idx <- match(subjectid, subjects)
    subject_data <- model_data %>% filter(subject == subjectid)
    pars <- medians[s_idx, ] %>% select(-subjectid)
    sim_data <- rmodel_wrapper(pars,
                               subject_data,
                               sample_func,
                               contaminant_prob = vars$p_contam,
                               min_rt = vars$min_rt,
                               max_rt = vars$max_rt)
    setTxtProgressBar(pb, s_idx)
    sim_data
})
  close(pb)
  test_data <- test_data %>%
    bind_rows()
}
