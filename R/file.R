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

get_estimation_data <- function(filename) {
  load(here::here("data", "output", filename), ex <- new.env())
  ex$sampler$data %>%
    tibble()
}


