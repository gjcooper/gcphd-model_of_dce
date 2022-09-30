# Modelling is controlled using environment variables. Extracting the variables,
# testing their values and other aspects are controlled here.

#' Known environment variables
known_vars <- c("MCCE_EST_EXP", "VDCE_DISPLAY", "NCPUS", "PBS_JOBID",
                "MCCE_TAG", "MCCE_EXP_DATA",
                "MCCE_MIN_RT", "MCCE_MAX_RT", "MCCE_CONTAM",
                "RANDOM_SEED", "MCCE_MODEL", "MCCE_METHOD",
                "MCCE_ORIG_JOB_DATA", "MCCE_STAGES",
                "MCCE_REC_MODEL", "MCCE_REC_MED", "MCCE_MODEL_FILE",
                "MCCE_REC_DATA")

known_experiments <- c("NumericVDCE", "SymbolicVDCE", "PrefDCE")
known_models <- c("std", "reduced", "reduced_more")
known_methods <- c("test", "model", "continue", "recovery")

msg_parts <- c("System Environment Variable",
               "not defined or unknown value",
               "must be provided and numeric",
               "should be a value between",
               "and")

test_cond <- function(name, condition) {
  if (!condition) {
    stop(paste(msg_parts[1], name, msg_parts[2]))
  }
}

test_num <- function(var, name) {
  if (is.na(var)) {
    stop(paste(msg_parts[1], name, msg_parts[3]))
  }
}

test_range <- function(var, name, lower, upper) {
  if (var < lower || var > upper) {
    stop(paste(msg_parts[1], name, msg_parts[4], lower, msg_parts[5], upper))
  }
}

method_tests <- function(vars) {
  if (vars$method == "recovery") {
    test_cond("MCCE_REC_MODEL", vars$test_model %in% names_ll())
  } else if (vars$method == "continue") {
    test_cond("MCCE_ORIG_JOB_DATA", vars$early_data != "")
    test_cond("MCCE_STAGES", vars$stages_to_run != "")
  }
}


#' Check the system environment variables and extract the relevant pieces for
#' modelling
#'
#' @return A list with the elements extracted from the sys env
#'
#' @export
check_env <- function() {
  envars <- as.list(Sys.getenv(known_vars))
  print_env(envars)

  vars <- exp_env(envars)
  vars <- c(vars, job_env(envars))
  vars <- c(vars, data_env(envars))
  vars <- c(vars, model_env(envars))

  method_tests(vars)
  vars
}

print_env <- function(envars) {
  envars <- t(data.frame(envars))
  message(paste0(c("Currently set environment variables",
                   "___________________________________",
                   paste0(capture.output(envars), collapse = "\n")),
                 collase = "\n"))
}

exp_env <- function(envars) {
  experiment <- envars$MCCE_EST_EXP
  min_rt <- as.numeric(envars$MCCE_MIN_RT)
  max_rt <- as.numeric(envars$MCCE_MAX_RT)
  p_contam <- as.numeric(envars$MCCE_CONTAM)
  displaytype <- envars$VDCE_DISPLAY

  # Tests for envars
  test_cond("MCCE_EST_EXP", experiment %in% known_experiments)
  if (experiment == "SymbolicVDCE") {
    test_cond("VDCE_DISPLAY", displaytype %in% c("Absent", "Greyed"))
  }
  test_num(min_rt, "MCCE_MIN_RT")
  test_num(max_rt, "MCCE_MAX_RT")
  test_num(p_contam, "MCCE_CONTAM")
  test_range(p_contam, "MCCE_CONTAM", 0, 1)
  list(experiment = experiment, min_rt = min_rt, max_rt = max_rt,
       p_contam = p_contam, displaytype = displaytype)
}

rand_str <- function(len) {
  paste(sample(c(strsplit("0123456789", "")[[1]],
                 letters),
               len, replace=TRUE),
        collapse="")
}

job_env <- function(envars) {
  cores <- ifelse(envars$NCPUS == "", 1, as.numeric(envars$NCPUS))
  jobid <- ifelse(envars$PBS_JOBID == "", rand_str(16), envars$PBS_JOBID)
  tag <- ifelse(envars$MCCE_TAG == "", "untagged", envars$MCCE_TAG)
  if (envars$RANDOM_SEED != "") {
    set.seed(as.numeric(envars$RANDOM_SEED))
  }
  list(cores = cores, jobid = jobid, tag = tag)
}

data_env <- function(envars) {
  early_data <- envars$MCCE_ORIG_JOB_DATA
  experimental_data <- envars$MCCE_EXP_DATA
  test_model <- envars$MCCE_REC_MODEL
  median_file <- envars$MCCE_REC_MED
  model_file <- envars$MCCE_MODEL_FILE
  recovery_data <- envars$MCCE_REC_DATA

  list(early_data = early_data,
       experimental_data = experimental_data,
       test_model = test_model,
       median_file = median_file,
       model_file = model_file,
       recovery_data = recovery_data)
}

model_env <- function(envars) {
  model <- envars$MCCE_MODEL
  method <- envars$MCCE_METHOD
  stages_to_run <- envars$MCCE_STAGES

  test_cond("MCCE_MODEL", model %in% known_models)
  test_cond("MCCE_METHOD", method %in% known_methods)

  list(model = model, method = method, stages_to_run = stages_to_run)
}
