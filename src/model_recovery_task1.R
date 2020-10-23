require(pmwg)
require(rtdists)
library(MCMCpack)
require(tidyverse)
source(here::here("src", "read_expyriment.R"))
source(here::here("src", "loglike.R"))

load(here::here("data", "output", "reboot_full.RData"), ex <- new.env())
medians <- readRDS(here::here("data", "output", "median_alpha_exp1.RDS"))

model_data <- ex$sampled$data %>%
  tibble()

test_model <- "IST"

sample_func <- rll_funcs[[match(test_model, ll_names)]]

subjects <- unique(model_data$subject)


test_data <- lapply(subjects, FUN = function(subjectid) {
  s_idx <- match(subjectid, subjects)
  subject_data <- model_data %>% filter(subject == subjectid)
  pars <- medians[s_idx, ]

  sample_func(pars, subject_data)
})

test_data <- test_data %>%
  bind_rows()
