library(mcce)
library(ggplot2)
library(tidyr)
library(dplyr)

median_file <- "median_alpha_pref.RDS"
model_file <- "PrefDCE_1896523.rcgbcm_Estimation5Model.RData"

get_estimation_data <- function(filename) {
  load(here::here("data", "output", filename), ex <- new.env())
  ex$sampled$data %>%
    tibble()
}

medians <- readRDS(here::here("data", "output", median_file))
model_data <- get_estimation_data(model_file)
subjects <- unique(model_data$subject)

all_data <- sapply(names_ll(), simplify=FALSE, USE.NAMES=TRUE, FUN=function(arch) {
  print(paste("Simulating from", arch))
  sample_func <- select_ll(arch, sample=TRUE)
  test_data <- sapply(subjects, USE.NAMES=TRUE, simplify=FALSE, FUN = function(subjectid) {
    s_idx <- match(subjectid, subjects)
    subject_data <- model_data %>% filter(subject == subjectid)
    pars <- medians[s_idx, ]
    t0 <- exp(pars$t0)
    minrt <- min(subject_data$rt)
    rmodel_wrapper(pars, subject_data, sample_func)
  })
  test_data <- test_data %>%
    bind_rows()
})

all_data <- all_data %>% bind_rows(.id="arch")

# Add in original data
all_data <- all_data %>% group_by(subject, arch) %>% mutate(trial = row_number())
model_data <- model_data %>% group_by(subject) %>% mutate(trial = row_number())

model_data %>%
  ungroup() %>%
  right_join(
    all_data,
    by = c("subject", "trial", "price", "rating", "group", "hand"),
    suffix = c(".real", ".sim")
  ) %>%
  select(-c(starts_with("v_"), hand, group)) %>%
  group_by(subject) %>% 
  nest() %>%
  ungroup() %>%
  slice_sample(n=4) %>%
  unnest(cols=c(data)) %>%
  mutate(arch = factor(arch)) %>%
  mutate(cell = paste0(price, rating)) %>%
  filter(rt.sim < 5) %>%
  filter(rt.real < 5) %>%
  ggplot(aes(x = rt.sim, colour = arch)) +
  geom_density() + #(bins = 100, alpha = 0.2) +
  geom_density(mapping=aes(x=rt.real), colour="black") +
  facet_grid(vars(subject), vars(cell), scales="free") 

median_aug <- bind_cols(medians, distinct(model_data, subject)) %>%
  mutate(across(where(is.numeric), exp))


quants_pars <- model_data %>%
  group_by(subject) %>%
  summarise(rt = quantile(rt), q = seq(0,1,0.25)) %>%
  ungroup() %>%
  left_join(median_aug, by = "subject")

# Do any model t0 parameters fail the t0 > minimum rt test?
quants_pars %>%
  filter(q == 0) %>%
  select(c(subject, rt, t0)) %>%
  mutate(failed = t0 > rt) %>%
  filter(failed == TRUE)

# Do any simulated rt's fail the rt < t0 test
sim_quants_pars <- all_data %>%
  group_by(subject, arch) %>%
  summarise(rt = quantile(rt), q = seq(0,1,0.25)) %>%
  ungroup() %>%
  left_join(median_aug, by = "subject")

# Do any model t0 parameters fail the t0 > minimum rt test?
sim_quants_pars %>%
  filter(q == 0) %>%
  select(c(subject, arch, rt, t0)) %>%
  mutate(failed = t0 < rt) %>%
  filter(failed == TRUE)
