library(tidyverse)
library(pmwg)
library(mcce)

pref_file <- here::here("data", "output", "PrefDCE_2506730.rcgbcm_Estimation5Model.RData")

pref_samples <- get_samples(pref_file)

set.seed(2468)

participant <- sample(pref_samples$subjects, 1)
participant_data <- pref_samples$data %>% filter(subject == participant)
participant_random_effects <- pref_samples %>%
  as_mcmc(selection = "alpha", filter = "sample") %>%
  pluck(participant) %>%
  as_tibble

means <- participant_random_effects %>%
  summarise(across(everything(), mean))

winning_arch <- means %>%
  select(starts_with('alpha')) %>%
  rowwise() %>%
  mutate(max_alpha = names(.)[which.max(across(everything()))]) %>%
  pull(max_alpha) %>%
  str_split("_") %>%
  unlist %>%
  `[`(2)

nonarch_parameters <- means %>%
  select(!starts_with("alpha")) %>%
  unlist


ll_func <- select_ll(winning_arch)

ll_vals <- sapply(names(nonarch_parameters), FUN = function(par) {
  par_val <- nonarch_parameters[par]
  alt_vals <- seq(par_val - 1, par_val + 1, 0.05)
  x <- nonarch_parameters
  tibble(
    par_val = alt_vals,
    ll = sapply(alt_vals, FUN = function(sub) {
      x[par] <- sub
      y <- transform_pars(x)
      trial_ll <- model_wrapper(y, participant_data, ll_func)
      sum(log(pmax(trial_ll, 1e-10)))
    })
  )
}, simplify=FALSE, USE.NAMES=TRUE) %>%
bind_rows(.id="parameter")

par_tbl <- tibble(parameter=names(nonarch_parameters), value=nonarch_parameters)

ll_vals %>%
  ggplot(aes(x = par_val, y = ll)) +
  geom_line() +
  geom_point() +
  geom_vline(data=par_tbl, aes(xintercept=value)) +
  facet_wrap(~ parameter, scales="free")

group_level <- pref_samples %>%
  as_mcmc(selection = "theta_mu", filter = "sample") %>%
  as_tibble %>%
  summarise(across(everything(), mean)) %>%
  select(-starts_with("alpha")) %>%
  pivot_longer(cols=everything()) %>%
  print
