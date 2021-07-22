library(mcce)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(pmwg)

original_data <- get_samples(
  here::here(
    "data",
    "output",
    "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData"
  )
)

sample_df <- original_data %>%
  as_mcmc(filter = "sample") %>%
  as_tibble()

sample_df %>%
  select(starts_with("v_")) %>%
  pivot_longer(
    everything(),
    names_to = c("drift", "response", "attribute", "salience"),
    names_transform = list(
      response = ~ readr::parse_factor(.x, levels = c("acc", "rej")),
      attribute = ~ readr::parse_factor(.x, levels = c("p", "r")),
      salience = ~ readr::parse_factor(.x, levels = c("H", "L", "D"))
    ),
    names_sep = "_"
  ) %>%
  select(-drift) %>%
  rename(drift = value) %>%
  mutate(
    response = fct_recode(response, Accept = "acc", Reject = "rej"),
    attribute = fct_recode(attribute, Price = "p", Rating = "r")
  ) %>%
  ggplot(mapping = aes(x = salience, y = drift)) +
  geom_boxplot(aes(fill = salience)) +
  facet_grid(vars(response), vars(attribute))

dev.new()
sample_df %>%
  select(starts_with("v_")) %>%
  pivot_longer(
    everything(),
    names_to = c("drift", "response", "attribute", "salience"),
    names_transform = list(
      response = ~ readr::parse_factor(.x, levels = c("acc", "rej")),
      attribute = ~ readr::parse_factor(.x, levels = c("p", "r")),
      salience = ~ readr::parse_factor(.x, levels = c("H", "L", "D"))
    ),
    names_sep = "_"
  ) %>%
  select(-drift) %>%
  rename(drift = value) %>%
  mutate(
    response = fct_recode(response, Accept = "acc", Reject = "rej"),
    attribute = fct_recode(attribute, Price = "p", Rating = "r")
  ) %>%
  ggplot(mapping = aes(x = salience, y = drift)) +
  geom_boxplot(aes(fill = salience)) +
  facet_grid(vars(attribute), vars(response))
