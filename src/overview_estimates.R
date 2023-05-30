library(mcce)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(watercolours)
library(readr)
library(forcats)
library(tidyr)
library(pmwg)
library(purrr)

datadir <- here::here("data", "output")

pref_samples <- tibble(fname = list.files(datadir, full.names=TRUE, pattern="^Pref.*.RData"),
                       sampler_name = c("sampled", rep("sampler", 13)))

get_theta <- function(fname, sampler_name) {
  print(paste("Extracting from", fname))
  get_samples(fname, final_obj = sampler_name) %>%
    pluck("samples", "theta_mu") %>%
    t %>%
    as_tibble
}

pref_samples <- pref_samples %>%
  mutate(sampler = map2(fname, sampler_name, get_theta)) %>%
  mutate(name = str_extract(fname, "PrefDCE.*")) %>%
  select(-fname)

pref_samples <- pref_samples %>%
  filter(!str_detect(name, "test")) %>%
  filter(!str_detect(name, "CB")) %>%
  filter(!str_detect(name, "staged"))


pref_samples <- pref_samples %>%
  unnest()

pdf(file = here::here("results", "Preferential", "theta_mu_overview.pdf"), width = 14.1, height = 7.53)
for (par in colnames(pref_samples)[2:27]) {
  g <- pref_samples %>%
    select(c(name, one_of(par))) %>%
    group_by(name) %>%
    mutate(x = row_number()) %>%
    ggplot(aes_string(x = "x", y = par)) +
    geom_line() +
    facet_wrap(~ name)
  print(g)
}
dev.off()

pref_samples %>%
  group_by(name) %>%
  mutate(x = row_number()) %>%
  ggplot(aes(x = x, y = A)) +
  geom_line() +
  facet_wrap(~ name)
