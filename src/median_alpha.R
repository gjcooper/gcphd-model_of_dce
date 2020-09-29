require(pmwg)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
load('data/output/reboot_full.RData')

median_theta_mu <- sampled %>%
  as_mcmc(filter = "sample") %>%
  as_tibble() %>%
  setNames(sampled$par_names) %>%
  summarise_all(median)

convert <- function(mcmc_list, par_names) {
  converted_list <- list()
  for (idx in seq_along(mcmc_list)) {
    mcmc_name <- names(mcmc_list)[idx]
    converted_list[[mcmc_name]] <- mcmc_list[[idx]] %>%
      as_tibble() %>%
      setNames(par_names) %>%
      summarise_all(median)
  }
  converted_list
}

median_alpha <- sampled %>%
  as_mcmc(selection = "alpha", filter = "sample") %>%
  convert(sampled$par_names) %>%
  bind_rows()
median_alpha

median_alpha %>%
  rownames_to_column("SubjectID") %>%
  pivot_longer(!SubjectID) %>%
  ggplot(mapping=aes(x=name, y=value, color=SubjectID)) +
    geom_point() +
    geom_point(data = pivot_longer(median_theta_mu, everything()),
               mapping = aes(x=name, y=value),
               size = 3,
               colour = "black") +
    theme(axis.text.x=element_text(angle = -90, hjust = 0))

saveRDS(median_alpha, file = "median_alpha.RDS")
