require(pmwg)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

tmu <- function(sampler) {
  sampler %>%
    as_mcmc(filter = "sample") %>%
    as_tibble() %>%
    setNames(sampler$par_names) %>%
    summarise_all(median)
}

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

alph <- function(sampler) {
  sampler %>%
    as_mcmc(selection = "alpha", filter = "sample") %>%
    convert(sampler$par_names) %>%
    bind_rows()
}

plot_medians <- function(median_alpha, median_theta_mu, transform=identity) {
  median_alpha %>%
    transform %>%
    rownames_to_column("SubjectID") %>%
    pivot_longer(!SubjectID) %>%
    ggplot(mapping = aes(x = name, y = value)) +
    geom_point(colour = "#D0781C") +
    geom_point(
      data = pivot_longer(median_theta_mu %>% transform, everything()),
      mapping = aes(x = name, y = value),
      size = 3,
      colour = "black"
    ) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0))
}

run_all <- function(filename, output) {
  load(here::here("data", "output", filename))
  theta_median <- tmu(sampled)
  alpha_median <- alph(sampled)
  p <- plot_medians(alpha_median, theta_median, transform = exp) +
    scale_y_continuous(trans = "log2")
  print(p)
  saveRDS(alpha_median, file = here::here("data", "output", output))
}
  
run_all("NumericVDCE_1878182.rcgbcm_Estimation5Model.RData", "median_alpha_exp1.RDS")
run_all("PrefDCE_2506730.rcgbcm_Estimation5Model.RData", "median_alpha_pref.RDS")
