require(pmwg)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

tmu <- function(sampled) {
  median_theta_mu <- sampled %>%
    as_mcmc(filter = "sample") %>%
    as_tibble() %>%
    setNames(sampled$par_names) %>%
    summarise_all(median)
  median_theta_mu
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

alph <- function(sampled) {
  median_alpha <- sampled %>%
    as_mcmc(selection = "alpha", filter = "sample") %>%
    convert(sampled$par_names) %>%
    bind_rows()
  median_alpha
}

plot_medians <- function(median_alpha, median_theta_mu) {
  median_alpha %>%
    rownames_to_column("SubjectID") %>%
    pivot_longer(!SubjectID) %>%
    ggplot(mapping = aes(x = name, y = value, color = SubjectID)) +
    geom_point() +
    geom_point(
      data = pivot_longer(median_theta_mu, everything()),
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
  print(plot_medians(alpha_median, theta_median))
  saveRDS(alpha_median, file = here::here("data", "output", output))
}
  
run_all("NumericVDCE_1878182.rcgbcm_Estimation5Model.RData", "median_alpha_exp1.RDS")
run_all("Task2_Absent_1069903.rcgbcm_CorrectedTry1.RData", "median_alpha_exp2_abs.RDS")
run_all("Task2_Greyed_1069904.rcgbcm_CorrectedTry1.RData", "median_alpha_exp2_grey.RDS")
