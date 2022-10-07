library(mcce)
library(dplyr)
library(ggplot2)

run_all <- function(filename, output) {
  medians <- get_samples(here::here("data", "output", filename)) %>%
    extract_parameters() %>%
    group_by(subjectid, parameter) %>%
    summarise(value = median(value))

  p <- plot_summary(medians, transform = exp) +
    scale_y_continuous(trans = "log2")
  print(p)

  medians %>%
    filter(subjectid != "theta_mu") %>%
    pivot_wider(names_from=parameter, values_from=value) %>%
  saveRDS(file = here::here("data", "output", output))
}

run_all("NumericVDCE_1878182.rcgbcm_Estimation5Model.RData", "median_alpha_exp1.RDS")
run_all("PrefDCE_2506730.rcgbcm_Estimation5Model.RData", "median_alpha_pref.RDS")
