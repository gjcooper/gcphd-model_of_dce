library(mcce)
library(ggplot2)

run_all <- function(filename, output) {
  load(here::here("data", "output", filename))
  theta_median <- tmu(sampled)
  alpha_median <- alph(sampled)
  p <- plot_summary(alpha_median, theta_median, transform = exp) +
    scale_y_continuous(trans = "log2")
  print(p)
  saveRDS(alpha_median, file = here::here("data", "output", output))
}

run_all("NumericVDCE_1878182.rcgbcm_Estimation5Model.RData", "median_alpha_exp1.RDS")
run_all("PrefDCE_2506730.rcgbcm_Estimation5Model.RData", "median_alpha_pref.RDS")
