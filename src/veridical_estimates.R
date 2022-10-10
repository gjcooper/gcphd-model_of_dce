library(mcce)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(watercolours)
library(readr)

# Load all the samples
accept_file <- here::here("data", "output", "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData")
reject_file <- here::here("data", "output", "NumericVDCE_2117450.rcgbcm_RejectEstimation5Model.RData")

accept_samples <- get_samples(accept_file)
reject_samples <- get_samples(reject_file)

samples <- list(accept = accept_samples, reject = reject_samples)

model_medians <- sapply(samples, function(x) {
  extract_parameters(x, str_subset(x$par_names, "alpha")) %>%
    get_summary(tform = exp) %>%
    arch_medians()
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

model_medians <- bind_rows(model_medians, .id = "source") %>%
  mutate(source = factor(source, levels = c("accept", "reject")))

model_plot <- function(medians) {
  Par_order <- c("IST", "CB", "FPP", "MW", "IEX")
  medians %>%
    filter(subjectid != "Group") %>%
    mutate(parameter = factor(parameter, Par_order)) %>%
    ggplot(aes(x = subjectid, y = rel_val, fill = parameter)) +
    geom_col() +
    scale_fill_watercolour() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(y = "Relative Evidence") +
    scale_y_continuous(labels = NULL, breaks = NULL)
}

model_medians %>% mutate(parameter = factor(parameter, Par_order)) %>% arch_plot + scale_fill_watercolour() + facet_grid(rows = vars(source), scales="free", space = "free", shrink=TRUE, drop=TRUE) + coord_flip() + theme(axis.title.x = element_text())

ggsave(
  filename = here::here(
    "results",
    "Veridical",
    paste0("EstimatedArch_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)


par_medians <- sapply(samples, function(x) {
  extract_parameters(x, str_subset(x$par_names, "alpha", negate = TRUE)) %>%
    filter(subjectid != "theta_mu") %>%
    get_summary()
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

par_medians <- bind_rows(par_medians, .id = "source")

par_colours <- c("t0" = "grey",
                 "A" = "grey", "b_acc" = "#73842E", "b_rej" = "#D0781C",
                 "v_acc_p_H" = "#569F72", "v_acc_p_L" = "#407755", "v_acc_p_D" = "#2B5039",
                 "v_rej_p_H" = "#F29F40", "v_rej_p_L" = "#EF8C1A", "v_rej_p_D" = "#BF6C0D",
                 "v_acc_r_H" = "#85C0FF", "v_acc_r_L" = "#47A0FF", "v_acc_r_D" = "#0A81FF",
                 "v_rej_r_H" = "#BBA9A0", "v_rej_r_L" = "#A2887C", "v_rej_r_D" = "#83685D")

par_medians %>%
  mutate(source = recode(source, accept = "IEX", reject = "IST")) %>%
  mutate(colour = par_colours[parameter]) %>%
  ggplot(aes(x = parameter, y = exp(value), fill = colour, color = source)) +
  scale_colour_watercolour() +
  scale_fill_identity() +
  geom_boxplot()

ggsave(
  filename = here::here(
    "results",
    "Veridical",
    paste0("EstimatedPars_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)
