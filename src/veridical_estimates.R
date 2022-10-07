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
    get_summary() %>%
    group_by(subjectid) %>%
    mutate(rel_val = value / sum(value))
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

model_medians <- bind_rows(model_medians, .id = "source") %>%
  mutate(source = factor(source, levels = c("accept", "reject"))) %>%
  mutate(subjectid = case_when(
    subjectid == "theta_mu" ~ "Group",
    TRUE ~ str_pad(subjectid, 2, pad = "0")
  ))

model_plot <- function(medians) {
  Par_order <- c("IST", "CB", "FPP", "MW", "IEX")
  medians %>%
    filter(subjectid != "Group") %>%
    mutate(Parameter = factor(Parameter, Par_order)) %>%
    ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    scale_fill_watercolour() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(y = "Relative Evidence") +
    scale_y_continuous(labels = NULL, breaks = NULL)
}

design <- "
  11111
  22333
"
accept_plot <- model_medians %>% filter(source == "accept")  %>% model_plot
reject_plot <- model_medians %>% filter(source == "reject")  %>% model_plot
accept_plot + reject_plot + guide_area() + plot_layout(design = design, guides = "collect")

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

par_medians <- bind_rows(par_medians, .id = "source") %>% mutate(value = log(value))

par_colours <- c("t0" = "grey",
                 "A" = "grey", "b_acc" = "#73842E", "b_rej" = "#D0781C",
                 "v_acc_p_H" = "#569F72", "v_acc_p_L" = "#407755", "v_acc_p_D" = "#2B5039",
                 "v_rej_p_H" = "#F29F40", "v_rej_p_L" = "#EF8C1A", "v_rej_p_D" = "#BF6C0D",
                 "v_acc_r_H" = "#85C0FF", "v_acc_r_L" = "#47A0FF", "v_acc_r_D" = "#0A81FF",
                 "v_rej_r_H" = "#BBA9A0", "v_rej_r_L" = "#A2887C", "v_rej_r_D" = "#83685D")

par_medians %>%
  mutate(source = recode(source, accept = "IEX", reject = "IST")) %>%
  mutate(colour = par_colours[Parameter]) %>%
  ggplot(aes(x = Parameter, y = exp(value), fill = colour, color = source)) +
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
