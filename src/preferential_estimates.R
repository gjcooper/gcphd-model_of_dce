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

# Load all the samples
pref_file <- here::here("data", "output", "PrefDCE_2506730.rcgbcm_Estimation5Model.RData")

pref_samples <- get_samples(pref_file)


pdf(file = here::here("results", "Preferential", paste0("theta_mu_trace_", Sys.Date(), ".pdf")), width = 14.1, height = 7.53)
for (par in pref_samples$par_names) {
  g <- pref_samples %>%
    as_mcmc %>%
    data.frame %>%
    tibble %>%
    mutate(sample_id = row_number()) %>%
    mutate(stage = pref_samples$samples$stage) %>%
    pivot_longer(cols = -c(sample_id, stage), names_to = "parameter") %>%
    filter(parameter == par) %>%
    ggplot(aes(x = sample_id, y = value, colour = stage)) +
    geom_line() +
    scale_colour_watercolour() +
    labs(title = par) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank())
    print(g)
}
dev.off()

random_effects <- pref_samples %>%
  as_mcmc(selection = "alpha") %>%
  lapply(FUN = function(x) {
      x %>% data.frame() %>% tibble()
  }) %>%
  bind_rows(.id = "subjectid") %>%
  tibble()

for (par in pref_samples$par_names) {
  print(paste("Generating pdf for", par))
  pdf(file = here::here("results", "Preferential", paste0(par, "_randeff_trace_", Sys.Date(), ".pdf")), width = 14.1, height = 7.53)
  all_subjs <- unique(random_effects$subjectid)
  for (subj_arr in split(all_subjs, ceiling(seq_along(all_subjs) / 4))) {
    g <- random_effects %>%
      filter(subjectid %in% subj_arr) %>%
      group_by(subjectid) %>%
      mutate(sample_id = row_number()) %>%
      mutate(stage = pref_samples$samples$stage) %>%
      ungroup() %>%
      pivot_longer(cols = -c(subjectid, sample_id, stage), names_to = "parameter") %>%
      filter(parameter == par) %>%
      ggplot(aes(x = sample_id, y = value, colour = stage)) +
      geom_line() +
      scale_colour_watercolour() +
      facet_wrap(~ subjectid, nrow=2, ncol=2) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank())
    print(g)
  }
  dev.off()
}

model_medians <- pref_samples %>%
  extract_parameters(str_subset(.$par_names, "alpha")) %>%
  get_medians() %>%
  group_by(subjectid) %>%
  mutate(rel_val = value / sum(value)) %>%
  mutate(subjectid = case_when(
    subjectid == "theta_mu" ~ "Group",
    TRUE ~ str_pad(subjectid, 2, pad = "0")
  ))

subject_order <- model_medians %>%
  filter(Parameter == "MW") %>%
  arrange(desc(rel_val)) %>%
  pull(subjectid)

Par_order <- model_medians %>%
  group_by(Parameter) %>%
  summarise(mean_val = mean(rel_val)) %>%
  arrange(mean_val) %>%
  pull(Parameter)

model_plot <- function(medians) {
  Par_order <- c("FPP", "IST", "IEX", "CB", "MW")
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

model_medians %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(subjectid = fct_relevel(subjectid, "Group", after=Inf)) %>%
  filter(subjectid != "Group") %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  model_plot

ggsave(
  filename = here::here(
    "results",
    "Preferential",
    paste0("EstimatedArch_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)

par_medians <- pref_samples %>%
  extract_parameters(str_subset(.$par_names, "alpha", negate = TRUE)) %>%
  filter(subjectid != 'theta_mu') %>%
  get_medians(alpha = FALSE) %>%
  mutate(value = log(value))

group_pars <- pref_samples %>%
  extract_parameters(str_subset(.$par_names, "alpha", negate = TRUE)) %>%
  filter(subjectid == 'theta_mu') %>%
  get_medians(alpha = FALSE) %>%
  mutate(value = log(value))


par_colours <- c("t0" = "grey",
                 "A" = "grey", "b_acc" = "#73842E", "b_rej" = "#D0781C",
                 "v_acc_p_H" = "#569F72", "v_acc_p_L" = "#407755", "v_acc_p_D" = "#2B5039",
                 "v_rej_p_H" = "#F29F40", "v_rej_p_L" = "#EF8C1A", "v_rej_p_D" = "#BF6C0D",
                 "v_acc_r_H" = "#85C0FF", "v_acc_r_L" = "#47A0FF", "v_acc_r_D" = "#0A81FF",
                 "v_rej_r_H" = "#BBA9A0", "v_rej_r_L" = "#A2887C", "v_rej_r_D" = "#83685D")

par_medians %>%
  mutate(colour = par_colours[Parameter]) %>%
  ggplot(aes(x = Parameter, y = exp(value), fill = colour)) +
  geom_boxplot() +
  ylim(c(0, 7.5)) +
  scale_fill_identity()

ggsave(
  filename = here::here(
    "results",
    "Preferential",
    paste0("EstimatedRE_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)

group_pars %>%
  mutate(colour = rep(par_colours, n_distinct(.$subjectid))) %>%
  ggplot(aes(x = Parameter, y = exp(value), fill=colour)) +
  geom_col() +
  scale_fill_identity()

ggsave(
  filename = here::here(
    "results",
    "Preferential",
    paste0("EstimatedPars_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)
