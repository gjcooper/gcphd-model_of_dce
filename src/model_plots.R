library(readr)
library(watercolours)
library(mcce)
library(dplyr)
library(pmwg)
library(ggplot2)
library(tidyr)
library(forcats)
library(stringr)
library(patchwork)



task <- "Preferential" ## Other option, something like Preferential

if (task == "Veridical") {
  odata <- "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData"
} else {
  odata <- "PrefDCE_2506730.rcgbcm_Estimation5Model.RData"
}

sample_df <- get_samples(here::here("data", "output", odata)) %>%
  extract_parameters(filter = "sample")

#Plot drift rates at group level for attribute x choice (2 x 2) for each salience level.
sample_df %>%
  filter(subjectid == "theta_mu") %>%
  get_drifts %>%
  ggplot(mapping = aes(x = salience, y = drift)) +
  geom_boxplot(aes(fill = salience)) +
  facet_grid(vars(response), vars(attribute)) +
  scale_fill_watercolour()

model_medians <- sample_df %>%
  get_summary() %>%
  group_by(subjectid) %>%
  mutate(rel_val = value / sum(value)) %>%
  mutate(Parameter = as.character(Parameter)) %>%
  mutate(subjectid = case_when(
    subjectid == "theta_mu" ~ "Group",
    TRUE ~ str_pad(subjectid, 2, pad = "0")
  ))

Par_order <- model_medians %>%
  group_by(Parameter) %>%
  summarise(median_val = median(rel_val)) %>%
  arrange(median_val) %>%
  pull(Parameter)

subject_order <- model_medians %>%
  filter(subjectid != "Group") %>%
  filter(Parameter == Par_order[5]) %>%
  group_by(subjectid) %>%
  summarise(maxRel = max(rel_val)) %>%
  group_by(subjectid) %>%
  summarise(meanMax = mean(maxRel)) %>%
  arrange(desc(meanMax)) %>%
  pull(subjectid)



subject_arch <- model_medians %>%
  filter(subjectid != "Group") %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    labs(y = "Relative Evidence") +
    scale_fill_watercolour() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank()) +
    scale_y_continuous(labels = NULL, breaks = NULL) +
    ggtitle("Indidividual Subjects")

group_arch <- model_medians %>%
  filter(subjectid == "Group") %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    scale_fill_watercolour() +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_y_continuous(labels = NULL, breaks = NULL) +
    ggtitle("Group Level")

subject_arch + group_arch +
  plot_layout(widths = c(6, 1), guides = "collect") +
  plot_annotation(
    title = "Relative Evidence for each of the 5 architectures",
    subtitle = "Each row corresponds to the generating architecture",
    caption = "The area of each stacked bar is the relative proportion of each of the five possible modelled architectures as their dirichlet process"
    )

par_medians <- sample_df %>%
  select(-starts_with("alpha")) %>%
  get_summary() %>%
  mutate(value = log(value))

par_medians %>%
  ggplot(aes(y = value, x = Parameter)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(alpha = 0.2) +
  geom_point(data = par_medians %>% filter(subjectid == "theta_mu"), size=5, shape=4)
