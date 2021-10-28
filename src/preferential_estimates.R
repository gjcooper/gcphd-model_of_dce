library(mcce)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(paletti)
library(readr)
library(forcats)

frankwebb_cols <- read_lines(file = "palette.txt")
viz_palette(frankwebb_cols, "Frank Webb palette")
fill_palette <- get_scale_fill(get_pal(frankwebb_cols))
col_palette <- get_scale_colour(get_pal(frankwebb_cols))

# Load all the samples
pref_file <- here::here("data", "output", "PrefDCE_1896523.rcgbcm_Estimation5Model.RData")

pref_samples <- get_samples(pref_file)

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
  filter(Parameter == "CB") %>%
  arrange(desc(rel_val)) %>%
  pull(subjectid)

Par_order <- model_medians %>%
  group_by(Parameter) %>%
  summarise(mean_val = mean(rel_val)) %>%
  arrange(mean_val) %>%
  pull(Parameter)

model_plot <- function(medians) {
  Par_order <- c("IST", "CB", "FPP", "MW", "IEX")
  medians %>%
    filter(subjectid != "Group") %>%
    mutate(Parameter = factor(Parameter, Par_order)) %>%
    ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    fill_palette() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(y = "Relative Evidence") +
    scale_y_continuous(labels = NULL, breaks = NULL)
}

model_medians %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(subjectid = fct_relevel(subjectid, "Group", after=Inf)) %>%
#  mutate(subjectid = factor(as.numeric(subjectid), levels = 1:27, labels = c(letters, "Group"))) %>%
  filter(subjectid != "Group") %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  model_plot

par_medians <- pref_samples %>%
  extract_parameters(str_subset(.$par_names, "alpha", negate = TRUE)) %>%
  get_medians(alpha = FALSE) %>%
  mutate(value = log(value))

par_medians %>%
  ggplot(aes(x = Parameter, y = exp(value))) +
  geom_boxplot() +
  ylim(c(0, 7.5))
