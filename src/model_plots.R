library(readr)
library(paletti)
library(mcce)
library(dplyr)
library(pmwg)
library(ggplot2)
library(tidyr)
library(forcats)
library(stringr)
library(patchwork)


frankwebb_cols <- read_lines(file = "palette.txt")
viz_palette(frankwebb_cols, "Frank Webb palette")
fill_palette <- get_scale_fill(get_pal(frankwebb_cols))
col_palette <- get_scale_colour(get_pal(frankwebb_cols))

task <- "Veridical" ## Other option, something like Preferential

if (task == "Veridical") {
  odata <- "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData"
} else {
  odata <- "PrefDCE_1896523.rcgbcm_Estimation5Model.RData"
}

pmwg_samples <- get_samples(
  here::here(
    "data",
    "output",
    odata
  )
)

sample_df <- pmwg_samples %>%
  as_mcmc(filter = "sample") %>%
  as_tibble()

#Plot drift rates at group level for attribute x choice (2 x 2) for each salience level.
sample_df %>%
  select(starts_with("v_")) %>%
  pivot_longer(
    everything(),
    names_to = c("drift", "response", "attribute", "salience"),
    names_transform = list(
      response = ~ readr::parse_factor(.x, levels = c("acc", "rej")),
      attribute = ~ readr::parse_factor(.x, levels = c("p", "r")),
      salience = ~ readr::parse_factor(.x, levels = c("H", "L", "D"))
    ),
    names_sep = "_"
  ) %>%
  select(-drift) %>%
  rename(drift = value) %>%
  mutate(
    response = fct_recode(response, Accept = "acc", Reject = "rej"),
    attribute = fct_recode(attribute, Price = "p", Rating = "r")
  ) %>%
  ggplot(mapping = aes(x = salience, y = drift)) +
  geom_boxplot(aes(fill = salience)) +
  facet_grid(vars(response), vars(attribute))

model_medians <- pmwg_samples %>%
  extract_parameters(str_subset(.$par_names, "alpha")) %>%
  get_medians() %>%
  group_by(subjectid) %>%
  mutate(rel_val = value / sum(value)) %>%
  mutate(Parameter = as.character(Parameter)) %>%
  mutate(subjectid = case_when(
    subjectid == "theta_mu" ~ "Group",
    TRUE ~ str_pad(subjectid, 2, pad = "0")
  ))

subject_order <- model_medians %>%
  filter(subjectid != "Group") %>%
  group_by(subjectid) %>%
  summarise(maxRel = max(rel_val)) %>%
  group_by(subjectid) %>%
  summarise(meanMax = mean(maxRel)) %>%
  arrange(desc(meanMax)) %>%
  pull(subjectid)


Par_order <- model_medians %>%
  group_by(Parameter) %>%
  summarise(median_val = median(rel_val)) %>%
  arrange(median_val) %>%
  pull(Parameter)

subject_arch <- model_medians %>%
  filter(subjectid != "Group") %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    labs(y = "Relative Evidence") +
    fill_palette() +
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
    fill_palette() +
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

par_medians <- pmwg_samples %>%
  extract_parameters(., str_subset(.$par_names, "alpha", negate = TRUE)) %>%
  get_medians(alpha = FALSE) %>%
  mutate(value = log(value))

par_medians %>%
  ggplot(aes(y = value, x = Parameter)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.3)
