library(tidyr)
library(dplyr)
library(forcats)
library(ggplot2)
library(mcce)

# Global options
task <- "Task1"  # Other option is Task2 for Symbolic task
condition <- "Accept"  # Other option is Reject
display <- "Absent"  # Other option is Greyed

raw_data <- read_expyriment_data(here::here("data", "input", task), "S*")

if (task == "Task1") {
  reformatted <- raw_data %>%
    reformat() %>%
    mutate(Accept = case_when(
      (ResponseCounterbalancing == "RIGHT (/)") & (Button == "RIGHT") ~ TRUE,
      !(ResponseCounterbalancing == "RIGHT (/)") & !(Button == "RIGHT") ~ TRUE,
      TRUE ~ FALSE)
  )
} else {
  reformatted <- raw_data %>%
    reformat() %>%
    mutate(Accept = Button == "RIGHT")
}

# Plot NA per subject
reformatted %>%
  filter_all(any_vars(is.na(.))) %>%
  group_by(subject_id) %>%
  summarise(n = n(), accept=all(acceptAND)) %>%
  ggplot(map = aes(x = subject_id, y = n, fill = accept)) + geom_col()

responded_trials <- reformatted %>%
  drop_na()

responded_trials

group_ex1_categories <- function(x) {
  x %>% group_by(subject_id, trial_cat)
}
group_ex1_outer <- function(x) {
  x %>% group_by(subject_id)
}
group_ex2_categories <- function(x) {
  x %>% group_by(subject_id, trial_cat, Display)
}
group_ex2_outer <- function(x) {
  x %>% group_by(subject_id, Display)
}

if (task == "Task1") {
  group_categories <- group_ex1_categories
  group_outer <- group_ex1_outer
} else {
  group_categories <- group_ex2_categories
  group_outer <- group_ex2_outer
}


# Matching actual process used in SFT paper
filtered_by_trial_category <- responded_trials %>%
  filter(RT > 300, RT < quantile(RT, .95)) %>%
  group_categories %>%
  mutate(pc_correct = mean(Correct)) %>%
  group_outer %>%
  filter(min(pc_correct) >= 0.8)


filtered_by_global_correct <- responded_trials %>%
  group_outer %>%
  filter(mean(Correct) >= 0.8)

global_correct_plot <- responded_trials %>%
  group_outer %>%
  mutate(pc_correct = mean(Correct)) %>%
  ungroup() %>%
  mutate(subject_id = fct_reorder(subject_id, pc_correct)) %>%
  ggplot(map = aes(x = subject_id, y = pc_correct, colour = acceptAND)) +
  scale_x_discrete() +
  geom_hline(yintercept = 0.8) +
  geom_hline(yintercept = 0.75, linetype = "dashed")

if (task == "Task1") {
  global_correct_plot + geom_point(size = 3)
} else {
  global_correct_plot + geom_point(aes(shape = Display), size = 3)
}

category_split_plot <- responded_trials %>%
  group_categories %>%
  mutate(pc_correct = mean(Correct)) %>%
  ungroup() %>%
  group_outer %>%
  mutate(min_pc_correct = min(pc_correct)) %>%
  ungroup() %>%
  mutate(subject_id = fct_reorder(subject_id, min_pc_correct)) %>%
  ggplot(map = aes(
    x = subject_id,
    y = pc_correct,
    colour = acceptAND,
    shape = trial_cat)
  ) +
  geom_point(size = 3) +
  scale_x_discrete() +
  geom_hline(yintercept = 0.8) +
  geom_hline(yintercept = 0.75, linetype = "dashed")

if (task == "Task1") {
  category_split_plot
} else {
  category_split_plot + facet_wrap(~Display)
}

model_data_format <- function(cleaned_data, accept_trials = TRUE) {
  # Only accept trials
  # Create simplifed data for modelling with rtdists
  # add drift parameter names
  cleaned_data %>%
    filter(acceptAND == accept_trials) %>%
    transmute(
      rt = RT / 1000,
      subject = subject_id,
      accept = as.numeric(Accept) + 1,
      price = PriceSalience,
      rating = RatingSalience) %>%
    mutate(
      v_acc_p = paste0("v_acc_p_", price),
      v_rej_p = paste0("v_rej_p_", price),
      v_acc_r = paste0("v_acc_r_", rating),
      v_rej_r = paste0("v_rej_r_", rating)
    )
}

if (task == "Task1") {
  savefile <- paste0("Task1_preprocessed_", condition, ".RDS")
} else {
  savefile <- paste0("Task2_preprocessed_", condition, "_", display, ".RDS")
  filtered_by_trial_category <- filtered_by_trial_category %>%
    filter(Display == display)
}

filtered_by_trial_category %>%
  model_data_format(accept_trials = (condition == "Accept")) %>%
  saveRDS(file = here::here("data", "output", savefile))
