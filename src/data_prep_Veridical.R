library(tidyr)
library(dplyr)
library(forcats)
library(ggplot2)
library(mcce)

# Global options
task <- "Task1"  # Other option is Task2 for Symbolic task
condition <- "Accept"  # Other option is Reject

cfix <- function(x, end = 1) {
  substr(x, 3, nchar(x) - end)
}

short_codes <- c(H = "High", L = "Low", D = "OutOfBounds")

# Read in data,
# Remove practice blocks
# Clean column names/values
# Remove unused columns (TrialID, Price, Rating, PriceRatingOrder
# Turn Correct column into logical
# Get double target, single target and double distractor trials
# Clean based on minimum % correct in all of four trial categories,
# Drop rows with no response (NA values)
# 

raw_data <- read_expyriment_data(here::here("data", "input", task), "S*")

reformat <- function(x) {
  factor_cols <- c("PriceRatingOrder", "ResponseCounterbalancing",
                   "AcceptRejectFocus", "GreyedItemDisplay")
y <- raw_data %>%
  tibble() %>%
  filter(BlockName != "Practice Block") %>%
  rename_with(.fn = cfix, .cols = starts_with("b'"), end=2) %>%
  mutate(across(any_of(factor_cols), .fns = cfix)) %>%
  mutate(
    PriceSalience = fct_recode(PriceSalience, !!!short_codes),
    RatingSalience = fct_recode(RatingSalience, !!!short_codes)
  ) %>%
  select(-c(TrialID, Price, Rating, PriceRatingOrder)) %>%
  mutate(Correct = as.logical(Correct)) %>%
  mutate(
    price_match = PriceSalience %in% c("H", "L"),
    rating_match = RatingSalience %in% c("H", "L")
  ) %>%
  mutate(trial_cat = case_when(
    price_match & rating_match ~ "both",
    price_match & !rating_match ~ "psing",
    rating_match & !price_match ~ "rsing",
    !(price_match | rating_match) ~ "neither"
  )) %>%
  rename(Price = PriceSalience, Rating = RatingSalience) %>%
  mutate(subject_id = factor(subject_id)) %>%
  mutate(acceptAND = AcceptRejectFocus == "Accept")
}

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


filtered_by_trial_category <- responded_trials %>%
    group_categories %>%
    mutate(pc_correct = mean(Correct)) %>%
    group_outer %>%
    filter(min(pc_correct) >= 0.8) %>%
    filter(acceptAND)

filtered_by_global_correct <- responded_trials %>%
  group_outer %>%
  filter(mean(Correct) >= 0.8) %>%
  filter(acceptAND)

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

model_data_format <- function(cleaned_data) {
  # Only accept trials
  # Create simplifed data for modelling with rtdists
  # add drift parameter names
  cleaned_data %>%
    transmute(
      rt = RT / 1000,
      subject = subject_id,
      accept = as.numeric(Accept) + 1,
      price = Price,
      rating = Rating) %>%
    mutate(
      v_acc_p = paste0("v_acc_p_", price),
      v_rej_p = paste0("v_rej_p_", price),
      v_acc_r = paste0("v_acc_r_", rating),
      v_rej_r = paste0("v_rej_r_", rating)
    )
}

if (task == "Task1") {
  filtered_by_trial_category %>%
    model_data_format %>%
    saveRDS(file = here::here("data", "output", "Task1_preprocessed.RDS"))
} else {
  filtered_by_trial_category %>%
    filter(Display == "Absent") %>%
    model_data_format %>%
    saveRDS(file = here::here("data", "output", "Task2_preprocessed_Absent.RDS"))

  filtered_by_trial_category %>%
    filter(Display == "Greyed") %>%
    model_data_format %>%
    saveRDS(file = here::here("data", "output", "Task2_preprocessed_Greyed.RDS"))
}
