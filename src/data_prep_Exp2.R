library(tidyr)
library(dplyr)
library(forcats)
library(ggplot2)
devtools::load_all()


cfix <- function(x) {
  substr(x, 3, nchar(x) - 1)
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

raw_data <- read_expyriment_data(here::here("data", "input", "Task2"), "S*")
reformatted_data <- raw_data %>%
  filter(BlockName != "Practice Block") %>%
  rename(
    PriceRatingOrder = `b'PriceRatingOrder' `,
    GreyedItemDisplay = `b'GreyedItemDisplay' `,
    AcceptRejectFocus = `b'AcceptRejectFocus' `
  ) %>%
  mutate(
    PriceRatingOrder = fct_relabel(PriceRatingOrder, cfix),
    GreyedItemDisplay = fct_relabel(GreyedItemDisplay, cfix),
    AcceptRejectFocus = fct_relabel(AcceptRejectFocus, cfix)
  ) %>%
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
  tibble() %>%
  mutate(
    acceptAND = AcceptRejectFocus == "Accept",
    Accept = Button == "RIGHT"
  ) %>%
  rename(Price = PriceSalience, Rating = RatingSalience) %>%
  select(
    -c(
      AcceptRejectFocus,
      price_match, rating_match,
    )
  ) %>%
  mutate(subject_id = factor(subject_id))

# Plot NA per subject
reformatted_data %>%
  filter_all(any_vars(is.na(.))) %>%
  count(subject_id) %>%
  ggplot(map = aes(x = subject_id, y = n)) + geom_col()

responded_trials <- reformatted_data %>%
  drop_na()

responded_trials

trial_cat_filter <- function(trials, display) {
  trials %>%
    filter(Display == display) %>%
    group_by(subject_id, trial_cat) %>%
    mutate(pc_correct = mean(Correct)) %>%
    group_by(subject_id) %>%
    filter(min(pc_correct) >= 0.8) %>%
    filter(acceptAND)
}

filtered_by_trial_category <- responded_trials %>%
  group_by(subject_id, trial_cat, Display) %>%
  mutate(pc_correct = mean(Correct)) %>%
  group_by(subject_id, Display) %>%
  filter(min(pc_correct) >= 0.8) %>%
  filter(acceptAND)

filtered_by_trial_category_absent <- trial_cat_filter(responded_trials, "Absent")
filtered_by_trial_category_greyed <- trial_cat_filter(responded_trials, "Greyed")

filtered_by_trial_category
filtered_by_trial_category_absent
filtered_by_trial_category_greyed

filtered_by_global_correct <- responded_trials %>%
  group_by(subject_id) %>%
  filter(mean(Correct) >= 0.8) %>%
  filter(acceptAND)

global_filter <- function(trials, display) {
  trials %>%
    filter(Display == display) %>%
    group_by(subject_id) %>%
    filter(mean(Correct) >= 0.8) %>%
    filter(acceptAND)
}

filtered_by_global_correct_absent <- global_filter(responded_trials, "Absent")
filtered_by_global_correct_greyed <- global_filter(responded_trials, "Greyed")

filtered_by_global_correct
filtered_by_global_correct_absent
filtered_by_global_correct_greyed

responded_trials %>%
  group_by(subject_id, Display) %>%
  mutate(pc_correct = mean(Correct)) %>%
  ungroup() %>%
  mutate(subject_id = fct_reorder(subject_id, pc_correct)) %>%
  ggplot(map = aes(x=subject_id, y=pc_correct, colour=acceptAND, shape=Display)) +
  geom_point(size=3) +
#  facet_grid(~Display) +
  scale_x_discrete() +
  geom_hline(yintercept=0.8) + geom_hline(yintercept=0.75, linetype = "dashed")

responded_trials %>%
  group_by(subject_id, trial_cat, Display) %>%
  mutate(pc_correct = mean(Correct)) %>%
  ungroup() %>%
  group_by(subject_id, Display) %>%
  mutate(min_pc_correct = min(pc_correct)) %>%
  ungroup() %>%
  mutate(subject_id = fct_reorder(subject_id, min_pc_correct)) %>%
  ggplot(map = aes(x=subject_id, y=pc_correct, colour=acceptAND, shape=trial_cat)) +
  geom_point(size = 4) +
  scale_x_discrete() +
  facet_wrap(~Display) +
  geom_hline(yintercept=0.8) + geom_hline(yintercept=0.75, linetype = "dashed")

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

filtered_by_trial_category_absent %>%
  model_data_format %>%
  saveRDS(file = here::here("data", "output", "Task2_preprocessed_Absent.RDS"))

filtered_by_trial_category_greyed %>%
  model_data_format %>%
  saveRDS(file = here::here("data", "output", "Task2_preprocessed_Greyed.RDS"))

prevAbs <- readRDS(file=here::here("data", "output", "old_T2_prepr_Abs.RDS"))
prevGrey <- readRDS(file=here::here("data", "output", "old_T2_prepr_Gre.RDS"))

compare <- function(prev, now) {
  t1 <- prev %>% filter(acceptAND)
  t2 <- now %>% select(-c(BlockName, pc_correct))
  all(t1 == t2)
}

compare(prevAbs, filtered_by_trial_category_absent)
compare(prevGrey, filtered_by_trial_category_greyed)
