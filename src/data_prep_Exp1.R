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

raw_data <- read_expyriment_data(here::here("data", "input", "Task1"), "S*")
reformatted_data <- raw_data %>%
  filter(BlockName != "Practice Block") %>%
  rename(
    PriceRatingOrder = `b'PriceRatingOrder' `,
    ResponseCounterbalancing = `b'ResponseCounterbalancing' `,
    AcceptRejectFocus = `b'AcceptRejectFocus' `
  ) %>%
  mutate(
    PriceRatingOrder = fct_relabel(PriceRatingOrder, cfix),
    ResponseCounterbalancing = fct_relabel(ResponseCounterbalancing, cfix),
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
    acceptRight = ResponseCounterbalancing == "RIGHT (/)",
    RespRight = Button == "RIGHT"
  ) %>%
  rename(Price = PriceSalience, Rating = RatingSalience) %>%
  mutate(Accept = case_when(
    acceptRight & RespRight ~ TRUE,
    !acceptRight & !RespRight ~ TRUE,
    TRUE ~ FALSE)
  ) %>%
  select(
    -c(
      ResponseCounterbalancing, AcceptRejectFocus,
      price_match, rating_match,
      acceptRight, RespRight
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

filtered_by_trial_category <- responded_trials %>%
  group_by(subject_id, trial_cat) %>%
  mutate(pc_correct = mean(Correct)) %>%
  group_by(subject_id) %>%
  filter(min(pc_correct) >= 0.8) %>%
  filter(acceptAND)

filtered_by_trial_category

filtered_by_global_correct <- responded_trials %>%
  group_by(subject_id) %>%
  filter(mean(Correct) >= 0.8) %>%
  filter(acceptAND)

filtered_by_global_correct

responded_trials %>%
  group_by(subject_id) %>%
  mutate(pc_correct = mean(Correct)) %>%
  ungroup() %>%
  mutate(subject_id = fct_reorder(subject_id, pc_correct)) %>%
  ggplot(map = aes(x=subject_id, y=pc_correct, colour=acceptAND)) +
  geom_point() +
  scale_x_discrete() +
  geom_hline(yintercept=0.8) + geom_hline(yintercept=0.75, linetype = "dashed")

responded_trials %>%
  group_by(subject_id, trial_cat) %>%
  mutate(pc_correct = mean(Correct)) %>%
  ungroup() %>%
  group_by(subject_id) %>%
  mutate(min_pc_correct = min(pc_correct)) %>%
  ungroup() %>%
  mutate(subject_id = fct_reorder(subject_id, min_pc_correct)) %>%
  ggplot(map = aes(x=subject_id, y=pc_correct, colour=acceptAND, shape=trial_cat)) +
  geom_point(size = 4) +
  scale_x_discrete() +
  geom_hline(yintercept=0.8) + geom_hline(yintercept=0.75, linetype = "dashed")


# Only accept trials
# Create simplifed data for modelling with rtdists
# add drift parameter names
preprocessed_data <- filtered_by_trial_category %>%
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

saveRDS(preprocessed_data, file=here::here("data", "output", "Task1_preprocessed.RDS"))
prev2 <- readRDS(file=here::here("data", "output", "old_T1_prepr.RDS"))

t1 <- prev2 %>% filter(acceptAND) %>% select(-c(acceptRight, RespRight))
t2 <- filtered_by_trial_category %>% select(-c(BlockName, pc_correct))
all(t1 == t2)
