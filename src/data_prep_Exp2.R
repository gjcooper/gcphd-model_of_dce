library(tidyr)
library(dplyr)
library(forcats)
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
task2_data <- read_expyriment_data(here::here("data", "input", "Task2"), "S*") %>%
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
  )  %>%
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
  group_by(subject_id, trial_cat, Display) %>%
  mutate(pc_correct = mean(Correct)) %>%
  group_by(subject_id, Display) %>%
  filter(min(pc_correct) >= 0.8) %>%
  drop_na() %>%
  mutate(
    acceptAND = AcceptRejectFocus == "Accept",
  ) %>%
  select(-c(BlockName, AcceptRejectFocus, pc_correct,
            price_match, rating_match)) %>%
  rename(Price = PriceSalience, Rating = RatingSalience) %>%
  mutate(Accept = Button == "RIGHT")

task2_data %>%
  filter(Display == "Absent") %>%
  saveRDS(file=here::here("data", "output", "Task2_preprocessed_Absent.RDS"))

task2_data %>%
  filter(Display == "Greyed") %>%
  saveRDS(file=here::here("data", "output", "Task2_preprocessed_Greyed.RDS"))

