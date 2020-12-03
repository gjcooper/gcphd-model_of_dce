require(tidyverse)
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
task1_data <- read_expyriment_data(here::here("data", "input", "Task1"), "S*") %>%
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
  group_by(subject_id, trial_cat) %>%
  mutate(pc_correct = mean(Correct)) %>%
  group_by(subject_id) %>%
  filter(min(pc_correct) >= 0.8) %>%
  drop_na() %>%
  mutate(
    acceptAND = AcceptRejectFocus == "Accept",
    acceptRight = ResponseCounterbalancing == "RIGHT (/)"
  ) %>%
  select(-c(BlockName, ResponseCounterbalancing, AcceptRejectFocus, pc_correct,
            price_match, rating_match)) %>%
  rename(Price = PriceSalience, Rating = RatingSalience) %>%
  mutate(RespRight = Button == "RIGHT") %>%
  mutate(Accept = case_when(
    acceptRight & RespRight ~ TRUE,
    !acceptRight & !RespRight ~ TRUE,
    TRUE ~ FALSE)
  )

# task1_data %>%
#   filter(trial_cat != "both") %>%
#   sample_n(1) %>%
#   head() %>%
#   data.frame()

saveRDS(task1_data, file=here::here("data", "output", "Task1_preprocessed.RDS"))
