library(dplyr)
library(mcce)
library(tidyr)
library(ggplot2)
library(readr)
library(watercolours)

fac_cols <- c("cell_name", "sonaID", "group", "hand")
sft_columns <- c("rt", "group", "hand", "sonaID", "accept", "cell_name")
unused_cols <- c("view_history", "trial_type", "trial_index", "time_elapsed", "internal_node_id", "stimulus", "type", "staircase_point", "staircase_result", "button_pressed", "responses", "response")

datafile <- here::here("data", "input", "PrefSFT", "results_20201116155132_sem2.txt")
raw <- get_JATOS_data(datafile)
extracted <- extract_data(raw)

raw_sft <- extracted %>%
  filter(type == "sft") %>%
  mutate_at(fac_cols, as.factor) %>%
  select(-all_of(unused_cols))

# Remove duplicate 100591 trials (715/15) - somehow they ran part of the task
# twice
nodupes_sft <- raw_sft  %>%
  filter(sonaID != 100591 | hand == "right_accept")

condition_split <- function(df) {
  df %>%
    group_by(group, hand) %>%
    summarise(condition_count = n_distinct(sonaID))
}

# Original split by conditions
condition_split(nodupes_sft)

# Remove participants with < 3 blocks worth of data (540 trials)
too_few_trials <- nodupes_sft %>%
  group_by(sonaID) %>%
  add_count() %>%
  filter(n < 540)

drop_subjects <- too_few_trials %>%
  distinct(sonaID) %>%
  pull(sonaID)

enough_data_sft <- nodupes_sft %>%
  filter(!sonaID %in% drop_subjects)

# Drop NA in rt
no_na_sft <- enough_data_sft %>%
  mutate(rt = as.numeric(rt)) %>%
  filter(!is.na(rt))

# Look at consistency by RT

consistent_sft <- no_na_sft %>%
  mutate(should_accept = cell_name %in% c("HH", "HL", "LH", "LL")) %>%
  mutate(did_accept = accept == "true") %>%
  mutate(
    consistent = case_when(
      should_accept & did_accept ~ TRUE,
      !should_accept & !did_accept ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  select(-c(should_accept, did_accept))

short_rts <- consistent_sft %>%
  mutate(too_short = rt < 350, too_short_old = rt < 200) %>%
  group_by(sonaID) %>%
  summarise(
    num_short = sum(too_short),
    percent_short_350 = sum(too_short) / n(),
    percent_short_200 = sum(too_short_old) / n()
  )

# Drop subjects with more than 20% short RT's (< 350ms)
drop_subjects <- short_rts %>%
  filter(percent_short_350 > 0.2) %>%
  pull(sonaID)

enough_long_sft <- consistent_sft %>%
  filter(!sonaID %in% drop_subjects)

# Drop short RT trials
short_removed_sft <- enough_long_sft %>%
  filter(rt >= 350)

# Drop long RT trials
long_removed_sft <- short_removed_sft %>%
  filter(rt < 10000)

# Both implicit and explicit has a min threshold point of 0.2 and a max of 1

# From the javascript code, take a point on the line between 0.2 and 1.0,
# which is the scale for threshold, and using values from the survey of hotel
# prices transform than value into a price and a rating (inv_rating) that
# conforms to the relationship between price and ratingn in the surveyed
# hotels. lin_rating is the rating if there was a linear correspondence
# between the threshold point and the rating range from surveyed hotels.
#
# invrel takes a rating (as returned by rel function in the inv_rating
# component) and transforms it to a point on the threshold scale that matches
# the linear relationship for rating. (lin_point_match)
#
# rel <- function(point) {
#   min_price = 100, max_price = 450, min_rating = 6.5, max_rating = 9
#   rate_i = 9.73, rate_s = -319.06
#   price = min_price + point * (max_price - min_price)
#   inv_rating = rate_i + rate_s * (1 / price)
#   lin_rating = min_rating + point * (max_rating - min_rating)
#   list(price=price, inv_rating = inv_rating, lin_rating = lin_rating)
# }
#
# invrel <- function(rating) {
#   min_price = 100, max_price = 450, min_rating = 6.5, max_rating = 9
#   rate_i = 9.73, rate_s = -319.06
#   price_match = rate_s / (rating - rate_i)
#   point_match = (price_match - min_price) / (max_price - min_price)
#   lin_point_match = (rating - min_rating) / (max_rating - min_rating)
#   list( pr = price_match, pt = point_match, lpt = lin_point_match )
# }

# Also from the javascript are the scales used to calculate the price and rating
# for high, low and distractor cell types on the threshold point scale.
# So for price, a high salience price would be sampled on the threshold scale
# from a truncated gaussian 0.15 below the threshold value, with sd 0.025 and
# truncated at 0.1 to 0.2 below.
#
# psft.scales = {
#   H: [0.15, 0.025, 0.1, 0.2],
#   L: [0.05, 0.025, 0, 0.1],
#   D: [-0.1, 0.025, -0.05, -0.15],
# }

expl_cols <- c("rt", "sonaID", "price_threshold", "rating_threshold", "response")
explicit <- extracted %>%
  filter(type == "slider") %>%
  mutate(sonaID = factor(sonaID, levels=levels(raw_sft$sonaID))) %>%
  select(all_of(expl_cols))

expl_summ <- explicit %>%
  mutate(response = as.numeric(response)) %>%
  group_by(sonaID) %>%
  summarise(diff = max(response) - min(response))

drop_subjects <- expl_summ %>%
  filter(diff > 0.1) %>%
  pull(sonaID)

explicit_ok <- long_removed_sft %>%
  filter(!sonaID %in% drop_subjects)

impl_cols <- c("accept", "price", "rating", "staircase_point", "staircase_result", "rt", "sonaID", "price_threshold", "rating_threshold")
implicit <- extracted %>%
  filter(type == "staircase") %>%
  mutate(sonaID = factor(sonaID, levels=levels(raw_sft$sonaID))) %>%
  select(all_of(impl_cols))

impl_summ <- implicit %>%
  mutate(staircase_point = as.numeric(staircase_point)) %>%
  filter(staircase_result == "reversal") %>%
  group_by(sonaID) %>%
  slice_tail(n = 8) %>%
  summarise(mean = mean(staircase_point), sd = sd(staircase_point))

drop_subjects <- impl_summ %>%
  filter(sd > 0.025) %>%
    pull(sonaID)

final_sft <- explicit_ok %>%
  filter(!sonaID %in% drop_subjects)


# Consistency by easiest cell varieties
final_sft %>%
  mutate(cell_type = case_when(
    cell_name %in% c("HH", "HL", "LH", "LL") ~ "DoubleTarget",
    cell_name %in% c("HD", "LD") ~ "PriceTarget",
    cell_name %in% c("DH", "DL") ~ "QualityTarget",
    TRUE  ~ "DoubleDistractor")) %>%
  group_by(sonaID, cell_type) %>%
  summarise(avg_consistency = mean(consistent)) %>%
  group_by(sonaID) %>%
  mutate(sort_col = mean(avg_consistency)) %>%
  group_by(sort_col, sonaID) %>%
  mutate(idx = cur_group_id()) %>%
  arrange(idx) %>%
  ggplot(aes(x = idx, y = avg_consistency, colour = cell_type)) +
  geom_point() +
  geom_line(aes(y = sort_col), colour = "black") +
  scale_colour_watercolour(palette="durer") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
  labs(
    title = "Consistency by cell type by participant",
    y = "Percentage of trials consistent with expected response",
    x = "Participants ordered by avg consistency",
    caption = "The line represents an average of all cell_type means",
    )

ggsave(filename = "Consistency_cell_type.png", dpi = 200)

#Final condition split after removing inconsistent participants
condition_split(final_sft)

pref_data <- final_sft %>%
  select(all_of(sft_columns)) %>%
  mutate(rt = as.numeric(rt) / 1000) %>%
  mutate(subject = as.factor(sonaID)) %>%
  mutate(accept = as.numeric(as.logical(accept)) + 1) %>%
  separate(cell_name, into = c("price", "rating"), sep = 1) %>%
  relocate(subject, rt, accept, price, rating) %>%
  mutate(price = factor(price, levels = c("H", "L", "D")),
         rating = factor(rating, levels = c("H", "L", "D"))) %>%
  select(-sonaID) %>%
  drop_na %>%
  mutate(
    v_acc_p = paste0("v_acc_p_", price),
    v_rej_p = paste0("v_rej_p_", price),
    v_acc_r = paste0("v_acc_r_", rating),
    v_rej_r = paste0("v_rej_r_", rating)
  )

saveRDS(pref_data, file = here::here("data", "output", "Pref_threshold_preprocessed.RDS"))
