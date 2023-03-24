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

consistent_responders <- long_removed_sft %>%
  group_by(sonaID, cell_name) %>%
  summarise(avg_consistency = mean(consistent)) %>%
  group_by(sonaID) %>%
  slice_min(order_by = avg_consistency, n = 1) %>%
  arrange(desc(avg_consistency)) %>%
  filter(avg_consistency > 0.33) %>%
  pull(sonaID)

final_sft <- long_removed_sft %>%
  filter(sonaID %in% consistent_responders)

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

saveRDS(pref_data, file = here::here("data", "output", "Pref_performance_preprocessed.RDS"))
