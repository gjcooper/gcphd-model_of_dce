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

# Condition split after removing too few trials
condition_split(enough_data_sft)

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

# Look at consistency by duration
consistency_rt_plot <- function(df) {
  bar_col = watercolour$zorn$discrete['chair']
  df %>%
    mutate(
      rt_group = cut(rt, c(seq(0, 700, 50), Inf), include.lowest = TRUE)
    ) %>%
    group_by(rt_group) %>%
    summarise(avg_consistency = mean(consistent), ntrials = n()) %>%
    mutate(
      cumtrials = cumsum(ntrials),
      pctrials = ntrials / sum(ntrials),
      cumpctrials = cumtrials / sum(ntrials),
      txtlabels = paste(
        round(pctrials, 2),
        paste0("(", round(cumpctrials, 2), ")")
      )
    ) %>%
    ggplot(aes(x = rt_group, y = avg_consistency)) +
    geom_bar(stat = "identity", fill = bar_col) +
    geom_text(aes(label = txtlabels, y = avg_consistency + 0.02)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .05)) +
    labs(
      title = "Consistency across different RT windows",
      caption = "Text labels above each bar reflect the proportion of trials the bar represents, and the cumulative proportion from the left, up to and including this bar in brackets",
      x = "Response Time Window (ms)",
      y = "Average Consistency",
    )
}

consistent_sft %>%
  consistency_rt_plot +
  labs(subtitle = "Collapsed across all cells")

ggsave(filename = "Consistency_by_RT.png", dpi = 200)

consistent_sft %>%
  filter(cell_name == "HH") %>%
  consistency_rt_plot +
  labs(subtitle = "For HH cells only")

ggsave(filename = "Consistency_by_RT_in_HH.png", dpi = 200)

consistent_sft %>%
  filter(cell_name == "DD") %>%
  consistency_rt_plot +
  labs(subtitle = "For DD cells only")

ggsave(filename = "Consistency_by_RT_in_DD.png", dpi = 200)

short_rts <- consistent_sft %>%
  mutate(too_short = rt < 350, too_short_old = rt < 200) %>%
  group_by(sonaID) %>%
  summarise(
    num_short = sum(too_short),
    percent_short_350 = sum(too_short) / n(),
    percent_short_200 = sum(too_short_old) / n()
  )

# Visualise subject x percent_short
# 350ms chosen as responses > 350ms show differences between baseline HH accept
# rates. (That is consistency increased significantly above baseline)
cmap <- c("Below 350ms" = watercolour$zorn$discrete['chair'],
          "Below 200ms" = watercolour$zorn$discrete['boat'])
short_rts %>%
  arrange(percent_short_350) %>%
  mutate(idx = row_number(), ) %>%
  ggplot(aes(x = idx, y = percent_short_350, colour = "Below 350ms")) +
  geom_point() +
  geom_point(
    data = short_rts %>%
      arrange(percent_short_200) %>%
      mutate(idx = row_number()),
    mapping = aes(x = idx, y = percent_short_200, colour = "Below 200ms"),
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .05)) +
  geom_hline(yintercept = 0.2, colour = cmap[1], linetype = 2) +
  geom_hline(yintercept = 0.1, colour = cmap[2], linetype = 2) +
  scale_colour_manual(name = "Cut off point", values = cmap) +
  labs(
    title = "Fast responses by participants",
    y = "Percentage of trials below cutoff",
    x = "Participants ordered by percent under cutoff",
    caption = "NOTE: points along the x-axis do not necessarily represent the same participant as they are ordered separately"
    ) +
  theme(legend.position = "top")
    
ggsave(filename = "Percent_Fast_responses.png", dpi = 200)

# Drop subjects with more than 20% short RT's (< 350ms)
drop_subjects <- short_rts %>%
  filter(percent_short_350 > 0.2) %>%
  pull(sonaID)

enough_long_sft <- consistent_sft %>%
  filter(!sonaID %in% drop_subjects)

# Condition split after dropping subjects with too many short RT's
condition_split(enough_long_sft)

# Drop short RT trials
short_removed_sft <- enough_long_sft %>%
  filter(rt >= 350)

#Check trial numbers at the lower end
short_removed_sft %>%
  count(sonaID) %>%
  arrange(n) %>%
  head

# Participants ordered by percentage of RT's greater than 10s
short_removed_sft %>%
  mutate(too_long = rt > 10000) %>%
  group_by(sonaID) %>%
  summarise(percent_long = sum(too_long) / n()) %>%
  mutate(percent_long = round(percent_long, 3)) %>%
  arrange(desc(percent_long))


# Drop long RT trials
long_removed_sft <- short_removed_sft %>%
  filter(rt < 10000)

#Check trial numbers at the lower end
long_removed_sft %>%
  count(sonaID) %>%
  arrange(n) %>%
  head

# Consistency by cell
long_removed_sft %>%
  group_by(cell_name) %>%
  summarise(num_consistent = sum(consistent) / n()) %>%
  ggplot(aes(x = cell_name, y = num_consistent)) +
  geom_bar(stat = "identity")

pt_col <- watercolour$zorn$discrete['chair']
# Overall consistency (collapsed across cell types)
long_removed_sft %>%
  group_by(sonaID) %>%
  summarise(avg_consistency = mean(consistent)) %>%
  arrange(avg_consistency) %>%
  mutate(idx = row_number()) %>%
  ggplot(aes(x = idx, y = avg_consistency)) +
  geom_point(colour = pt_col) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
  geom_hline(yintercept = 0.6, colour = pt_col, linetype = 2) +
  labs(
    title = "Overall consistency by participant",
    y = "Percentage of trials consistent with expected response",
    x = "Participants ordered by consistency",
    )

ggsave(filename = "Consistency_Overall.png", dpi = 200)

# Consistency by easiest cell types
long_removed_sft %>%
  filter(cell_name %in% c("HH", "DD")) %>%
  group_by(sonaID, cell_name) %>%
  summarise(avg_consistency = mean(consistent)) %>%
  group_by(sonaID) %>%
  mutate(sort_col = mean(avg_consistency)) %>%
  group_by(sort_col, sonaID) %>%
  mutate(idx = cur_group_id()) %>%
  arrange(idx) %>%
  ggplot(aes(x = idx, y = avg_consistency, colour = cell_name)) +
  geom_point() +
  geom_line(aes(y = sort_col), colour = pt_col) +
  scale_colour_watercolour() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
  labs(
    title = "Consistency by cell by participant",
    y = "Percentage of trials consistent with expected response",
    x = "Participants ordered by HH+DD",
    caption = "The line represents HH + DD / 2",
    )

ggsave(filename = "Consistency_HH_DD.png", dpi = 200)

cell_consistency <- function(df, cell) {
  df %>%
    filter(cell_name == cell) %>%
    group_by(sonaID) %>%
    summarise(avg_consistency = mean(consistent)) %>%
    arrange(avg_consistency) %>%
    mutate(idx = row_number())
}

long_removed_sft %>%
  cell_consistency("HH") %>%
  bind_rows(cell_consistency(long_removed_sft, "DD"), .id = "cell") %>%
  mutate(cell = ifelse(cell == 1, "HH", "DD")) %>%
  ggplot(aes(x = idx, y = avg_consistency, colour = cell)) +
  geom_point() +
  scale_colour_watercolour() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
  labs(
    title = "Consistency by cell by participant",
    y = "Percentage of trials consistent with expected response",
    x = "Participants ordered by HH or DD separately",
    )

ggsave(filename = "Consistency_HH_DD_separate.png", dpi = 200)

drop_subjects <- long_removed_sft %>%
  cell_consistency("HH") %>%
  filter(avg_consistency < 0.75) %>%
  pull(sonaID)

drop_subjects <- long_removed_sft %>%
  cell_consistency("DD") %>%
  filter(avg_consistency < 0.75) %>%
  pull(sonaID) %>%
  c(drop_subjects) %>%
  unique

final_sft <- long_removed_sft %>%
  filter(!sonaID %in% drop_subjects)

#Final condition split after removing inconsistent participants
condition_split(final_sft)
condition_split(long_removed_sft)

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

# Note: Of the final participants, 8 ran out of time before the task finished

# For Guy's categorisation modelling?

final_sft %>%
  mutate(rt = as.numeric(rt) / 1000,
         subject = as.factor(sonaID),
         accept = as.logical(accept),
         price = round(as.numeric(price), 2),
         rating = round(as.numeric(rating), 2)) %>%
  separate(cell_name, into = c("price_lvl", "rating_lvl"), sep = 1) %>%
  relocate(subject, rt, accept, price, rating) %>%
  mutate(price_lvl = factor(price_lvl, levels = c("H", "L", "D")),
         rating_lvl = factor(rating_lvl, levels = c("H", "L", "D"))) %>%
  select(-c(sonaID, error, price_threshold, rating_threshold, key_press)) %>%
  saveRDS(file = here::here("data", "output", "Gavins_PrefData_for_Guy.RDS"))

saveRDS(pref_data, file = here::here("data", "output", "Pref_preprocessed.RDS"))
