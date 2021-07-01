library(dplyr)
library(mcce)
library(tidyr)

fac_cols <- c("cell_name", "sonaID", "group", "hand")
sft_columns <- c("rt", "group", "hand", "sonaID", "accept", "cell_name")

datafile <- here::here("data", "input", "PrefSFT", "results_20201116155132_sem2.txt")
raw <- get_JATOS_data(datafile)
extracted <- extract_data(raw)

pref_data <- extracted %>%
  filter(type == "sft") %>%
  mutate_at(fac_cols, as.factor) %>%
  select(all_of(sft_columns)) %>%
  mutate(rt = as.numeric(rt) / 1000) %>%
  mutate(subject = as.factor(sonaID)) %>%
  mutate(accept = as.numeric(as.logical(accept)) + 1) %>%
  separate(cell_name, into = c("price", "rating"), sep = 1) %>%
  relocate(subject, rt, accept, price, rating) %>%
  mutate(price = factor(price, levels = c("H", "L", "D")),
         rating = factor(rating, levels = c("H", "L", "D"))) %>%
  select(-sonaID) %>%
  mutate(
    v_acc_p = paste0("v_acc_p_", price),
    v_rej_p = paste0("v_rej_p_", price),
    v_acc_r = paste0("v_acc_r_", rating),
    v_rej_r = paste0("v_rej_r_", rating)
  )

saveRDS(pref_data, file = here::here("data", "output", "Pref_preprocessed.RDS"))
