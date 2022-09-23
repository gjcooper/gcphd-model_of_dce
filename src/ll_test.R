library(dplyr)
library(tidyr)
pkg <- asNamespace("mcce")
library(mcce)
devtools::load_all()
library(rtdists)
library(ggplot2)
library(patchwork)
library(parallel)
library(cubature)

cells <- c("HH", "HL", "LH", "LL", "HD", "LD", "DH", "DL", "DD")

get_plausible_vals <- function() {
  medians <- readRDS(here::here("data", "output", "median_alpha_exp1.RDS"))
  medians %>%
    select(-subjectid) %>%
    summarise_all(mean) %>%
    round(1) %>%
    mutate(t0 = -Inf)
}

get_modded_vals <- function() {
  original <- get_plausible_vals()
  original %>%
    exp %>%
    mutate(across(c(v_rej_p_H, v_rej_r_H, v_rej_p_L, v_rej_r_L), ~ .x + 1)) %>%
    mutate(b_rej = b_rej / 4) %>%
    mutate(v_acc_p_H = 1) %>%
    log
}

get_rearranged_vals <- function(example_vals) {
  test_vals <- example_vals %>%
    dplyr::select(!(starts_with("v") | starts_with("alpha")))

  drifts <- sapply(cells, function(cell) {
    drifts <- example_vals %>%
      dplyr::select(starts_with("v")) %>%
      pivot_longer(
        cols = everything(),
        names_sep = "_",
        names_to = c("drift", "accept", "attr", "level")
      ) %>%
      mutate(
        attr = ifelse(attr == "p", "Price", "Rating"),
        accept = ifelse(accept == "acc", "Acc", "Rej")
      ) %>%
      mutate(accattr = paste0(accept, attr)) %>%
      dplyr::select(-drift) %>%
      filter(if_else(
        attr == "Price" & level == substring(cell, 1, 1) |
        attr == "Rating" & level == substring(cell, 2, 2),
        TRUE,
        FALSE
      )) %>%
      dplyr::select(accattr, value) %>%
      pivot_wider(names_from = accattr, values_from = value)
  }, simplify = FALSE)
  list(
    all = example_vals,
    fixed = test_vals,
    drifts = drifts
  )
}

# Generate samples using rll_func
get_predicted_data <- function(trials_per_cell, pars, func) {
  test_trials <- tibble(
    price = as.factor(rep(c("H", "L", "D"), each = trials_per_cell * 3)),
    rating = as.factor(rep(c("H", "L", "D"), each = trials_per_cell, times = 3))
    ) %>%
    mutate(
      v_acc_p = paste0("v_acc_p_", price),
      v_rej_p = paste0("v_rej_p_", price),
      v_acc_r = paste0("v_acc_r_", rating),
      v_rej_r = paste0("v_rej_r_", rating)
    )

  test_data <- rmodel_wrapper(pars, test_trials, func) %>%
    mutate(cell = paste0(price, rating)) %>%
    mutate(cell = factor(cell, levels = cells))
}

run_integrate <- function(drifts, pars, func) {
  accept_int <- sapply(drifts, function(drifts_tbl) {
    integrate_args <- c(
      as.list(pars),
      drifts = list(drifts_tbl),
      f = simple_model_wrapper,
      lower = 0,
      upper = 100,
      accept = TRUE,
      model = func,
      relTol = 1e-15
      )
    int_res <- do.call(cubintegrate, integrate_args)
    int_res$integral
  })
  reject_int <- sapply(drifts, function(drifts_tbl) {
    integrate_args <- c(
      as.list(pars),
      drifts = list(drifts_tbl),
      f = simple_model_wrapper,
      lower = 0,
      upper = 100,
      accept = FALSE,
      model = func,
      relTol = 1e-15
      )
    int_res <- do.call(cubintegrate, integrate_args)
    int_res$integral
  })
  bind_cols(accept = accept_int, reject = reject_int, cell = names(accept_int))
}

get_proportions <- function(pp_data, int_data) {
  sample_prop <- pp_data %>%
    count(cell, accept) %>%
    group_by(cell) %>%
    mutate(prop = prop.table(n)) %>%
    filter(accept == 2) %>%
    dplyr::select(-c(accept, n))

  combined <- int_data %>%
    full_join(sample_prop, by = "cell") %>%
    relocate(c(cell, accept, reject, prop)) %>%
    mutate(adiff = prop - accept, rdiff = prop - (1 - reject))

  combined
}

new_IST <- function(rt, A, b_acc, b_rej, t0, drifts, accept) { # nolint
  if (accept) {
    ll <- dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)  *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)  *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1))
  }
  ll
}

new_IEX <- function(rt, A, b_acc, b_rej, t0, drifts, accept) { # nolint
  if (accept) {
    ll <- dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1))
  }
  ll
}

full_MW <- function(rt, A, b_acc, b_rej, t0, drifts, accept) { # nolint
  if (accept) {
    ll <- dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)  *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) +
          dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)  *
          (1 - plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1))
  } else {
    ll <- dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)  *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)  *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)) *
          plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1) +
          dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1) *
          plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)  *
          plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1) *
          (1 - plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1))
  }
  ll
}

function_list <- list(
  old_IST = list(ll = pkg$ll_IST, rsample = pkg$rll_IST),
  new_IST = list(ll = new_IST, rsample = pkg$rll_IST),
  old_IEX = list(ll = pkg$ll_IEX, rsample = pkg$rll_IEX),
  new_IEX = list(ll = new_IEX, rsample = pkg$rll_IEX),
  old_CB = list(ll = pkg$ll_CB, rsample = pkg$rll_CB),
  new_FPP = list(ll = pkg$ll_FPP, rsample = pkg$rll_FPP),
  new_MW = list(ll = pkg$ll_MW, rsample = pkg$rll_MW),
  full_MW = list(ll = full_MW, rsample = pkg$rll_MW)
)


# Start the clock!
# 2 cores, 10k samples, 31 minutes
ptm <- proc.time()
all_runs <- lapply(function_list, function(func_pair) {
  pars <- get_modded_vals() %>% get_rearranged_vals()
  pred <- get_predicted_data(1e2, pars$all, func_pair$rsample)
  intf <- run_integrate(pars$drifts, pars$fixed, func_pair$ll)
  res <- get_proportions(pred, intf)
  res
})
# Stop the clock
proc.time() - ptm

combined <- bind_rows(all_runs, .id = "function_set") %>%
  mutate(
    cell = factor(cell, levels = cells),
    function_set = factor(function_set, levels = names(function_list))
  )
acc_diff <- combined %>%
  ggplot(aes(
    x = cell,
    y = adiff,
    colour = function_set,
    group = function_set
  )) +
  geom_line() +
  labs(
    title = "% Accept in samples - integrate using accept = TRUE",
    x = "Cell from Design",
    y = "Difference"
  ) +
  ylim(-1, 1)
rej_diff <- combined %>%
  ggplot(aes(
    x = cell,
    y = rdiff,
    colour = function_set,
    group = function_set
  )) +
  geom_line() +
  labs(
    title = "% Accept in samples - (1 - integrate using accept = FALSE)",
    x = "Cell from Design",
    y = "Difference"
  ) +
  ylim(-1, 1)
acc_diff / rej_diff

acc_abs <- combined %>%
  ggplot(aes(
    x = cell,
    y = accept,
    colour = function_set,
    group = function_set
  )) +
  geom_line() +
  labs(
    title = "Integrate Accept",
    x = "Cell from Design",
    y = "Accept percentage"
  ) +
  ylim(0, 1)
rej_abs <- combined %>%
  ggplot(aes(
    x = cell,
    y = reject,
    colour = function_set,
    group = function_set
  )) +
  geom_line() +
  labs(
    title = "Integrate Reject",
    x = "Cell from Design",
    y = "Reject percentage"
  ) +
  ylim(0, 1)
acc_abs / rej_abs

saveRDS(combined, "combined.RDS")
