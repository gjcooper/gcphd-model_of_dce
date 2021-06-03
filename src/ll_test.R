library(dplyr)
library(tidyr)
pkg <- asNamespace("mcce")


medians <- readRDS(here::here("data", "output", "median_alpha_exp1.RDS"))
example_vals <- medians[5, ] %>% round(1)
test_vals <- example_vals %>%
  dplyr::select(!(starts_with("v") | starts_with("alpha")))

cell <- "HH"
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
  pivot_wider(names_from=accattr, values_from=value)
  
drifts

# Generate samples using rll_func
n_rsamples <- 1e6 + 8
n_outer <- n_rsamples / 3
n_inner <- n_outer / 3
test_trials <- tibble(
  price = as.factor(rep(c("H", "L", "D"), each = n_outer)),
  rating = as.factor(rep(c("H", "L", "D"), each = n_inner, times = 3))
  ) %>%
  mutate(
    v_acc_p = paste0("v_acc_p_", price),
    v_rej_p = paste0("v_rej_p_", price),
    v_acc_r = paste0("v_acc_r_", rating),
    v_rej_r = paste0("v_rej_r_", rating)
  )

test_data <- rmodel_wrapper(test_vals, test_trials, pkg$rll_IST)

integrate_args <- c(
  as.list(test_vals),
  drifts = list(drifts),
  f = simple_model_ll,
  lower = unname(unlist(exp(test_vals["t0"]))) + 1e-3,
  upper = Inf,
  accept = TRUE,
  model = pkg$ll_IST
  )

inttest <- c(
  as.list(test_vals),
  drifts = list(drifts)
  )

do.call(integrate, integrate_args)

integrate(f, lower, upper, ..., subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)
