library(mcce)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(paletti)
library(readr)
library(forcats)
library(tidyr)
library(pmwg)
library(purrr)
library(ggcorrplot)

fw_cols <- read_lines(file = "palette.txt")

frank_colmap <- scale_fill_manual(name = "Architecture", values = fw_cols, breaks = names_ll())

# Load all the samples
pref_file <- here::here("data", "output", "PrefDCE_3Yaob1kCATGm_staged_burn.RData")
pref_file <- here::here("data", "output", "PrefDCE_ACUGlaCbumcx_reduced_continue.RData")

viewstage <- "burn"

pref_samples <- get_samples(pref_file, final_obj="sampler")


pdf(file = here::here("results", "Reduced", paste0("theta_mu_trace_", Sys.Date(), ".pdf")), width = 14.1, height = 7.53)
for (par in pref_samples$par_names) {
  g <- pref_samples %>%
    as_mcmc %>%
    data.frame %>%
    tibble %>%
    mutate(sample_id = row_number()) %>%
    mutate(stage = pref_samples$samples$stage) %>%
    pivot_longer(cols = -c(sample_id, stage), names_to = "parameter") %>%
    filter(parameter == par) %>%
    ggplot(aes(x = sample_id, y = value, colour = stage)) +
    geom_line() +
    scale_colour_manual(values = fw_cols) +
    labs(title = par) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank())
    print(g)
}
dev.off()

random_effects <- pref_samples %>%
  as_mcmc(selection = "alpha") %>%
  lapply(FUN = function(x) {
      x %>% data.frame() %>% tibble()
  }) %>%
  bind_rows(.id = "subjectid") %>%
  tibble()

for (par in pref_samples$par_names) {
  print(paste("Generating pdf for", par))
  pdf(file = here::here("results", "Reduced", paste0(par, "_randeff_trace_", Sys.Date(), ".pdf")), width = 14.1, height = 7.53)
  all_subjs <- unique(random_effects$subjectid)
  for (subj_arr in split(all_subjs, ceiling(seq_along(all_subjs) / 4))) {
    g <- random_effects %>%
      filter(subjectid %in% subj_arr) %>%
      group_by(subjectid) %>%
      mutate(sample_id = row_number()) %>%
      mutate(stage = pref_samples$samples$stage) %>%
      ungroup() %>%
      pivot_longer(cols = -c(subjectid, sample_id, stage), names_to = "parameter") %>%
      filter(parameter == par) %>%
      ggplot(aes(x = sample_id, y = value, colour = stage)) +
      geom_line() +
      scale_colour_manual(values = fw_cols) +
      facet_wrap(~ subjectid, nrow=2, ncol=2) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank())
    print(g)
  }
  dev.off()
}

model_medians <- pref_samples %>%
  extract_parameters(str_subset(.$par_names, "alpha"), filter = viewstage) %>%
  get_medians() %>%
  group_by(subjectid) %>%
  mutate(rel_val = value / sum(value)) %>%
  mutate(subjectid = case_when(
    subjectid == "theta_mu" ~ "Group",
    TRUE ~ str_pad(subjectid, 2, pad = "0")
  ))

subject_order <- model_medians %>%
  filter(Parameter == "IEX") %>%
  arrange(desc(rel_val)) %>%
  pull(subjectid)

Par_order <- model_medians %>%
  group_by(Parameter) %>%
  summarise(mean_val = mean(rel_val)) %>%
  arrange(mean_val) %>%
  pull(Parameter)

model_plot <- function(medians) {
  Par_order <- c("FPP", "IST", "IEX", "CB", "MW")
  medians %>%
    filter(subjectid != "Group") %>%
    mutate(Parameter = factor(Parameter, Par_order)) %>%
    ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    frank_colmap +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(y = "Relative Evidence") +
    scale_y_continuous(labels = NULL, breaks = NULL)
}

model_medians %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(subjectid = fct_relevel(subjectid, "Group", after=Inf)) %>%
  filter(subjectid != "Group") %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  model_plot

ggsave(
  filename = here::here(
    "results",
    "Reduced",
    paste0("EstimatedArch_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)

par_medians <- pref_samples %>%
  extract_parameters(str_subset(.$par_names, "alpha", negate = TRUE), filter=viewstage) %>%
  filter(subjectid != 'theta_mu') %>%
  get_medians(alpha = FALSE) %>%
  mutate(value = log(value))

group_pars <- pref_samples %>%
  extract_parameters(str_subset(.$par_names, "alpha", negate = TRUE), filter=viewstage) %>%
  filter(subjectid == 'theta_mu') %>%
  get_medians(alpha = FALSE) %>%
  mutate(value = log(value))


fw_cols %>%
  as_tibble %>%
  ggplot(aes(x = seq_len(length(fw_cols)), y = 1, fill=value)) +
  geom_col() + scale_fill_identity() +
  labs(x = "Index", y = NULL) + theme(axis.text.y = element_blank())

par_colours <- c("t0" = "grey",
                 "A" = "grey", "b_acc" = fw_cols[3], "b_rej" = fw_cols[2],
                 "v_acc_p_H" = "#3b8114", "v_acc_p_L" = "#67b13d", "v_acc_p_D" = "#95db6d",
                 "beta0_p" = "#008252", "beta1_p" = "#008252",
                 "v_acc_r_H" = "#402c8e", "v_acc_r_L" = "#6853bb", "v_acc_r_D" = "#957df2",
                 "beta0_r" = "#57057e", "beta1_r" = "#57057e")

par_medians %>%
  mutate(colour = par_colours[Parameter]) %>%
  mutate(value = ifelse(startsWith(Parameter, "v"), value, exp(value))) %>%
  mutate(Parameter = factor(Parameter, levels = names(par_colours))) %>%
  ggplot(aes(x = Parameter, y = value, fill = colour)) +
  geom_boxplot() +
  ylim(c(-5, 15)) +
  scale_fill_identity()

ggsave(
  filename = here::here(
    "results",
    "Reduced",
    paste0("EstimatedRE_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)

set.seed(2468)

participant <- sample(pref_samples$subjects, 1)
participant_data <- pref_samples$data %>% filter(subject == participant)
participant_random_effects <- pref_samples %>%
  as_mcmc(selection = "alpha", filter = viewstage) %>%
  pluck(participant) %>%
  data.frame() %>%
  tibble() %>%
  print

means <- participant_random_effects %>%
  summarise(across(everything(), mean))

winning_arch <- means %>%
  select(starts_with('alpha')) %>%
  rowwise() %>%
  mutate(max_alpha = names(.)[which.max(across(everything()))]) %>%
  pull(max_alpha) %>%
  str_split("_") %>%
  unlist %>%
  `[`(2)

nonarch_parameters <- means %>%
  select(!starts_with("alpha")) %>%
  unlist


ll_func <- select_ll(winning_arch)
stim_levels <- c("H", "L", "D")

ll_vals <- sapply(names(nonarch_parameters), FUN = function(par) {
  par_val <- nonarch_parameters[par]
  alt_vals <- seq(par_val - 1, par_val + 1, 0.05)
  x <- nonarch_parameters
  tibble(
    par_val = alt_vals,
    ll = sapply(alt_vals, FUN = function(sub) {
      x[par] <- sub
      force_pos_mask <- !startsWith(names(x), "v")
      y = x
      y[force_pos_mask] <- exp(x[force_pos_mask])
      trial_ll <- reduced_model_wrapper(y, participant_data, ll_func)
      sum(log(pmax(trial_ll, 1e-10)))
    })
  )
}, simplify=FALSE, USE.NAMES=TRUE) %>%
bind_rows(.id="parameter")

par_tbl <- tibble(parameter=names(nonarch_parameters), value=nonarch_parameters)

ll_vals %>%
  ggplot(aes(x = par_val, y = ll)) +
  geom_line() +
  geom_point() +
  geom_vline(data=par_tbl, aes(xintercept=value)) +
  facet_wrap(~ parameter, scales="free")

ggsave(
  filename = here::here(
    "results",
    "Reduced",
    paste0("ll_test_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 10,
  height = 10,
  units = "in"
)


group_level <- pref_samples %>%
  as_mcmc(selection = "theta_mu", filter = viewstage) %>%
  data.frame() %>%
  tibble() %>%
  summarise(across(everything(), mean)) %>%
  select(-starts_with("alpha")) %>%
  mutate(
    v_rej_p_H_derived = exp(beta0_p) - exp(beta1_p)*v_acc_p_H,
    v_rej_p_L_derived = exp(beta0_p) - exp(beta1_p)*v_acc_p_L,
    v_rej_p_D_derived = exp(beta0_p) - exp(beta1_p)*v_acc_p_D,
    v_rej_r_H_derived = exp(beta0_r) - exp(beta1_r)*v_acc_r_H,
    v_rej_r_L_derived = exp(beta0_r) - exp(beta1_r)*v_acc_r_L,
    v_rej_r_D_derived = exp(beta0_r) - exp(beta1_r)*v_acc_r_D,
  ) %>%
  pivot_longer(cols=everything()) %>%
  mutate(value = ifelse(startsWith(name ,"v"), value, exp(value))) %>%
  print

group_level %>%
  ggplot(aes(x = name, y = value)) +
  geom_col()

covmat <- pref_samples %>%
  as_mcmc(selection = "theta_sig", filter = viewstage) %>%
  data.frame() %>%
  tibble %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(everything(), names_to = c("X", "Y"), names_sep = "\\.") %>%
  mutate(across(c(X, Y), factor)) %>%
  xtabs(value ~ X + Y, .)

ggcorrplot(covmat)
ggcorrplot(cov2cor(covmat))
ggcorrplot(cov2cor(covmat), hc.order = TRUE, type = "lower", outline.col = "white", lab=TRUE)
ggcorrplot(cov2cor(covmat), type = "lower", outline.col = "white", lab=TRUE)


pref_samples %>%
  `[[`("samples") %>%
  `[[`("subj_ll") %>%
  t %>%
  data.frame %>%
  tibble %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols=-idx, names_to="subject", values_to='ll') %>%
  ggplot(aes(x = idx, y = ll, group=subject)) +
  geom_line(colour = "black", alpha=0.15)

sigmas <- pref_samples %>%
  as_mcmc(selection = "theta_sig") %>%
  data.frame %>%
  tibble %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols=-idx, names_sep="\\.", names_to=c("X", "Y"), values_to="sigma")

sigmas %>%
  mutate(XY = paste0(X, Y)) %>%
  mutate(col = ifelse(X == Y, "blue", "black")) %>%
  ggplot(aes(x = idx, y = log(sigma), group=XY, colour = col)) +
  geom_line(alpha=0.25) +
  geom_line(data = sigmas %>% filter(X == Y), aes(group=X), colour="blue", alpha=0.25) +
  scale_colour_identity(labels=c("covariance", "variance"), guide = "legend")
