library(mcce)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(watercolours)
library(readr)
library(forcats)
library(tidyr)
library(pmwg)
library(purrr)
library(ggcorrplot)


# Load all the samples
pref_file <- here::here("data", "output", "PrefDCE_3Yaob1kCATGm_staged_burn.RData")
pref_file <- here::here("data", "output", "PrefDCE_S6I7q4nycmRv_short_burn_cont.RData")
pref_file <- here::here("data", "output", "PrefDCE_S6I7q4nycmRv_short_burn_cont_2.RData")
pref_file <- here::here("data", "output", "PrefDCE_ACUGlaCbumcx_reduced_continue.RData")

viewstage <- "sample"

pref_sampler <- get_samples(pref_file, final_obj="sampler")
parameters <- pref_sampler$par_names
archs <- parameters[startsWith(parameters, "alpha")]
subjects <- pref_sampler$subjects
pref_samples <- extract_parameters(pref_sampler)


mean_window <- 150
pdf(file = here::here("results", "Reduced", paste0("theta_mu_trace_", Sys.Date(), ".pdf")), width = 14.1, height = 7.53)
for (par in parameters) {
  g <- pref_samples %>%
    filter(subjectid == "theta_mu", parameter == par) %>%
    mutate(mean_so_far = ifelse(sampleid < mean_window, NaN, RcppRoll::roll_meanr(value, n=mean_window))) %>%
    trace_plot(par) +
    geom_line(aes(y = mean_so_far), colour = "black") +
    scale_colour_watercolour()
  print(g)

  if (interactive()) {
    readline(prompt="Press [enter] to continue")
  }
}
dev.off()

for (par in parameters) {
  print(paste("Generating pdf for", par))
  pdf(file = here::here("results", "Reduced", paste0(par, "_randeff_trace_", Sys.Date(), ".pdf")), width = 14.1, height = 7.53)
  for (subj_arr in split(subjects, ceiling(seq_along(subjects) / 4))) {
    g <- pref_samples %>%
      filter(subjectid %in% subj_arr) %>%
      trace_plot(par) +
      scale_colour_watercolour() +
      facet_wrap(~ subjectid, nrow=2, ncol=2)
    print(g)

    if (interactive()) {
      readline(prompt="Press [enter] to continue")
    }
  }
  dev.off()
}

model_medians <- pref_samples %>%
  filter(stageid == viewstage, parameter %in% archs) %>%
  get_summary(tform = exp) %>%
  arch_medians()

arch_order <- model_medians %>%
  group_by(parameter) %>%
  summarise(mean_val = mean(rel_val)) %>%
  arrange(mean_val) %>%
  pull(parameter)

most_common_arch <- arch_order[length(arch_order)]

subject_order <- model_medians %>%
  filter(parameter == most_common_arch) %>%
  arrange(desc(rel_val)) %>%
  pull(subjectid)

model_medians %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(subjectid = fct_relevel(subjectid, "Group", after=Inf)) %>%
  filter(subjectid != "Group") %>%
  mutate(parameter = factor(parameter, arch_order)) %>%
  arch_plot() +
  scale_fill_watercolour()

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
  filter(stageid == viewstage, !(parameter %in% archs)) %>%
  get_summary()

par_colours <- c("t0" = "grey",
                 "A" = "grey", "b_acc" = "#73842E", "b_rej" = "#D0781C",
                 "v_acc_p_H" = "#3b8114", "v_acc_p_L" = "#67b13d", "v_acc_p_D" = "#95db6d",
                 "beta0_p" = "#008252", "beta1_p" = "#008252",
                 "v_acc_r_H" = "#402c8e", "v_acc_r_L" = "#6853bb", "v_acc_r_D" = "#957df2",
                 "beta0_r" = "#57057e", "beta1_r" = "#57057e")

par_medians %>%
  filter(subjectid != "theta_mu") %>%
  mutate(colour = par_colours[parameter]) %>%
  mutate(value = ifelse(startsWith(parameter, "v"), value, exp(value))) %>%
  mutate(parameter = factor(parameter, levels = names(par_colours))) %>%
  ggplot(aes(x = parameter, y = value, fill = colour)) +
  geom_boxplot() +
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

participant <- sample(subjects, 1)
participant_data <- pref_sampler$data %>% filter(subject == participant)
participant_random_effects <- pref_samples %>%
  filter(stageid == viewstage, subjectid == participant) %>%
  print
participant_arch <- model_medians %>% filter(subjectid == participant)

winning_arch <- participant_arch %>%
  arrange(desc(value)) %>%
  pull(parameter) %>%
  head(1)

nonarch_parameters <- par_medians %>%
  filter(subjectid == participant) %>%
  ungroup() %>%
  select(-subjectid) %>%
  pivot_wider(names_from=parameter, values_from=value) %>%
  unlist()

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
      y <- transform_pars(x, tforms = "reduced")
      trial_ll <- model_wrapper(y, participant_data, ll_func)
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


group_level <- par_medians %>%
  filter(subjectid == "theta_mu") %>%
  ungroup() %>%
  select(-subjectid) %>%
  pivot_wider(names_from=parameter, values_from=value) %>%
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
  geom_col() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

cor_df <- pref_sampler %>%
  extract_cov(filter = viewstage) %>%
  apply(c(1,2), mean) %>%
  cov2cor %>%
  matrix_to_df

cor_df %>%
  matdf_plot() +
  scale_fill_watercolour(type = "diverging")

pref_sampler %>%
  `[[`("samples") %>%
  `[[`("subj_ll") %>%
  t %>%
  data.frame %>%
  tibble %>%
  mutate(idx = row_number()) %>%
  pivot_longer(cols=-idx, names_to="subject", values_to='ll') %>%
  ggplot(aes(x = idx, y = ll, group=subject)) +
  geom_line(colour = "black", alpha=0.15)

sigmas <- pref_sampler %>%
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
