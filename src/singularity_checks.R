library(tidyverse)
library(pmwg)
library(mcce)

Sys.setenv(MCCE_EST_EXP="PrefDCE",
           NCPUS=1,
           MCCE_EXP_DATA="Pref_preprocessed.RDS",
           MCCE_MIN_RT=0.35,
           MCCE_MAX_RT=10,
           MCCE_CONTAM=0.02,
           MCCE_TAG="singularity_test",
           RANDOM_SEED=101,
           MCCE_MODEL="reduced",
           MCCE_METHOD="model")

vars <- check_env()

model_run <- "stdneg"
pref_file <- here::here("data", "output", "EstimationMarchProblems", "StdNeg", "pmwg_obj_2c262061c5.RDS")

model_run <- "std"
pref_file <- here::here("data", "output", "EstimationMarchProblems", "Normal", "pmwg_obj_6af320275c9a.RDS")

model_run <- "reduced"
pref_file <- here::here("data", "output", "EstimationMarchProblems", "Reduced", "pmwg_obj_49d242101de7.RDS")

pref_sampler <- readRDS(pref_file)

parameters <- pref_sampler$par_names
archs <- parameters[startsWith(parameters, "alpha")]
subjects <- pref_sampler$subjects
pref_samples <- extract_parameters(pref_sampler)

par_list <- c("A", "alpha_CB", "b_acc", "v_acc_p_H", "v_rej_r_L", "v_acc_r_D")
for (par in par_list) {
  g <- pref_samples %>%
    filter(subjectid == "theta_mu") %>%
    trace_plot(par)
  print(g)
  if (interactive()) {
    readline(prompt = "Press [enter] to continue")
  }
}

subj_arr <- sample(split(subjects, ceiling(seq_along(subjects) / 4)), 1) %>%
  pluck(1)

dev.new()

for (par in par_list) {
  g <- pref_samples %>%
    filter(subjectid %in% subj_arr) %>%
    trace_plot(par) +
    facet_wrap(~ subjectid, nrow = 2, ncol = 2)
  print(g)
  if (interactive()) {
    readline(prompt = "Press [enter] to continue")
  }
}


par_medians <- pref_samples %>%
  filter(!(parameter %in% archs)) %>%
  get_summary() %>%
  mutate(tvalue = ifelse(startsWith(parameter, "v"), value, exp(value)))

par_colours <- c("t0" = "grey",
                 "A" = "grey", "b_acc" = "#73842E", "b_rej" = "#D0781C",
                 "v_acc_p_H" = "#569F72", "v_acc_p_L" = "#407755", "v_acc_p_D" = "#2B5039",
                 "v_rej_p_H" = "#F29F40", "v_rej_p_L" = "#EF8C1A", "v_rej_p_D" = "#BF6C0D",
                 "v_acc_r_H" = "#85C0FF", "v_acc_r_L" = "#47A0FF", "v_acc_r_D" = "#0A81FF",
                 "v_rej_r_H" = "#BBA9A0", "v_rej_r_L" = "#A2887C", "v_rej_r_D" = "#83685D")

dev.new()
par_medians %>%
  filter(subjectid != "theta_mu") %>%
  mutate(colour = par_colours[parameter]) %>%
  ggplot(aes(x = parameter, y = value, fill = colour)) +
  geom_boxplot() +
  scale_fill_identity()

dev.new()
par_medians %>%
  filter(subjectid == "theta_mu") %>%
  mutate(colour = par_colours[parameter]) %>%
  ggplot(aes(x = parameter, y = value, fill = colour)) +
  geom_col() +
  scale_fill_identity()


ll <- pref_sampler %>%
  pluck("samples", "subj_ll") %>%
  t() %>%
  as_tibble() %>%
  mutate(idx = row_number()) %>%
  pivot_longer(-idx, names_to="subject", values_to="ll")

dev.new()
ll %>%
  ggplot(aes(x = idx)) +
  geom_line(aes(y = ll, group=subject))

eps_ll <- pref_sampler %>%
  pluck("samples", "epsilon") %>%
  t() %>%
  as_tibble() %>%
  mutate(idx = row_number()) %>%
  pivot_longer(-idx, names_to="subject", values_to="epsilon") %>%
  left_join(ll, by=c("idx", "subject"))
  
dev.new()
eps_ll %>%
  ggplot(aes(x = idx)) +
  geom_line(aes(y = epsilon, group=subject))

covars <- pref_sampler %>% pluck("samples", "theta_sig")

covar_df <- apply(covars, 3, function(x) {
    as_tibble(data.frame(x), rownames="variable") %>%
      pivot_longer(-variable)}) %>%
  bind_rows(.id = "sample") %>%
  mutate(sample = as.integer(sample))

for (par in par_list) {
  dev.new()
  p <- covar_df %>%
    filter(variable == par) %>%
    ggplot(aes(x = sample, y = value)) +
    geom_line() +
    labs(x = "iteration", y = "variance") +
    facet_wrap(~ name, scale = "free") +
    labs(title = par)
  print(p)

  if (interactive()) {
    readline(prompt = "Press [enter] to continue")
  }
}


input_data <- pref_sampler$data

sub_splits <- split(subjects, ceiling(seq_along(subjects) / 4))

pc_acc <- input_data %>%
  count(subject, accept) %>%
  group_by(subject) %>%
  mutate(pc = n/sum(n)) %>%
  select(subject, accept, pc)

for (sub_arr in sub_splits) {
  print(sub_arr)
  p <- input_data %>%
    left_join(pc_acc, by = c("subject", "accept")) %>%
    filter(subject %in% sub_arr) %>%
    select(-starts_with("v")) %>%
    mutate(accept = factor(accept, levels = c(2, 1), labels = c("Accept", "Reject"))) %>%
    group_by(subject, accept) %>%
    arrange(rt) %>%
    mutate(cdfsf = seq(from = 0, to = 1, length.out = n())) %>%
    mutate(adjcdf = cdfsf * pc) %>%
    ggplot(aes(x = rt, y = adjcdf, colour = accept)) +
    geom_hline(yintercept = 0, colour = "grey", linetype = 2) +
    geom_hline(yintercept = 4/9, colour = "#10bdcb", linetype = 2) +
    geom_hline(yintercept = 6/9, colour = "#1063cb", linetype = 2) +
    geom_hline(yintercept = 5/9, colour = "#88069f", linetype = 2) +
    geom_line() +
    scale_color_manual(values = c("Accept" = "#138138", "Reject" = "#b10606")) +
    facet_wrap(~ subject)
  print(p)

  if (interactive()) {
    readline(prompt = "Press [enter] to continue")
  }
}

for (sub in unique(input_data$subject)) {
  print(sub)
  pc_acc <- input_data %>%
    filter(subject == sub) %>%
    count(price, rating, accept) %>%
    group_by(price, rating) %>%
    mutate(pc = n/sum(n)) %>%
    select(accept, price, rating, pc)

  p <- input_data %>%
    filter(subject == sub) %>%
    left_join(pc_acc, by = c("price", "rating", "accept")) %>%
    select(-starts_with("v")) %>%
    mutate(accept = factor(accept, levels = c(2, 1), labels = c("Accept", "Reject"))) %>%
    group_by(price, rating, accept) %>%
    arrange(rt) %>%
    mutate(cdfsf = seq(from = 0, to = 1, length.out = n())) %>%
    mutate(adjcdf = cdfsf * pc) %>%
    ggplot(aes(x = rt, y = adjcdf, colour = accept)) +
    geom_hline(yintercept = 0, colour = "grey", linetype = 2) +
    geom_hline(yintercept = 1, colour = "grey", linetype = 2) +
    geom_line() +
    scale_color_manual(values = c("Accept" = "#138138", "Reject" = "#b10606")) +
    facet_grid(rows = vars(price), cols = vars(rating))
  print(p)

  if (interactive()) {
    readline(prompt = "Press [enter] to continue")
  }
}


  pull(rt) %>%
  ecdf %>%
  plot


## Run stage function expanded
# Func Args
pmwgs <- pref_sampler
stage <- "burn"
iter <- 10
particles <- 20
n_cores <- 1
n_unique <- NA
epsilon <- NULL
p_accept <- .8
mix <- NULL
pdist_update_n <- NA

subj_epsilon <- pmwgs$samples$epsilon[, pmwgs$samples$idx]
if (is.null(subj_epsilon)) {
  message("ERROR: no subject specific epsilon values found in sampler object")
  stop("Try running augment_sampler_epsilon(sampler) first")
}

if (is.na(subj_epsilon[1])) {
  epsilon <- set_epsilon(pmwgs$n_pars, epsilon)
  subj_epsilon <- rep(epsilon, pmwgs$n_subjects)
}

mix <- pmwg:::set_mix(stage, mix)

# Hyper parameters for epsilon tuning.
# See Garthwaite, P. H., Fan, Y., & Sisson, S. A. (2016).
alpha_star <- -stats::qnorm(p_accept / 2)
n0 <- round(5 / (p_accept * (1 - p_accept)))
# Set necessary local variables
.unique_inc <- 20
apply_fn <- lapply
# Set stable (fixed) new_sample argument for this run
stable_args <- list(
  X = 1:pmwgs$n_subjects,
  FUN = pmwg:::new_sample,
  data = pmwgs$data,
  num_particles = particles,
  # parameters argument  will be generated each iteration below
  # efficient arguments will be generated if needed below
  mix_proportion = mix,
  likelihood_func = pmwgs$ll_func,
  subjects = pmwgs$subjects
)
if (n_cores > 1) {
  apply_fn <- parallel::mclapply
  stable_args$mc.cores <- n_cores # nolint
}

# Display stage to screen
msgs <- list(
  burn = "Phase 1: Burn in\n",
  adapt = "Phase 2: Adaptation\n",
  sample = "Phase 3: Sampling\n"
)
cat(msgs[[stage]])

# Build new sample storage
pmwgs <- pmwg:::extend_sampler(pmwgs, iter, stage)
# create progress bar
start_iter <- pmwgs$samples$idx

collected_msgs <- list()

# Main iteration loop
i <- 1
# Create/update efficient proposal distribution if we are in sampling phase.
stable_args <- utils::modifyList(
  stable_args,
  pmwg:::set_proposal(i, stage, pmwgs, pdist_update_n)
)

# Problem in Gibbs step -> going down one level
sampler <- pmwgs

  # Get single iter versions, tmu = theta_mu, tsig = theta_sig
  last <- pmwg:::last_sample(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  # Here mu is group mean, so we are getting mean and variance
  var_mu <- MASS::ginv(
    sampler$n_subjects * last$tsinv + prior$theta_mu_invar
  )
  levelplot(var_mu)
  mean_mu <- as.vector(
    var_mu %*% (last$tsinv %*% apply(last$alpha, 1, sum) +
                prior$theta_mu_invar %*% prior$theta_mu_mean)
  )
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  levelplot(chol_var_mu)
  levelplot(last$tsinv)
  corrplot::corrplot(last$tsig, is.corr=FALSE, method="color", tl.col="black", cl.pos='b', cl.length=3)
  # New sample for mu.
  tmu <- mvtnorm::rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names

  # New values for group var
  theta_temp <- last$alpha - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
  corrplot::corrplot(B_half, is.corr=FALSE, method="color", tl.col="black", cl.pos='b', cl.length=3)
  tsig <- MCMCpack::riwish(hyper$k_half, B_half) # New sample for group variance  *********HERE IS THE FAILURE****************
  tsinv <- MASS::ginv(tsig)

  # Sample new mixing weights.
  a_half <- 1 / stats::rgamma(
    n = sampler$n_pars,
    shape = hyper$v_shape,
    scale = 1 / (hyper$v_half * diag(tsinv) + hyper$A_half)
  )
  list(
    tmu = tmu,
    tsig = tsig,
    tsinv = tsinv,
    a_half = a_half,
    alpha = last$alpha
  )



tryCatch(
  pars <- pmwg:::gibbs_step(pmwgs),
  error = function(err_cond) {
    pmwg:::gibbs_step_err(pmwgs, err_cond)
  }
)

  iter_args <- list(
    parameters = pars,
    epsilon = subj_epsilon
  )

  tmp <- do.call(apply_fn, c(stable_args, iter_args))
  lapply(tmp, function(x) {
    if (inherits(x, "try-error")) {
      cat("ERROR: At least 1 call to log likelihood method caused an error\n")
      traceback(x)
      stop()
    }
  })

  ll <- unlist(lapply(tmp, attr, "ll"))
  alpha <- array(unlist(tmp), dim = dim(pars$alpha))

  # Store results locally.
  j <- start_iter + i
  pmwgs$samples$theta_mu[, j] <- pars$tmu
  pmwgs$samples$theta_sig[, , j] <- pars$tsig
  pmwgs$samples$last_theta_sig_inv <- pars$tsinv
  pmwgs$samples$alpha[, , j] <- alpha
  pmwgs$samples$idx <- j
  pmwgs$samples$subj_ll[, j] <- ll
  pmwgs$samples$a_half[, j] <- pars$a_half
  pmwgs$samples$epsilon[, j] <- subj_epsilon

  # Epsilon tuning. See Garthwaite, P. H., Fan, Y., & Sisson, S. A. (2016).
  if (!is.null(p_accept)) {
    if (j > n0) {
      acc <-  pmwgs$samples$alpha[1, , j] != pmwgs$samples$alpha[1, , (j - 1)]
      subj_epsilon <- update_epsilon(subj_epsilon, acc, p_accept, j,
                                     pmwgs$n_pars, alpha_star)
    }
  }

  if (stage == "adapt") {
    res <- test_sampler_adapted(pmwgs, n_unique, i)
    if (res[1] == "success") {
      collected_msgs <- c(collected_msgs, res[2])
      break
    } else if (res[1] == "increase") {
      n_unique <- n_unique + .unique_inc
      collected_msgs <- c(collected_msgs, res[2])
    }
  }
}
if (length(collected_msgs) > 0) lapply(collected_msgs, cat)
if (stage == "adapt") {
  if (i == iter) {
    message(paste(
      "WARNING:",
      "Particle Metropolis within Gibbs Sampler did not",
      "finish adaptation phase early (all", i, "iterations were",
      "run).\nYou should examine your samples and perhaps start",
      "a longer adaptation run."
    ))
  } else {
    pmwgs <- trim_na(pmwgs)
  }
}
pmwgs

