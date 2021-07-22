get_data <- function(final_obj = "sampled") {
  # Load in the data into the global environment
  f <- tkgetOpenFile(
    title = "RData file",
    initialdir = here::here("data/output")
  ) %>%
    as.character()

  load(f, envir = (e <- new.env()))

  e[[final_obj]]
}

extract_alphas <- function(plot_obj) {
  tmu_alphas <- plot_obj %>%
    as_mcmc(filter = "sample") %>%
    data.frame() %>%
    tibble() %>%
    select(starts_with("alpha")) %>%
    mutate(subjectid = "theta_mu")

  as_mcmc(plot_obj, selection = "alpha", filter = "sample") %>%
    lapply(FUN = function(x) {
      x %>%
        data.frame() %>%
        tibble() %>%
        select(starts_with("alpha"))
    }) %>%
    bind_rows(.id = "subjectid") %>%
    tibble() %>%
    rbind(tmu_alphas)
}

get_medians <- function(res) {
  res %>%
    group_by(subjectid) %>%
    summarise(across(.fns = median)) %>%
    mutate_at(vars(contains("alpha")), exp) %>%
    pivot_longer(-subjectid) %>%
    rename(Parameter = name) %>%
    group_by(subjectid) %>%
    mutate(allvals = sum(value)) %>%
    mutate(rel_val = value / allvals) %>%
    mutate(Parameter = as.factor(Parameter)) %>%
    mutate(Parameter = fct_relabel(Parameter, ~ sub(".*_", "", .x)))
}

#' Plot the evidence for each model by subject
#'
#' @param plot_obj - the pmwg sampler object containing estimates
#' @param relative - Whether to plot the evidence as relative or absolute
#'
#' @return None - side effect is the creation of the plot
plot_alphas <- function(plot_obj, relative = TRUE) {
  medians <- extract_alphas(plot_obj) %>%
    get_medians()

  if (relative) {
    ggplot(medians, aes(x = subjectid, y = rel_val, fill = Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Relative Evidence") +
      scale_fill_brewer(palette = "Set2", name = "Model")
  } else {
    ggplot(medians, aes(x = subjectid, y = value, fill = Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Absolute Evidence") +
      scale_fill_brewer(palette = "Set2", name = "Model")
  }
}



compare <- function(original, recovery, relative = TRUE) {
  original <- extract_alphas(original) %>%
    get_medians()
  recovery <- extract_alphas(recovery) %>%
    get_medians()

  medians <- bind_rows(
    estimates = original,
    recovery = recovery,
    .id = "source"
  )

  if (relative) {
    ggplot(medians, aes(x = subjectid, y = rel_val, fill = Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Relative Evidence") +
      scale_fill_brewer(palette = "Set2", name = "Model") +
      facet_grid(rows = vars(source))
  } else {
    ggplot(medians, aes(x = subjectid, y = value, fill = Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Absolute Evidence") +
      scale_fill_brewer(palette = "Set2", name = "Model") +
      facet_grid(rows = vars(source))
  }
}

compare_data <- function(original, recovery, comparison = "rt") {
  odata <- original$data %>%
    as_tibble()
  rdata <- recovery$data %>%
    as_tibble()
  alldata <- bind_rows(original = odata, recovery = rdata, .id = "source")
  if (comparison == "rt") {
    alldata %>%
      filter(rt < 5) %>%
      unite("cell", price:rating, sep='') %>%
      ggplot(mapping = aes(x = rt, colour = source)) +
      geom_density() +
      facet_wrap(~cell)
  } else if (comparison == "rtbox") {
    alldata %>%
      filter(rt < 5) %>%
      unite("cell", price:rating, sep='') %>%
      ggplot(mapping = aes(x = cell, y = rt, colour = source)) +
      geom_boxplot()
  } else if (comparison == "response") {
    alldata %>%
      mutate(accept=factor(accept, labels=c('reject', 'accept'))) %>%
      ggplot(mapping = aes(x = accept, group = source, fill = source)) +
      geom_bar(position = "dodge") +
      facet_wrap(~subject)
  }
}

# Look at the estimates in absolute and relative terms, then use coda mcmc plot for theta_mu
recovery <- get_data()
original <- get_data()

(g+ggtitle("Recovery from Coactive")) / (go+ggtitle("Original Data"))
compare(original, recovery)
compare(original, recovery, relative = FALSE)
compare_data(original, recovery, "rt")
compare_data(original, recovery, "rtbox")
compare_data(original, recovery, "response")

# Looking at Task1
ex1mct3 <- get_data()

sample_data <- ex1mct3
mcmcplot(as_mcmc(sample_data, filter = "sample"))
mcmcplot(as_mcmc(sample_data, select = "alpha", filter = "sample"))
sample_df <- sample_data %>%
  as_mcmc() %>%
  as_tibble()

sample_data %>%
  as_mcmc(filter="sample") %>%
  data.frame() %>%
  tibble() %>%
  summarise_all(mean) %>%
  data.frame() %>%
  exp %>%
  round(2)

sample_df %>%
  select(starts_with("v_")) %>%
  pivot_longer(
    everything(),
    names_to = c("drift", "response", "attribute", "salience"),
    names_transform = list(
      response = ~ readr::parse_factor(.x, levels = c("acc", "rej")),
      attribute = ~ readr::parse_factor(.x, levels = c("p", "r")),
      salience = ~ readr::parse_factor(.x, levels = c("H", "L", "D"))
    ),
    names_sep = "_"
  ) %>%
  select(-drift) %>%
  rename(drift = value) %>%
  mutate(
    response = fct_recode(response, Accept = "acc", Reject = "rej"),
    attribute = fct_recode(attribute, Price = "p", Rating = "r")
  ) %>%
  ggplot(mapping = aes(x = salience, y = drift)) +
  geom_boxplot(aes(fill = salience)) +
  facet_grid(vars(response), vars(attribute))

dev.new()
sample_df %>%
  select(starts_with("v_")) %>%
  pivot_longer(
    everything(),
    names_to = c("drift", "response", "attribute", "salience"),
    names_transform = list(
      response = ~ readr::parse_factor(.x, levels = c("acc", "rej")),
      attribute = ~ readr::parse_factor(.x, levels = c("p", "r")),
      salience = ~ readr::parse_factor(.x, levels = c("H", "L", "D"))
    ),
    names_sep = "_"
  ) %>%
  select(-drift) %>%
  rename(drift = value) %>%
  mutate(
    response = fct_recode(response, Accept = "acc", Reject = "rej"),
    attribute = fct_recode(attribute, Price = "p", Rating = "r")
  ) %>%
  ggplot(mapping = aes(x = salience, y = drift)) +
  geom_boxplot(aes(fill = salience)) +
  facet_grid(vars(attribute), vars(response))
