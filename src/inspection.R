library(coda)
require(tidyverse)
require(tcltk)
require(pmwg)

datafile <- as.character(tkgetOpenFile())
load(datafile)

alphas <- parameters[grepl("^a.*", parameters)]
base <- c("A", "b_pos", "b_neg", "t0")
alldrift <- parameters[grepl("^v.*", parameters)]
posdrift <- parameters[grepl("^v_pos.*", parameters)]
negdrift <- parameters[grepl("^v_neg.*", parameters)]
fourcell <- c("v_pos_HH", "v_neg_HH", "v_pos_HL", "v_neg_HL", "v_pos_LH", "v_neg_LH", "v_pos_LL", "v_neg_LL")
# plot(burned2, pars=base)
# dev.new()
# plot(burned3, pars=base)
# plot(burned2, pars=base)
# dev.new()
#base plot(burned3, pars=base)

#Look at phases etc
#knitr::kable(table(adaptmk1$samples$stage))
plot_obj <- sampled
knitr::kable(table(plot_obj$samples$stage))


plot_theta_mu <- function(plot_obj) {
  samples <- as_mcmc(plot_obj, filter="sample")
  dimnames(samples) <- list(NULL,parameters)
  plot(samples, smooth=TRUE)
}

look_at_alphas <- function(plot_obj, relative=TRUE) {
  tmus <- as_mcmc(plot_obj, filter="sample")
  dimnames(tmus) <- list(NULL, parameters)
  tmu_alphas <- tmus[, 1:7] %>% data.frame() %>% tibble()
  tmu_alphas$subjectid <- 'theta_mu'

  res <- as_mcmc(plot_obj, selection="alpha", filter="sample") %>%
    lapply(FUN = function(X) {
      dimnames(X) <- list(NULL, parameters)
      data.frame(X[, 1:7])
    }) %>%
  bind_rows(.id="subjectid") %>%
  tibble() %>%
  rbind(tmu_alphas)


  medians <- res %>% group_by(subjectid) %>% summarise(across(.fns=median)) %>%
    mutate_at(vars(contains("alpha")), exp) %>%
    pivot_longer(-subjectid) %>%
    rename(Parameter=name) %>%
    group_by(subjectid) %>%
    mutate(allvals = sum(value)) %>%
    mutate(relative_evidence = value/allvals) %>%
    mutate(Parameter = as.factor(Parameter))
  
  par_labels = sapply(strsplit(levels(medians$Parameter), "_"), function(x) {x[2]})

  if (relative) {
    ggplot(medians, aes(x=subjectid, y=relative_evidence, fill=Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Relative Evidence") +
      scale_fill_brewer(palette="Set2", name = "Model", labels = par_labels)
  } else {
    ggplot(medians, aes(x=subjectid, y=value, fill=Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Absolute Evidence") +
      scale_fill_brewer(palette="Set2", name = "Model", labels = par_labels)
  }
}

look_at_alphas(plot_obj)
look_at_alphas(plot_obj, relative=FALSE)
plot_theta_mu(plot_obj)
