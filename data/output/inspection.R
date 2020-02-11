library(coda)
library(dplyr)
library(tidyr)
# load("cceburn5kfull.RData")
# burned1 <- burned
# load("cceburn5kfull-2.RData")
# burned2 <- burned
# load("cceburn5kfull-3.RData")
# burned3 <- burned
# load("cceadapt.RData")
source("~/PMwG/demo/plot.R")

load('cceadaptmore.RData')

load('cceburn5kfullN1.RData')
burned1 <- burned
load('cceburn5kfullN2.RData')
burned2 <- burned
load('cceburn5kfullN3.RData')
burned3 <- burned

load('ccesample.RData')

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
knitr::kable(table(plot_obj$samples$stage))

plot_obj <- sample

plot(plot_obj, pars=alphas)
plot(plot_obj, type='alpha', pars=alphas)
plot(plot_obj, pars = base)
plot(plot_obj, type='alpha', pars=base)
plot(plot_obj, pars = fourcell)
plot(plot_obj, type='alpha', pars=fourcell)
plot(plot_obj, pars=posdrift)
plot(plot_obj, type='alpha', pars=posdrift)
plot(plot_obj, pars=negdrift)
plot(plot_obj, type='alpha', pars=negdrift)

mcobj <- mcmc(t(plot_obj$samples$theta_mu))

sobj <- mcmc(t(exp(
  plot_obj$samples$theta_mu[,plot_obj$samples$stage == 'sample']
  )))

dimnames(sobj) <- list(NULL,parameters)

hd <- heidel.diag(sobj)
plot(sobj)
summary(sobj)$statistics

plot(sobj[,c('alpha_IST', 'alpha_IEX', 'A', 'b_pos', 'b_neg')])


means$run <- 1:nrow(means)
me2 <- pivot_longer(means, -run)
ggplot(me2) + geom_point()
ggplot(me2) + geom_point(aes(x=name, y=value))
ggplot(me2) + geom_point(aes(x=name, y=value, colour=run))
ggplot(me2) + geom_point(aes(x=name, y=value, colour=run)) + scale_colour_discrete()
ggplot(me2) + geom_point(aes(x=name, y=value, colour=run)) + scale_colour_continuous()
ggplot(me2) + geom_point(aes(x=name, y=value, colour=run)) + scale_colour_brewer()
table(me2$run)
ggplot(me2) + geom_point(aes(x=name, y=value, colour=run)) + scale_colour_manual()
me2$run <- as.factor(me2$run)
ggplot(me2) + geom_point(aes(x=name, y=value, colour=run))
me2
me2[me2$name != alpha_IST]
me2[me2$name != 'alpha_IST']
me2[me2$name != 'alpha_IST',]
me3 = me2[me2$name != 'alpha_IST',]
ggplot(me3) + geom_point(aes(x=name, y=value, colour=run))
ggplot(me2) + geom_point(aes(x=name, y=value, colour=run))
history
history()
