library(coda)
library(dplyr)
library(tidyr)
require(tcltk)

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

mcobj <- mcmc(t(plot_obj$samples$theta_mu))
dimnames(mcobj) <- list(NULL,parameters)
plot(mcobj, smooth=TRUE)


# hd <- heidel.diag(sobj)
# plot(sobj)
# summary(sobj)$statistics

# plot(sobj[,c('alpha_IST', 'alpha_IEX', 'A', 'b_pos', 'b_neg')])


# means$run <- 1:nrow(means)
# me2 <- pivot_longer(means, -run)
# ggplot(me2) + geom_point()
# ggplot(me2) + geom_point(aes(x=name, y=value))
# ggplot(me2) + geom_point(aes(x=name, y=value, colour=run))
# ggplot(me2) + geom_point(aes(x=name, y=value, colour=run)) + scale_colour_discrete()
# ggplot(me2) + geom_point(aes(x=name, y=value, colour=run)) + scale_colour_continuous()
# ggplot(me2) + geom_point(aes(x=name, y=value, colour=run)) + scale_colour_brewer()
# table(me2$run)
# ggplot(me2) + geom_point(aes(x=name, y=value, colour=run)) + scale_colour_manual()
# me2$run <- as.factor(me2$run)
# ggplot(me2) + geom_point(aes(x=name, y=value, colour=run))
# me2
# me2[me2$name != alpha_IST]
# me2[me2$name != 'alpha_IST']
# me2[me2$name != 'alpha_IST',]
# me3 = me2[me2$name != 'alpha_IST',]
# ggplot(me3) + geom_point(aes(x=name, y=value, colour=run))
# ggplot(me2) + geom_point(aes(x=name, y=value, colour=run))
# history
# history()
