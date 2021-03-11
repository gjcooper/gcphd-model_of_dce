source("read_data.R")
library(BayesFactor)
library(RColorBrewer)
library(varBF)

raw_data <- get_data()
messages <- lapply(raw_data, "[[", "messages")
events <- extract_data(raw_data)
events <- augment_data(events)
sections <- filter_data(events)
list2env(sections, environment())
survey <- get_survey(events)



write_to_file <- function() {
  write.csv(events, file = "all_data.csv")
  write.csv(complete, file = "complete_data.csv")
  write.csv(implicit, file = "implicit_data.csv")
  write.csv(explicit, file = "explicit_data.csv")
  write.csv(sft, file = "sft_data.csv")
  write.csv(survey, file = "survey_data.csv")

  save(complete, explicit, implicit, sft, survey, file = "psft2020.RData")
}

message("The following datasets were incomplete")
print(sft %>% group_by(sonaID) %>% count() %>% arrange(n) %>% filter(n<720))

# Show splits by conditions
sft %>%
  group_by(sonaID) %>%
  summarise(order = unique(attr_order), group = unique(group), hand = unique(hand)) %>%
  group_by(order, group, hand) %>%
  summarise(n = n()) %>%
  spread(group, n) %>%
  kable()

sft %>%
  distinct(sonaID) %>%
  count(name="Number of Datasets")


#Replicate Rmd

explicit.thresholds <- explicit%>% 
  group_by(sonaID) %>% 
  summarise(price=mean(price_thresh),
            rating=mean(rating_thresh))
implicit.thresholds <- implicit%>% 
  group_by(sonaID) %>% 
  summarise(price=unique(price_thresh),
            rating=unique(rating_thresh))

thresholds <- bind_rows(explicit.thresholds, implicit.thresholds, .id="condition")
thresholds$condition <- thresholds$condition %>% factor
levels(thresholds$condition) <- c("explicit", "implicit")

thresholds.summary <- thresholds %>% 
  group_by(condition) %>% 
  summarise(price.mean=mean(price), price.SD=sd(price),
            rating.mean=mean(rating), rating.SD=sd(rating))

# Price threshold plot
thresholds %>%
  ggplot(mapping = aes(x = sonaID, y = price, colour = condition)) +
  geom_point(size=4) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Participants") +
  ylab("Price Threshold") +
  geom_hline(mapping = aes(yintercept = price.mean, colour = condition), data = thresholds.summary, size=4) +
  scale_colour_brewer(palette="Set2")


ttestBF(explicit.thresholds$price, implicit.thresholds$price)
ttestBF(explicit.thresholds$rating, implicit.thresholds$rating)

indepvarBF(explicit.thresholds$price, implicit.thresholds$price)
indepvarBF(explicit.thresholds$rating, implicit.thresholds$rating)

explicit.times <- explicit %>% 
  group_by(sonaID) %>% 
  summarise(time=sum(rt))

implicit.times <- implicit %>% 
  group_by(sonaID) %>% 
  summarise(time=sum(rt), n.trials=length(rt))

implicit.times$n.trials %>% 
  sort %>% 
  plot(main="Implicit condition (Phase 1)", ylab="Number of staircase trials", xlab="Sorted participants")

implicit.times %>%
  summarise(median=median(n.trials), mean=mean(n.trials), 
            SD=sd(n.trials), IQR=IQR(n.trials), 
            minimum=min(n.trials), maximum=max(n.trials))

ttestBF(explicit.times$time, implicit.times$time)

bind_rows(explicit = explicit.times, implicit = implicit.times, .id = "condition") %>%
  mutate(condition = factor(condition, levels = c("explicit", "implicit"))) %>%
  group_by(condition) %>%
  ggplot(mapping = aes(x = sonaID, y = time, colour = condition)) +
  geom_point(size = 4) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Participants") +
  ylab("Time to Determine Threshold") +
  scale_colour_brewer(palette="Set2")


times <- bind_rows(explicit.times, implicit.times, .id="condition")
times$condition <- times$condition %>% factor
levels(times$condition) <- c("explicit", "implicit")

times %>% 
  group_by(condition) %>% 
  summarise(time.mean=mean(time), time.SD=sd(time))

choice.data <- sft

n.total.trials <- 720
choice.data %>% 
  count(sonaID) %>%
  pull(n) %>% 
  sort %>%
  plot(ylab="Number of trials completed", xlab="Sorted participants")
abline(h=n.total.trials * c(.25, .5, .75))

n.total.subjects <- choice.data$sonaID %>% unique %>% length
n.acceptable.trials <- .75 * n.total.trials
choice.data <- choice.data %>% 
  group_by(sonaID) %>% 
  filter(length(rt) >= n.acceptable.trials)
n.subjects <- choice.data$sonaID %>% unique %>% length
print(c(n.total.subjects, n.subjects))

choice.data <- choice.data %>% filter(!is.na(rt))
lower.cutoff <- .25
too.fast <- choice.data %>% 
  group_by(sonaID) %>% 
  summarise(y=sum(rt<lower.cutoff)/length(rt)*100)
too.fast$y %>%
  sort %>%
  plot(ylab="Percentage of very fast RTs", xlab="Sorted participants")
abline(h=20)

keep.subjects <- too.fast %>% 
  filter(y < 20) %>%
  pull(sonaID)
choice.data <- choice.data %>% 
  filter(sonaID %in% keep.subjects)
print(c(nrow(too.fast), length(keep.subjects)))

choice.data$rt %>% 
  sort %>% 
  plot(ylab="RT (s)", xlab="Sorted trials")

choice.data$rt %>% 
  sort %>%
  head(5000) %>%
  plot(ylab="RT (s)", xlab="Sorted trials")
  abline(h=lower.cutoff, v=1300)

upper.cutoff <- 60
choice.data$rt %>% 
  sort %>%
  tail(500) %>%
  plot(ylab="RT (s)", xlab="Sorted trials", ylim=c(0,200))
abline(h=upper.cutoff)

n.trials <- nrow(choice.data)
choice.data <- choice.data %>%
  filter(between(rt, lower.cutoff, upper.cutoff))
choice.data$rt %>% range
n.trials - nrow(choice.data)
(n.trials - nrow(choice.data)) / n.trials * 100

good.hotels <- c("HH","HL","LH","LL")
bad.hotels <- c("DD","DL","DH","LD","HD")
choice.data$correct <- NA
choice.data$correct[(choice.data$cell_name %in% good.hotels) & (choice.data$accept == TRUE)] <- TRUE
choice.data$correct[(choice.data$cell_name %in% good.hotels) & (choice.data$accept == FALSE)] <- FALSE
choice.data$correct[(choice.data$cell_name %in% bad.hotels) & (choice.data$accept == FALSE)] <- TRUE
choice.data$correct[(choice.data$cell_name %in% bad.hotels) & (choice.data$accept == TRUE)] <- FALSE

summary.data <- choice.data %>%
  group_by(sonaID, group, cell_name) %>%
  summarise(cor=mean(correct))

summary.data$cell_name <- factor(summary.data$cell_name, c(good.hotels, bad.hotels))

summary.data %>% 
  group_by(group, cell_name) %>% 
  summarise(correct=mean(cor)) %>% 
  ggplot(map=aes(x = cell_name, y = correct, fill = group)) +
  geom_col(position="dodge") + ylab("Consistency") + xlab("Cell") +
  scale_fill_brewer(palette="Set2")

  barplot(correct ~ group + cell_name, data=., ylim=c(0,1), beside=TRUE, legend=TRUE,
          ylab="Consistency", xlab="Cell")

summary.data$sonaID <- summary.data$sonaID %>% factor
summary.data$price_level <- summary.data$cell_name %>%
  substring(1,1) %>% 
  factor(levels=c("D","L","H"))
summary.data$rating_level <- summary.data$cell_name %>%
  substring(2,2) %>% 
  factor(levels=c("D","L","H"))

summary.data %>% 
  select(sonaID, group, price_level, rating_level, cor) %>%
  pivot_wider(names_from=price_level, values_from=cor, names_pref="price_") %>%
  pivot_wider(names_from=rating_level, values_from=c(price_D, price_H, price_L), names_pref="rating_")

m1 <- anovaBF(cor ~ group * price_level * rating_level + sonaID, whichRandom="sonaID", data=summary.data)
m1 %>% summary
