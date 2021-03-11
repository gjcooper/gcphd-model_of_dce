library(tidyverse)
library(ggforce)
# Load functions and imports from read_data.R file
raw_data <- get_data()
messages <- lapply(raw_data, "[[", "messages")
events <- extract_data(raw_data)
events <- augment_data(events)
events %>% select(price, rating, type, staircase_point) %>% filter(type == "staircase")
rel <- function(point) {
  min_price = 100
  max_price = 450
  min_rating = 6.5
  max_rating = 9
  rate_i = 9.73
  rate_s = -319.06
  price = min_price + point * (max_price - min_price)
  inv_rating = rate_i + rate_s * (1 / price)
  lin_rating = min_rating + point * (max_rating - min_rating)
  list(price=price, inv_rating = inv_rating, lin_rating = lin_rating)
}
diff(sapply( seq(from=0.2, to=1, by=0.02), FUN = function(val) { rel(val)$inv_rating}))
diff(sapply( seq(from=0.2, to=1, by=0.01), FUN = function(val) { rel(val)$inv_rating}))
plot(x=seq(from=0.2, to=.99, by=0.01), y=diff(sapply( seq(from=0.2, to=1, by=0.01), FUN = function(val) { rel(val)$inv_rating})), 'l', xlab="Raw Point", ylab="Rating change for raw change of 0.01")
dev.copy(png, "Rating_Change_0.01.png")
dev.off()
plot(x=seq(from=0.2, to=.98, by=0.02), y=diff(sapply( seq(from=0.2, to=1, by=0.02), FUN = function(val) { rel(val)$inv_rating})), 'l', xlab="Raw Point", ylab="Rating change for raw change of 0.02")
dev.copy(png, "Rating_Change_0.02.png")
dev.off()

events %>% select(price_thresh, rating_thresh, type, response) %>% filter(type == "slider")
events %>% select(price_threshold_rescaled, rating_threshold_rescaled, type, response) %>% filter(type == "slider")
H = c(0.15, 0.025, 0.1, 0.2)
L = c(0.05, 0.025, 0, 0.1)
D = c(-0.1, 0.025, -0.05, -0.1)

#Getting the vals object
vals <- events %>% select(sonaID, price, rating, price_thresh, rating_thresh, type, cell_name) %>% filter(type=="sft") %>% separate(cell_name, into=c("PriceSalience", "RatingSalience"), sep=1)

#Looking at all old sampled prices
gg = ggplot(vals, mapping=aes(x=price, color=PriceSalience)) + geom_histogram(binwidth=1) + geom_vline(aes(xintercept = price_thresh)) + facet_wrap_paginate(~ sonaID, nrow=3, ncol=3, page=1)
ggplot(vals, mapping=aes(x=price, color=PriceSalience)) + geom_histogram(binwidth=1) + geom_vline(aes(xintercept = price_thresh))
n_pages(gg)

# Looking at old sampled prices (to match new)
vals %>%
 group_by(price_thresh) %>%
 filter(n() == 720) %>%
 group_by(price_thresh) %>%
 filter(price_thresh %in% c(199.75, 415)) %>%
 ggplot(mapping=aes(x=price, color=PriceSalience)) +
   geom_histogram(binwidth=1) +
   geom_vline(aes(xintercept = price_thresh)) +
   facet_wrap(~ price_thresh)


#Looking at new sampled prices
vals %>%
 group_by(price_thresh) %>%
 filter(n() == 720) %>%
 ggplot(mapping=aes(x=price, color=PriceSalience)) +
   geom_histogram(binwidth=1) +
   geom_vline(aes(xintercept = price_thresh)) +
   facet_wrap(~ price_thresh)

pdf("AllPrices.pdf")
for (i in 1:n_pages(gg)) {
print(gg + facet_wrap_paginate(~ sonaID, nrow=3, ncol=3, page=i))
}
dev.off()
