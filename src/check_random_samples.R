library(dplyr)
library(ggplot2)

recovery_data <- "NumericVDCE_IEX_jDB1CRlMuY1Z_untagged_data.RDS"
real_data <- "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData"

get_data <- function(f, obj_name = "sampled") {
  # Load in the data into the global environment
  load(f, envir = (e <- new.env()))
  e[[obj_name]]$data
}

recovery <- readRDS(file = here::here("data", "output", recovery_data))
real <- get_data(here::here("data", "output", real_data))

combined <- bind_rows(real=real, recovery=recovery, .id="source")

combined %>%
  mutate(cell=paste0(price, rating)) %>%
  ggplot(aes(x = rt, colour = source)) +
  geom_density() +
  facet_wrap(~ cell)
