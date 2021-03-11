library(tidyverse)

sft.df <- read_csv('sft_data.csv')

sft.df <- sft.df %>%
  mutate(consistent = case_when(
    !str_detect(cell_name, "D") & accept ~ TRUE,
    str_detect(cell_name, "D") & !accept ~ TRUE,
    TRUE ~ FALSE
    )
  )

png(file="consistent.png", width=800, height=432)
g <- sft.df %>%
  group_by(sonaID) %>%
  summarise(pccorr=mean(consistent), group=first(group)) %>%
  mutate(sonaID = factor(sonaID)) %>%
  ggplot(aes(x=reorder(sonaID, pccorr), y=pccorr)) +
    geom_col(aes(fill=group)) +
    xlab("Participants") +
    ylab("Proportion of consistent responses") +
    geom_hline(yintercept=0.5) +
    theme_classic() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )

print(g)
dev.off()

# get_SIC <- function(sid) {
#   df <- sft.df[sft.df$sonaID == sid, ]
# 
# explicit <- read_csv("explicit_data.csv")
# 
# explicit <- explicit %>%
#   group_by(sonaID) %>%
#   mutate(trial_order = ifelse(row_number() == 1, "first", "second"))
# 
# explicit %>%
#   ggplot(aes(x = factor(sonaID), y = response, colour = factor(trial_order))) +
#     geom_point()
