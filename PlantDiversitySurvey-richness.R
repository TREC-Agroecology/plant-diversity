### This script organizes and evaluates the NEON Plant Diversity Survey
### conducted at TREC in Homestead, FL Fall 2017 - Spring 2018.
### Working directory is 'Source File Location'.


library(tidyverse)
library(agricolae)


## Data

surveys <- read_csv("data/PlantDiversitySurvey-surveys.csv")
surveys_w_plots <- surveys %>%
  separate(code, c("big_plot", "corner", "small_plot"), sep="\\.")


## Scale Summaries

ones <- surveys_w_plots %>%
  filter(small_plot == 1) %>%
  group_by(block, site, big_plot, corner) %>%
  summarize(ones = n())

tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(taxonID, taxonIDRemarks) %>%
  summarize(tens = n())

big_plots <- surveys_w_plots %>%  ## Big_plots = 100m
  group_by(block, site, big_plot) %>%
  distinct(taxonID, taxonIDRemarks) %>%
  summarize(hundreds = n())

big_plots <- mutate(big_plots, site_stat = paste(block, site, sep=""))

sites <- surveys %>%
  group_by(block, site) %>%
  distinct(taxonID, taxonIDRemarks) %>%
  summarize(records = n())

blocks <- surveys %>%
  group_by(block) %>%
  distinct(taxonID, taxonIDRemarks) %>%
  summarize(records = n())


## Combined Summaries

record_counts <- right_join(ones, tens)
record_counts <- inner_join(record_counts, big_plots)
record_counts <- mutate(record_counts, site_stat = paste(block, site, sep=""))


## Calculate and visualize average richness among scales

site_avg <- record_counts %>%
  group_by(block, site) %>%
  summarize(avg_ones = round(mean(ones), 0), sd_ones = round(sd(ones), 2),
            avg_tens = round(mean(tens), 0), sd_tens = round(sd(tens), 2),
            avg_hund = round(mean(hundreds), 0), sd_hund = round(sd(hundreds), 2))

site_avg_plot <- record_counts %>%
  group_by(block, site) %>%
  summarize("1" = round(mean(ones), 0),
            "10" = round(mean(tens), 0),
            "100" = round(mean(hundreds), 0))
site_avg_plot <- gather(site_avg_plot, scale, average, "1":"100")
site_avg_plot <- mutate(site_avg_plot, site_stat = paste(block, site, sep=""))

ggplot(site_avg_plot, aes(x=scale, y=average, group=site_stat)) +
  geom_line() +
  geom_point(size=3, aes(color=site_stat)) +
  labs(x="Scale [m2]", y="Average Record Count", color="Site") +
  theme_classic()
  

## ANOVA

hundreds_aov <- aov(hundreds ~ block + site_stat, data=big_plots)
summary(hundreds_aov)

tens_aov <- aov(tens ~ block + site_stat, data=record_counts)
summary(tens_aov)

ones_aov <- aov(ones ~ block + site_stat, data=record_counts)
summary(ones_aov)


## Tukey Post Hoc

HSD.test(hundreds_aov, "block")$groups
HSD.test(tens_aov, "block")$groups
HSD.test(ones_aov, "block")$groups
HSD.test(hundreds_aov, "site_stat")$groups
HSD.test(tens_aov, "site_stat")$groups
HSD.test(ones_aov, "site_stat")$groups
