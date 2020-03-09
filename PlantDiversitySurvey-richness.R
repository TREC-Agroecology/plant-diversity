### This script organizes and evaluates the NEON Plant Diversity Survey
### conducted at TREC in Homestead, FL Fall 2017 - Spring 2018.
### Working directory is 'Source File Location'.


library(tidyverse)
library(agricolae)


## Data

surveys <- read_csv("data/PlantDiversitySurvey-surveys.csv")
status <- read_csv("data/TRECstatus.csv")
invasive <- read_csv("data/TRECinvasive.csv")

surveys_w_plot <- surveys %>%
  filter(!is.na(genus)) %>%
  mutate(genus_species = paste(tolower(str_extract(genus, "...")),
                               tolower(str_trunc(species, 3, "right", "")), 
                               sep="")) %>%
  separate(code, c("big_plot", "corner", "small_plot"), sep="\\.")
     # Expect missing pieces for 100m2 plots
surveys_w_plots <- surveys_w_plot %>% 
  select(block, site, big_plot, corner, small_plot, genus_species)

echo <- read.csv("data/ECHO-surveys.csv") %>% 
  filter(UncertainId %in% c(NA, "Species")) %>%
  mutate(genus_species = tolower(genus_species)) %>%
  separate(code, c("big_plot", "corner", "small_plot"), sep="\\.") %>%
  select(block, site, big_plot, corner, small_plot, genus_species)

surveys_w_plots <- bind_rows(surveys_w_plots, echo)

pub_sites <- data.frame(block = c(1, 4, 14, 15, 31, 32), 
                     pub_site = c("TREC-NW", "TREC-NE", "TREC-SW", "TREC-SE", "ECHO-E", "ECHO-W"))

cluster <- data.frame(block = c(1, 1, 4, 4, 14, 14, 15, 15),
                      site = rep(c("N", "S"), 4),
                      status = c("high", "high", "low", "high",
                                 "high", "low", "low", "low"),
                      cluster = c("open", "open", "lawn", "open",
                                  "open", "hammock", "orchard", "orchard"))

## Scale Summaries

ones <- surveys_w_plots %>%
  filter(small_plot == 1) %>%
  group_by(block, site, big_plot, corner) %>%
  summarize(ones = n())

tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(genus_species) %>%
  summarize(tens = n())

big_plots <- surveys_w_plots %>%  ## Big_plots = 100m
  group_by(block, site, big_plot) %>%
  distinct(genus_species) %>%
  summarize(hundreds = n())

big_plots <- mutate(big_plots, site_stat = paste(block, site, sep=""))

sites <- surveys_w_plots %>%
  group_by(block, site) %>%
  distinct(genus_species) %>%
  summarize(richness = n())

blocks <- surveys_w_plots %>%
  group_by(block) %>%
  distinct(genus_species) %>%
  summarize(richness = n())


## Combined Summaries

record_counts <- ones %>% 
  right_join(tens) %>% 
  inner_join(big_plots) %>% 
  mutate(site_stat = paste(block, site, sep="")) %>% 
  left_join(cluster)

record_counts_sites <- ones %>% 
  right_join(tens) %>% 
  inner_join(big_plots) %>% 
  mutate(site_stat = paste(block, site, sep="")) %>% 
  left_join(pub_sites)


## Calculate and visualize richness among scales [[TREC ONLY]]

ggplot(tens, aes(x=tens, fill=site)) +  # tens visualization
  geom_histogram(binwidth = 5) +
  facet_grid(.~block) +
  labs(x="Richness", y="Count", fill="Site") +
  theme_bw(base_size=24, base_family="Helvetica")
ggsave("output/richness_ten.png", width = 12, height = 5)  

site_avg <- record_counts %>%
  group_by(block, site) %>%
  summarize(avg_ones = round(mean(ones, na.rm=TRUE), 0), sd_ones = round(sd(ones), 2),
            avg_tens = round(mean(tens), 0), sd_tens = round(sd(tens), 2),
            avg_hund = round(mean(hundreds), 0), sd_hund = round(sd(hundreds), 2))

site_avg_plot <- record_counts %>%
  group_by(status, cluster) %>%
  summarize("1" = round(mean(ones, na.rm=TRUE), 0),
            "10" = round(mean(tens), 0),
            "100" = round(mean(hundreds), 0)) %>% 
  filter(!is.na(cluster)) %>% 
  gather(scale, average, "1":"100")

ggplot(site_avg_plot, aes(x=scale, y=average, group=cluster)) +
  geom_line() +
  geom_point(size=3, alpha = 0.8,  aes(shape=status, color=as.factor(cluster))) +
  labs(x="Scale [m2]", y="Average Richness", shape="Soil Disturbance", color="Habitat") +
  theme_classic(base_size=14, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))
ggsave("output/avg_richness.png", width = 5, height = 4)  

## ANOVA and Tukey Post Hoc

sink("output/richness.txt")

hundreds_aov <- aov(hundreds ~ block + site_stat, data=big_plots)
cat("hundreds ~ block + site_stat\n") 
print(summary(hundreds_aov))
cat("\nblock\n")
HSD.test(hundreds_aov, "block")$groups
cat("\nsite\n")
HSD.test(hundreds_aov, "site_stat")$groups

tens_aov <- aov(tens ~ block + site_stat, data=record_counts)
cat("\ntens ~ block + site_stat\n") 
print(summary(tens_aov))
cat("\nblock\n")
HSD.test(tens_aov, "block")$groups
cat("\nsite\n")
HSD.test(tens_aov, "site_stat")$groups

ones_aov <- aov(ones ~ block + site_stat, data=record_counts)
cat("\nones ~ block + site_stat\n") 
summary(ones_aov)
cat("\nblock\n")
HSD.test(ones_aov, "block")$groups
cat("\nsite\n")
HSD.test(ones_aov, "site_stat")$groups

sink()

## Status Counts

block_status <- surveys_w_plot %>%
  left_join(cluster) %>% 
  group_by(cluster, block, site) %>%
  distinct(genus, species) %>% 
  left_join(status) %>% 
  left_join(invasive) %>% 
  mutate(assessment = str_replace(assessment, ".*", "1")) %>%
  summarize(richness = n(), natives = sum(native, na.rm=TRUE),
            established = sum(established_FL, na.rm=TRUE) - sum(native, na.rm=TRUE))

## Species Area Curves [[TREC ONLY]]

richness_results <- record_counts %>% 
  ungroup() %>% 
  select(cluster, "1" = ones, "10" = tens, "100" = hundreds) %>%
  filter(!is.na(cluster)) %>% 
  gather(scale, richness, "1":"100")

sink("output/species_area.txt")
for (habitat in unique(cluster$cluster)) {
  habitat_results <- filter(richness_results, cluster == habitat)
  test <- lm(log10(richness) ~ log10(as.numeric(scale)), data=habitat_results)
  print(paste(habitat, round(summary(test)$r.squared, 3), 
              round(summary(test)$coeff[2,1]), 3))
}
sink()