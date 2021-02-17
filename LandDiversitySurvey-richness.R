### This script organizes and evaluates the NEON Plant Diversity Survey
### conducted at TREC in Homestead, FL Spring 2020 with land history.
### Working directory is 'Source File Location'.


library(tidyverse)
library(agricolae)


## Data

surveys <- read_csv("data/NEONsurveyMar2020PlotData.csv")
status <- read_csv("data/TRECstatus.csv")
invasive <- read_csv("data/TRECinvasive.csv")

surveys_w_plot <- surveys %>%
  filter(Certainty == "ToSpp") %>%
  separate(Subplot, c("big_plot", "corner", "small_plot"), sep="\\.")
     # Expect missing pieces for 100m2 plots
surveys_w_plots <- surveys_w_plot %>% 
  select(block = Block, site = Land, big_plot, corner, small_plot, 
         genus_species = Code_CGM)

surveys_w_plots$block <- str_replace(surveys_w_plots$block, "B02", "B01")

pub_blocks <- data.frame(block = c("B01", "B04", "B14", "B15"), 
                      pub_block = factor(c("Old-Ag-NW", "Old-Ag-NE", "New-Ag-SW", "New-Ag-SE"), 
                                         levels = c("Old-Ag-NW", "Old-Ag-NE", "New-Ag-SW", "New-Ag-SE")))

pub_sites <- data.frame(site = c("cc", "gr"), pub_site = c("Ag", "Lawn"))

## Scale Summaries

ones <- surveys_w_plots %>%
  filter(small_plot == 1) %>%
  group_by(block, site, big_plot, corner) %>%
  summarize(ones = n())

tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(genus_species) %>%
  summarize(tens = n()) %>% 
  left_join(pub_sites) %>% 
  left_join(pub_blocks)

big_plots <- surveys_w_plots %>%  ## Big_plots = 100m
  group_by(block, site, big_plot) %>%
  distinct(genus_species) %>%
  summarize(hundreds = n())

big_plots <- mutate(big_plots, site_stat = paste(block, site, sep=""))

sites <- surveys_w_plots %>%
  group_by(block, site) %>%
  distinct(genus_species) %>%
  summarize(richness = n())
write_csv(sites, "output/block_site_total_richness.csv")

blocks <- surveys_w_plots %>%
  group_by(block) %>%
  distinct(genus_species) %>%
  summarize(richness = n())


## Combined Summaries

record_counts <- ones %>% 
  right_join(tens) %>% 
  inner_join(big_plots) %>% 
  mutate(site_stat = paste(block, site, sep=""))

record_counts_sites <- ones %>% 
  right_join(tens) %>% 
  inner_join(big_plots) %>% 
  mutate(site_stat = paste(block, site, sep="")) %>% 
  left_join(pub_sites) %>% 
  left_join(pub_blocks)


## Calculate and visualize richness among scales

ggplot(tens, aes(x=tens, fill=pub_site)) +  # tens visualization
  geom_histogram(binwidth = 5) +
  facet_grid(.~pub_block) +
  labs(x="Richness", y="Count", fill="Site") +
  theme_bw(base_size=24, base_family="Helvetica")
ggsave("output/2020_richness_ten.png", width = 12, height = 5)  

site_avg <- record_counts %>%
  group_by(block, site) %>%
  summarize(avg_ones = round(mean(ones, na.rm=TRUE), 0), sd_ones = round(sd(ones), 2),
            avg_tens = round(mean(tens), 0), sd_tens = round(sd(tens), 2),
            avg_hund = round(mean(hundreds), 0), sd_hund = round(sd(hundreds), 2)) 

site_avg_plot <- record_counts %>%
  group_by(block, site, site_stat) %>%
  summarize("1" = round(mean(ones, na.rm=TRUE), 0),
            "10" = round(mean(tens), 0),
            "100" = round(mean(hundreds), 0)) %>% 
  gather(scale, average, "1":"100") %>% 
  left_join(pub_sites) %>% 
  left_join(pub_blocks)

ggplot(site_avg_plot, aes(x=scale, y=average, group=site_stat)) +
  geom_line() +
  geom_point(size=3, alpha = 0.8,  aes(shape=pub_site, color=pub_block)) +
  labs(x="Scale [m2]", y="Average Richness", shape="Site", color="Block") +
  theme_classic(base_size=14, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))
ggsave("output/2020_avg_richness.png", width = 5, height = 4)  

## ANOVA and Tukey Post Hoc

sink("output/2020_richness.txt")

hundreds_aov <- aov(hundreds ~ block + site + block*site, data=big_plots)
cat("hundreds ~ block + site + block*site\n") 
print(summary(hundreds_aov))
cat("\nblock\n")
HSD.test(hundreds_aov, "block")$groups
cat("\nsite\n")
HSD.test(hundreds_aov, "site")$groups

tens_aov <- aov(tens ~ block + site + block*site, data=record_counts)
cat("\ntens ~ block + site + block*site\n") 
print(summary(tens_aov))
cat("\nblock\n")
HSD.test(tens_aov, "block")$groups
cat("\nsite\n")
HSD.test(tens_aov, "site")$groups

ones_aov <- aov(ones ~ block + site + block*site, data=record_counts)
cat("\nones ~ block + site + block*site\n") 
summary(ones_aov)
cat("\nblock\n")
HSD.test(ones_aov, "block")$groups
cat("\nsite\n")
HSD.test(ones_aov, "site")$groups

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

## Species Area Curves

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