### This script organizes and evaluates the NEON Plant Diversity Survey
### conducted at pine rocklands in Miami-Dade Summer 2021.
### Working directory is 'Source File Location'.


library(tidyverse)
library(agricolae)

## Data

data_files <- str_subset(list.files("data/"), "_NEON")
for (file in data_files) {
  site_data <- read_csv(paste0("data/", file))
  if (exists("surveys")){
    surveys <- bind_rows(surveys, site_data)
  } else {
    surveys <- site_data
  }
}

status <- read_csv("data/TRECstatus.csv") %>% ## This will need an update.
  mutate(GenusSpecies = paste(tolower(str_extract(genus, "...")),
                        tolower(str_trunc(species, 3, "right", "")), 
                        sep=""))
invasive <- read_csv("data/TRECinvasive.csv") %>% 
  mutate(invasive = 1, 
         GenusSpecies = paste(tolower(str_extract(genus, "...")),
                               tolower(str_trunc(species, 3, "right", "")), 
                               sep=""))

surveys_w_plot <- surveys %>%
  filter(!str_detect(GenusSpecies, "unknown.")) %>%
  filter(!is.na(Location)) %>% 
  separate(SiteCode, c("big_plot", "corner", "small_plot"), sep="\\.")
     # Expect missing pieces for 100m2 plots
surveys_w_plots <- surveys_w_plot %>% 
  select(Location, Site, big_plot, corner, small_plot, GenusSpecies)

## Scale Summaries

ones <- surveys_w_plots %>%
  filter(small_plot == 1) %>%
  group_by(Location, Site, big_plot, corner) %>%
  summarize(ones = n())

tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(Location, Site, big_plot, corner) %>%
  distinct(GenusSpecies) %>%
  summarize(tens = n())

big_plots <- surveys_w_plots %>%  ## Big_plots = 100m
  group_by(Location, Site, big_plot) %>%
  distinct(GenusSpecies) %>%
  summarize(hundreds = n())

big_plots <- mutate(big_plots, site_stat = paste(Location, Site, sep=""))

sites <- surveys_w_plots %>%
  group_by(Location, Site) %>%
  distinct(GenusSpecies) %>%
  summarize(richness = n())

blocks <- surveys_w_plots %>%
  group_by(block) %>%
  distinct(GenusSpecies) %>%
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
  group_by(block, site, cluster) %>%
  summarize(avg_ones = round(mean(ones, na.rm=TRUE), 0), sd_ones = round(sd(ones), 2),
            avg_tens = round(mean(tens), 0), sd_tens = round(sd(tens), 2),
            avg_hund = round(mean(hundreds), 0), sd_hund = round(sd(hundreds), 2))

site_avg_plot <- record_counts %>%
  group_by(status, cluster) %>%
  summarize("1" = round(mean(ones, na.rm=TRUE), 0),
            "10" = round(mean(tens), 0),
            "100" = round(mean(hundreds), 0)) %>% 
  gather(scale, average, "1":"100")  %>% 
  mutate(cluster = factor(cluster, levels=c("open", "lawn", "orchard", 
                                            "hammock", "flatwood")))

ggplot(site_avg_plot, aes(x=scale, y=average, group=cluster)) +
  geom_line() +
  geom_point(size=3, alpha = 0.8,  aes(shape=status, color=as.factor(cluster))) +
  labs(x="Scale [m2]", y="Average Richness", shape="Soil Disturbance", color="Habitat") +
  theme_classic(base_size=14, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00"))
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

big_plots_cluster <- left_join(big_plots, cluster)
record_counts_cluster <- left_join(record_counts, cluster)

sink("output/richness_cluster.txt")

hundreds_aov <- aov(hundreds ~ cluster, data=big_plots_cluster)
cat("hundreds ~ cluster\n") 
print(summary(hundreds_aov))
HSD.test(hundreds_aov, "cluster")$groups
hundreds_aov <- aov(hundreds ~ status, data=big_plots_cluster)
cat("\nhundreds ~ status\n") 
print(summary(hundreds_aov))
HSD.test(hundreds_aov, "status")$groups

tens_aov <- aov(tens ~ cluster, data=record_counts_cluster)
cat("\ntens ~ cluster\n") 
print(summary(tens_aov))
HSD.test(tens_aov, "cluster")$groups
tens_aov <- aov(tens ~ status, data=record_counts_cluster)
cat("\ntens ~ status\n") 
print(summary(tens_aov))
HSD.test(tens_aov, "status")$groups

ones_aov <- aov(ones ~ cluster, data=record_counts_cluster)
cat("\nones ~ cluster\n") 
summary(ones_aov)
HSD.test(ones_aov, "cluster")$groups
ones_aov <- aov(ones ~ status, data=record_counts_cluster)
cat("\nones ~ status\n") 
summary(ones_aov)
HSD.test(ones_aov, "status")$groups

sink()

## Status Counts

block_status <- surveys_w_plots %>%
  left_join(cluster) %>% 
  group_by(cluster, block, site) %>%
  distinct(GenusSpecies) %>% 
  left_join(status) %>% 
  left_join(invasive) %>% 
  summarize(richness = n(), natives = sum(native, na.rm=TRUE),
            established = sum(established_FL, na.rm=TRUE) - sum(native, na.rm=TRUE),
            invasive = sum(invasive, na.rm=TRUE))

cluster_status <- block_status %>% 
  group_by(cluster) %>% 
  summarize(avg_richness = mean(richness),
            avg_natives = mean(natives),
            avg_established = mean(established),
            avg_invasive = mean(invasive))



## Species Area Curves [[TREC ONLY]]

richness_results <- record_counts %>% 
  ungroup() %>% 
  select(cluster, "1" = ones, "10" = tens, "100" = hundreds) %>%
  filter(!is.na(cluster)) %>% 
  gather(scale, richness, "1":"100")

sink("output/species_area.txt")
print(paste("habitat", "r2", "slope","intercept"))
for (habitat in unique(cluster$cluster)) {
  habitat_results <- filter(richness_results, cluster == habitat)
  test <- lm(log10(richness) ~ log10(as.numeric(scale)), data=habitat_results)
  print(paste(habitat, round(summary(test)$r.squared, 3), 
              round(summary(test)$coeff[2,1], 3),
              round(summary(test)$coeff[1,1], 3)))
}
sink()
