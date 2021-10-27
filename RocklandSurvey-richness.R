### This script organizes and evaluates the NEON Plant Diversity Survey
### conducted at pine rocklands in Miami-Dade Summer 2021.
### Working directory is 'Source File Location'.


library(tidyverse)
library(agricolae)

## Data

data_files <- str_subset(list.files("data/"), "_NEON")
for (file in data_files) {
  Site_data <- read_csv(paste0("data/", file))
  if (exists("surveys")){
    surveys <- bind_rows(surveys, Site_data)
  } else {
    surveys <- Site_data
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
location_colors <-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                    "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")

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

big_plots <- mutate(big_plots, Site_stat = paste(Location, Site, sep=""))

Sites <- surveys_w_plots %>%
  group_by(Location, Site) %>%
  distinct(GenusSpecies) %>%
  summarize(richness = n())

Locations <- surveys_w_plots %>%
  group_by(Location) %>%
  distinct(GenusSpecies) %>%
  summarize(richness = n())


## Combined Summaries

record_counts <- ones %>% 
  right_join(tens) %>% 
  inner_join(big_plots) %>% 
  mutate(Site_stat = paste(Location, Site, sep=""))

record_counts_Sites <- ones %>% 
  right_join(tens) %>% 
  inner_join(big_plots) %>% 
  mutate(Site_stat = paste(Location, Site, sep=""))


## Calculate and visualize richness among scales

length(unique(surveys_w_plots$GenusSpecies))

ggplot(tens, aes(x=tens, fill=Site)) +  # tens visualization
  geom_histogram(binwidth = 5) +
  facet_wrap(vars(Location)) +
  labs(x="Richness", y="Count", fill="Site") +
  theme_bw(base_size=24, base_family="Helvetica")
ggsave("output/rockland_richness_ten.png", width = 12, height = 5)  

Site_avg <- record_counts %>%
  group_by(Location, Site) %>%
  summarize(avg_ones = round(mean(ones, na.rm=TRUE), 0), sd_ones = round(sd(ones), 2),
            avg_tens = round(mean(tens), 0), sd_tens = round(sd(tens), 2),
            avg_hund = round(mean(hundreds), 0), sd_hund = round(sd(hundreds), 2))

Site_avg_plot <- record_counts %>%
  group_by(Site, Location) %>%
  summarize("1" = round(mean(ones, na.rm=TRUE), 0),
            "10" = round(mean(tens), 0),
            "100" = round(mean(hundreds), 0)) %>% 
  gather(scale, average, "1":"100")

ggplot(Site_avg_plot, aes(x=scale, y=average)) +
  geom_jitter(size=3, alpha = 0.8,  aes(shape=Site, color=Location)) +
  labs(x="Scale [m2]", y="Average Richness", shape="Site", color="Location") +
  theme_classic(base_size=14, base_family="Helvetica") +
  scale_colour_manual(values=location_colors)
ggsave("output/rockland_avg_richness.png", width = 7.5, height = 6)  

## ANOVA and Tukey Post Hoc

sink("output/rockland_richness.txt")

hundreds_aov <- aov(hundreds ~ Location + Site + Location*Site, data=big_plots)
cat("hundreds ~ Location + Site + L*S\n") 
summary(hundreds_aov)
cat("\nLocation\n")
HSD.test(hundreds_aov, "Location")$groups
cat("\nSite\n")
HSD.test(hundreds_aov, "Site")$groups

tens_aov <- aov(tens ~ Location + Site + Location*Site, data=record_counts)
cat("\ntens ~ Location + Site + L*S\n") 
summary(tens_aov)
cat("\nLocation\n")
HSD.test(tens_aov, "Location")$groups
cat("\nSite\n")
HSD.test(tens_aov, "Site")$groups

ones_aov <- aov(ones ~ Location + Site + Location*Site, data=record_counts)
cat("\nones ~ Location + Site + L*S\n") 
summary(ones_aov)
cat("\nLocation\n")
HSD.test(ones_aov, "Location")$groups
cat("\nSite\n")
HSD.test(ones_aov, "Site")$groups

sink()


## Status Counts

Location_status <- surveys_w_plots %>%
  group_by(Location, Site) %>%
  mutate(GenusSpecies = str_remove(GenusSpecies, ".\\.")) %>% 
  distinct(GenusSpecies) %>% 
  left_join(status) %>% 
  left_join(invasive) %>% 
  summarize(richness = n(), natives = sum(native, na.rm=TRUE),
            established = sum(established_FL, na.rm=TRUE) - sum(native, na.rm=TRUE),
            invasive = sum(invasive, na.rm=TRUE))



## Species Area Curves

sink("output/rockland_species_area.txt")

print(paste("--location", "r2", "slope","intercept"))
for (location in unique(Site_avg_plot$Location)) {
  location_results <- filter(Site_avg_plot, Location == location)
  test <- lm(log10(average) ~ log10(as.numeric(scale)), data=location_results)
  print(paste(location, round(summary(test)$r.squared, 3), 
              round(summary(test)$coeff[2,1], 3),
              round(summary(test)$coeff[1,1], 3)))
}

print(paste("--site", "r2", "slope","intercept"))
for (site in unique(Site_avg_plot$Site)) {
  site_results <- filter(Site_avg_plot, Site == site)
  test <- lm(log10(average) ~ log10(as.numeric(scale)), data=site_results)
  print(paste(site, round(summary(test)$r.squared, 3), 
              round(summary(test)$coeff[2,1], 3),
              round(summary(test)$coeff[1,1], 3)))
}

sink()
