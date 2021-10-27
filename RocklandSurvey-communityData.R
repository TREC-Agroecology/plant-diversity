### Create 'Community data' table for vegan package for 2020 NEON sampling data

library(vegan)
library(betapart)
library(tidyverse)


build_community_data <- function(species, plots, survey){
  community_data <- matrix(nrow = nrow(plots), ncol = nrow(species))
  colnames(community_data) <- species$GenusSpecies
  for (c in 1:nrow(species)){
    print(c/nrow(species))
    target_species <- species$GenusSpecies[c]
    for (r in 1:nrow(plots)){
      survey_at_scale <- suppressMessages(left_join(plots[r, ], survey))
      survey_target <- survey_at_scale %>% 
        filter(GenusSpecies == target_species)
      community_data[r, c] <- nrow(survey_target)
    }
  }
  return(community_data)
}

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

surveys_w_plots <- surveys %>%
  filter(!str_detect(GenusSpecies, "unknown.")) %>%
  filter(!is.na(Location)) %>% 
  separate(SiteCode, c("big_plot", "corner", "small_plot"), sep="\\.")

all_species <- surveys_w_plots %>%
  distinct(GenusSpecies) %>%
  arrange(GenusSpecies)

tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(Location, Site, big_plot, corner) %>%
  distinct(GenusSpecies)

hundreds <- tens %>%
  group_by(Location, Site, big_plot) %>%
  distinct(GenusSpecies)

plots_tens <- distinct(tens, Location, Site, big_plot, corner) %>% 
  mutate(Site_code = paste(Location, Site, sep = ""))
plots_hundreds <- distinct(hundreds, Location, Site, big_plot) %>%
  mutate(Site_code = paste(Location, Site, sep = ""))
plots_Site <- distinct(surveys_w_plots, Location, Site) %>%
  mutate(Site_code = paste(Location, Site, sep = ""))

location_colors <-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                    "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")
## Species List

all_species_count <- surveys_w_plots %>%
  group_by(GenusSpecies) %>%
  summarize(count=n())
write_csv(all_species_count, "output/rockland_survey_species.csv")

## Shannon diversity and evenness at various scales

matrix_ten <- build_community_data(all_species, plots_tens, tens)
diversity_ten <- diversity(matrix_ten)
evenness_ten <- diversity_ten/log(specnumber(matrix_ten))

matrix_hundred <- build_community_data(all_species, plots_hundreds, hundreds)
diversity_hundred <- diversity(matrix_hundred)
evenness_hundred <- diversity_hundred/log(specnumber(matrix_hundred))

#matrix_Site <- build_community_data(all_species, plots_Site, tens)
#diversity_Site <- diversity(matrix_Site)
#evenness_Site <- diversity_Site/log(specnumber(matrix_Site))

## Plot Diversity and Evenness
diversity_ten_tbl <- plots_tens %>%
  bind_cols(diversity = diversity_ten)
ggplot(diversity_ten_tbl, aes(x=diversity, fill=Site)) +
  geom_histogram(binwidth = 0.25) +
  facet_grid(Location~.) +
  theme_bw()

## Non-metric Multidimensional Analysis [PCOA / metaMDS with cmdscale(), vegdist()]

### Tens

dist_plots_10 <- vegdist(matrix_ten, "bray")
nmds_plots_10 <- metaMDS(dist_plots_10, k=3, try=200, trace=TRUE)

nmds_plots_scores_10 <- plots_tens %>%
  bind_cols(NMDS1 = scores(nmds_plots_10)[,1], NMDS2 = scores(nmds_plots_10)[,2])

ggplot(nmds_plots_scores_10, aes(x=NMDS1, y=NMDS2, shape=Site, color=Location)) +
                            #label = paste(big_plot,corner))) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Location") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values = location_colors)
ggsave("output/rockland-nmds-10.png", width = 10, height = 6)

perm_plots_10a <- adonis2(dist_plots_10 ~ plots_tens$Location + plots_tens$Site,
                         permutations = 1000)

perm_plots_10b <- adonis2(dist_plots_10 ~ plots_tens$Site + plots_tens$Location,
                         permutations = 1000)

perm_plots_10c <- adonis(dist_plots_10 ~ plots_tens$Location + plots_tens$Site +
                           plots_tens$Site_code, permutations = 1000)

perm_plots_10d <- adonis2(dist_plots_10 ~ plots_tens$Site + plots_tens$Location +
                        plots_tens$Site_code, permutations = 1000)

sink("output/rockland-permanova-10.txt")
print(perm_plots_10a)
print(perm_plots_10b)
print(perm_plots_10c)
print(perm_plots_10d)
sink()


### Hundreds

dist_plots_100 <- vegdist(matrix_hundred, "bray")
nmds_plots_100 <- metaMDS(dist_plots_100, k=2, try=100, trace=TRUE)

nmds_plots_scores_100 <- plots_hundreds %>%
  bind_cols(NMDS1 = scores(nmds_plots_100)[,1], NMDS2 = scores(nmds_plots_100)[,2])

ggplot(nmds_plots_scores_100, aes(x=NMDS1, y=NMDS2, shape=Site, color=Location)) +
  #label = big_plot)) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Location") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values = location_colors)
ggsave("output/rockland-nmds-100.png", width = 10, height = 6)

perm_plots_100a <- adonis2(dist_plots_100 ~ plots_hundreds$Location + plots_hundreds$Site,
                          permutations = 1000)

perm_plots_100b <- adonis2(dist_plots_100 ~ plots_hundreds$Site + plots_hundreds$Location,
                          permutations = 1000)

perm_plots_100c <- adonis(dist_plots_100 ~ plots_hundreds$Location + plots_hundreds$Site +
                           plots_hundreds$Site_code, permutations = 1000)

perm_plots_100d <- adonis2(dist_plots_100 ~ plots_hundreds$Site + plots_hundreds$Location +
                            plots_hundreds$Site_code, permutations = 1000)

sink("output/rockland-permanova-100.txt")
print(perm_plots_100a)
print(perm_plots_100b)
print(perm_plots_100c)
print(perm_plots_100d)
sink()