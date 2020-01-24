### Create 'Community data' table for vegan package from NEON sampling data

library(vegan)
library(betapart)
library(tidyverse)


build_community_data <- function(species, plots, survey){
  community_data <- matrix(nrow = nrow(plots), ncol = nrow(species))
  colnames(community_data) <- species$genus_species
  for (c in 1:nrow(species)){
    target_species <- species$genus_species[c]
    for (r in 1:nrow(plots)){
      survey_at_scale <- suppressMessages(left_join(plots[r, ], survey))
      survey_target <- survey_at_scale %>% 
        filter(genus_species == target_species)
      community_data[r, c] <- nrow(survey_target)
    }
  }
  return(community_data)
}

## Data 

surveys <- read_csv("data/PlantDiversitySurvey-surveys.csv")
surveys_w_plots <- surveys %>%
  mutate(genus_species = paste(tolower(str_extract(genus, "....")),
                               tolower(str_extract(species, "...")), 
                               sep=".")) %>%
  separate(code, c("big_plot", "corner", "small_plot"), sep="\\.")

all_species <- surveys_w_plots %>%
  filter(!is.na(genus)) %>%
  distinct(genus_species, genus, species) %>%
  arrange(genus_species)
  
tens <- surveys_w_plots %>%
  filter(!is.na(genus)) %>%
  filter(!is.na(corner)) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(genus_species)

hundreds <- surveys_w_plots %>%
  filter(!is.na(genus)) %>%
  group_by(block, site, big_plot) %>%
  distinct(genus_species)

cluster <- data.frame(block = c(1, 1, 4, 4, 14, 14, 15, 15),
                      site = rep(c("N", "S"), 4),
                      status = c("high", "high", "low", "high",
                                 "high", "low", "low", "low"),
                      cluster = c("open", "open", "lawn", "open",
                                  "open", "hammock", "orchard", "orchard"))
plots_tens <- distinct(tens, big_plot, corner) %>%
  mutate(site_code = paste(block, site, sep = "")) %>%
  left_join(cluster, by = c("block", "site"))
plots_hundreds <- distinct(surveys_w_plots, block, site, big_plot) %>%
  mutate(site_code = paste(block, site, sep = "")) %>%
  left_join(cluster, by = c("block", "site"))
plots_site <- distinct(surveys_w_plots, block, site)

plots_tens_mixed <- read_csv("data/plots_tens_mixed.csv") # mixed habitat classificiation

## Shannon diversity and evenness at various scales

matrix_ten <- build_community_data(all_species, plots_tens, tens)
diversity_ten <- diversity(matrix_ten)
evenness_ten <- diversity_ten/log(specnumber(matrix_ten))

matrix_hundred <- build_community_data(all_species, plots_hundreds, hundreds)
diversity_hundred <- diversity(matrix_hundred)
evenness_hundred <- diversity_hundred/log(specnumber(matrix_hundred))

matrix_site <- build_community_data(all_species, plots_site, tens)
diversity_site <- diversity(matrix_site)
evenness_site <- diversity_site/log(specnumber(matrix_site))

## Plot Diversity and Evenness
diversity_ten_tbl <- plots_tens %>%
  bind_cols(diversity = diversity_ten)
ggplot(diversity_ten_tbl, aes(x=diversity, fill=site)) +
  geom_histogram(binwidth = 0.25) +
  facet_grid(block~.) +
  theme_bw()

## Non-metric Multidimensional Analysis [PCOA / metaMDS with cmdscale(), vegdist()]

### Tens

dist_plots_10 <- vegdist(matrix_ten, "bray")
nmds_plots_10 <- metaMDS(dist_plots_10, k=2, try=100, trace=TRUE)

nmds_plots_scores_10 <- plots_tens %>%
  bind_cols(NMDS1 = scores(nmds_plots_10)[,1], NMDS2 = scores(nmds_plots_10)[,2])

ggplot(nmds_plots_scores_10, aes(x=NMDS1, y=NMDS2, shape=site, color=as.factor(block))) +
                            #label = paste(big_plot,corner))) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Block") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))
ggsave("output/nmds-10.png", width = 8, height = 6)

perm_plots_10 <- adonis(dist_plots_10 ~ plots_tens$status + plots_tens$cluster +
                     plots_tens$block + plots_tens$site_code, permutations = 1000)

sink("output/permanova-10.txt")
print(perm_plots_10)
sink()


### Hundreds

dist_plots_100 <- vegdist(matrix_hundred, "bray")
nmds_plots_100 <- metaMDS(dist_plots_100, k=2, try=100, trace=TRUE)

nmds_plots_scores_100 <- plots_hundreds %>%
  bind_cols(NMDS1 = scores(nmds_plots_100)[,1], NMDS2 = scores(nmds_plots_100)[,2])

ggplot(nmds_plots_scores_100, aes(x=NMDS1, y=NMDS2, shape=site, color=as.factor(block))) +
  #label = big_plot)) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Block") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))
ggsave("output/nmds-100.png", width = 8, height = 6)

perm_plots_100 <- adonis(dist_plots_100 ~ plots_hundreds$status + plots_hundreds$cluster, permutations = 1000)

sink("output/permanova-100.txt")
print(perm_plots_100)
sink()

### Tens Habitat

dist_plots_10 <- vegdist(matrix_ten, "bray")
nmds_plots_10 <- metaMDS(dist_plots_10, k=2, try=100, trace=TRUE)

nmds_plots_scores_10 <- plots_tens_mixed %>%
  bind_cols(NMDS1 = scores(nmds_plots_10)[,1], NMDS2 = scores(nmds_plots_10)[,2])

ggplot(nmds_plots_scores_10, aes(x=NMDS1, y=NMDS2, shape=status, color=as.factor(cluster))) +
  #label = paste(big_plot,corner))) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Soil Disturbance", color="Habitat") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))
ggsave("output/nmds-10-habitat.png", width = 8, height = 6)

perm_plots_10_soil_habitat <- adonis2(dist_plots_10 ~ plots_tens_mixed$status + plots_tens_mixed$cluster,
                                      permutations = 1000)
perm_plots_10 <- adonis2(dist_plots_10 ~ plots_tens_mixed$status + plots_tens_mixed$cluster +
                          plots_tens_mixed$block + plots_tens_mixed$site_code, permutations = 1000)

sink("output/permanova-10-habitat.txt")
print(perm_plots_10_soil_habitat)
print(perm_plots_10)
sink()

### Hundreds Habitat

dist_plots_100 <- vegdist(matrix_hundred, "bray")
nmds_plots_100 <- metaMDS(dist_plots_100, k=2, try=100, trace=TRUE)

nmds_plots_scores_100 <- plots_hundreds %>%
  bind_cols(NMDS1 = scores(nmds_plots_100)[,1], NMDS2 = scores(nmds_plots_100)[,2])

ggplot(nmds_plots_scores_100, aes(x=NMDS1, y=NMDS2, shape=status, color=as.factor(cluster))) +
  #label = big_plot)) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Soil Disturbance", color="Habitat") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))
ggsave("output/nmds-100-habitat.png", width = 8, height = 6)

perm_plots_100_soil_habitat <- adonis2(dist_plots_100 ~ plots_hundreds$status + plots_hundreds$cluster,
                                       permutations = 1000)
perm_plots_100 <- adonis2(dist_plots_100 ~ plots_hundreds$status + plots_hundreds$cluster +
                            plots_hundreds$block + plots_hundreds$site_code, permutations = 1000)

sink("output/permanova-100-habitat.txt")
print(perm_plots_100_soil_habitat)
print(perm_plots_100)
sink()

### Hundreds Rank
rm(species_rank)
names_list <- c("genus", "species", "genus_species")
for (c in unique(plots_hundreds$cluster)){
  names_list <- c(names_list, c(paste("count", c, sep="_"), paste("rank", c, sep="_")))
  blocks <- filter(plots_hundreds, cluster == c)
  cluster_matrix <- build_community_data(all_species, blocks, hundreds)
  matrix_sum <- sort(colSums(cluster_matrix), decreasing = TRUE)
  if (exists("species_rank")){
    species_rank <- full_join(species_rank,
                              data.frame(genus_species = names(matrix_sum[matrix_sum>0]),
                                         count = matrix_sum[matrix_sum>0],
                                         rank = rank(-matrix_sum[matrix_sum>0], ties.method="min")),
                              by = "genus_species")
  } else {
    species_rank <- left_join(all_species, 
                              data.frame(genus_species = names(matrix_sum[matrix_sum>0]),
                                         count = matrix_sum[matrix_sum>0],
                                         rank = rank(-matrix_sum[matrix_sum>0], ties.method="min")),
                              by = "genus_species")
  }
}
names(species_rank) <- c(names_list)

rank_table <- species_rank %>%
  select(genus, species, starts_with("rank_"))

sink("output/top-species.txt")

top_open <- rank_table %>%
  filter(rank_open <= 3) %>%
  arrange(rank_open) %>%
  left_join(all_species) %>%
  select(genus, species, everything(), -rank_open, -genus_species)

cat("Open\n")
as.data.frame(top_open)

top_lawn <- rank_table %>%
  filter(rank_lawn == 1) %>%
  left_join(all_species) %>%
  select(genus, species, everything(), -rank_lawn, -genus_species)

cat("\nLawn\n")
as.data.frame(top_lawn)

top_orchard <- rank_table %>%
  filter(rank_orchard == 1) %>%
  left_join(all_species) %>%
  select(genus, species, everything(), -rank_orchard, -genus_species)

cat("\nOrchard\n")
as.data.frame(top_orchard)

top_hammock <- rank_table %>%
  filter(rank_hammock == 1) %>%
  left_join(all_species) %>%
  select(genus, species, everything(), -rank_hammock, -genus_species)

cat("\nHammock\n")
as.data.frame(top_hammock)
sink()

### EXTRA Tens Subset

species_rank <- tens %>%
  group_by(genus_species) %>%
  summarize(count=n()) %>%
  arrange(desc(count))

top_species <- filter(species_rank, count >5)

matrix_ten <- build_community_data(top_species, plots_tens, tens)
dist_plots <- vegdist(matrix_ten, "bray")
nmds_plots <- metaMDS(dist_plots, k=2, try=100, trace=TRUE)

nmds_plots_scores <- plots_tens %>%
  bind_cols(NMDS1 = scores(nmds_plots)[,1], NMDS2 = scores(nmds_plots)[,2])

ggplot(nmds_plots_scores, aes(x=NMDS1, y=NMDS2, shape=site, color=as.factor(block),
                              label = paste(big_plot,corner))) +
  geom_point(cex=5) +
  geom_text(color="black")+
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Block") +
  theme_bw(base_size=20, base_family="Helvetica")

perm_plots <- adonis(dist_plots ~ plots_tens$status + plots_tens$block +
                       plots_tens$site_code, permutations = 1000)
hist(perm_plots$f.perms)
anosim_plots <- anosim(dist_plots, plots_tens$block, permutations = 1000)

sink("output/permanova-10s.txt")
print(perm_plots)
sink()
  
