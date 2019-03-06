### Create 'Community data' table for vegan package from NEON sampling data

library(vegan)
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
  distinct(genus_species) %>%
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

status <- data.frame(block = c(1,4,14,15), status = c("old", "old", "new", "new"))
cluster <- data.frame(block = c(1, 1, 4, 4, 14, 14, 15, 15),
                      site = rep(c("N", "S"), 4),
                      cluster = c("open", "open", "open", "lawn",
                                  "lawn", "hammock", "orchard", "orchard"))
plots_tens <- distinct(tens, big_plot, corner) %>%
  mutate(site_code = paste(block, site, sep = "")) %>%
  left_join(status, by = c("block")) %>%
  left_join(cluster, by = c("block", "site"))
plots_hundreds <- distinct(surveys_w_plots, block, site, big_plot) %>%
  mutate(site_code = paste(block, site, sep = "")) %>%
  left_join(status, by = c("block")) %>%
  left_join(cluster, by = c("block", "site"))
plots_site <- distinct(surveys_w_plots, block, site)


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

#matrix_ten <- read.csv("data/matrix-ten.csv", row =1, header =T) 
dist_plots <- vegdist(matrix_ten, "bray")
#dist_plots_j <- vegdist(matrix_ten, "jaccard")
nmds_plots <- metaMDS(dist_plots, k=2, try=100, trace=TRUE)
#stressplot(nmds_plots)
#ordiplot(nmds_plots, display="sites", cex=1.25)

nmds_plots_scores <- plots_tens %>%
  bind_cols(NMDS1 = scores(nmds_plots)[,1], NMDS2 = scores(nmds_plots)[,2])

ggplot(nmds_plots_scores, aes(x=NMDS1, y=NMDS2, shape=site, color=as.factor(block),
                            label = paste(big_plot,corner))) +
  geom_point(cex=5) +
  geom_text(color="black")+
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Block") +
  theme_bw(base_size=20, base_family="Helvetica")
ggsave("output/nmds.png", width = 8, height = 6)

perm_plots <- adonis(dist_plots ~ plots_tens$status + plots_tens$block +
                     plots_tens$site_code, permutations = 1000)
hist(perm_plots$f.perms)
anosim_plots <- anosim(dist_plots, plots_tens$block, permutations = 1000)

sink("output/permanova.txt")
print(perm_plots)
sink()

### Tens Subset

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

### Hundreds

dist_plots <- vegdist(matrix_hundred, "bray")
nmds_plots <- metaMDS(dist_plots, k=2, try=100, trace=TRUE)

nmds_plots_scores <- plots_hundreds %>%
  bind_cols(NMDS1 = scores(nmds_plots)[,1], NMDS2 = scores(nmds_plots)[,2])

ggplot(nmds_plots_scores, aes(x=NMDS1, y=NMDS2, shape=site, color=as.factor(block),
                              label = big_plot)) +
  geom_point(cex=5) +
  geom_text(color="black")+
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Block") +
  theme_bw(base_size=20, base_family="Helvetica")
ggsave("output/nmds-100.png", width = 8, height = 6)

perm_plots <- adonis(dist_plots ~ plots_hundreds$status + plots_hundreds$cluster +
                       plots_hundreds$block + plots_hundreds$site_code, permutations = 1000)


sink("output/permanova-100.txt")
print(perm_plots)
sink()

rm(species_rank)
for (c in unique(plots_hundreds$cluster)){
  blocks <- filter(plots_hundreds, cluster == c)
  cluster_matrix <- build_community_data(all_species, blocks, hundreds)
  matrix_sum <- sort(colSums(cluster_matrix), decreasing = TRUE)
  if (exists("species_rank")){
    species_rank <- full_join(species_rank,
                              data.frame(genus_species = names(matrix_sum[matrix_sum>0]),
                                         rank = matrix_sum[matrix_sum>0]),
                              by = "genus_species")
  } else {
    species_rank <- left_join(all_species, 
                              data.frame(genus_species = names(matrix_sum[matrix_sum>0]),
                                         rank = matrix_sum[matrix_sum>0]),
                              by = "genus_species")
  }
}
names(species_rank) <- c("species_code", 
                         as.character(unique(clusters_100$cluster)))