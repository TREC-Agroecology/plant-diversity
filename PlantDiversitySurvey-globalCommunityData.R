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
  filter(!is.na(genus)) %>%
  mutate(genus_species = paste(tolower(str_extract(genus, "...")),
                               tolower(str_trunc(species, 3, "right", "")), 
                               sep="")) %>%
  separate(code, c("big_plot", "corner", "small_plot"), sep="\\.") %>% 
  # Expect missing pieces for 100m2 plots
  select(block, site, big_plot, corner, small_plot, genus_species)

echo <- read.csv("data/ECHO-surveys.csv") %>% 
  filter(UncertainId %in% c(NA, "Species")) %>%
  mutate(genus_species = tolower(genus_species)) %>%
  separate(code, c("big_plot", "corner", "small_plot"), sep="\\.") %>%
  select(block, site, big_plot, corner, small_plot, genus_species)

surveys_w_plots <- bind_rows(surveys_w_plots, echo)

all_species <- surveys_w_plots %>%
  distinct(genus_species) %>%
  arrange(genus_species)
  
tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(genus_species)

hundreds <- surveys_w_plots %>%
  group_by(block, site, big_plot) %>%
  distinct(genus_species)

status <- data.frame(block = c(1, 4, 14, 15, 31, 32), 
                     status = c("old", "old", "new", "new", "old", "new"))
cluster <- data.frame(block = c(1, 1, 4, 4, 14, 14, 15, 15, 31, 31, 32, 32),
                      site = rep(c("N", "S"), 6),
                      cluster = c("open", "open", "open", "lawn",
                                  "lawn", "hammock", "orchard", "orchard",
                                  "echo", "echo", "echo", "echo"))
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
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))
ggsave("output/nmds-10-global.png", width = 8, height = 6)

perm_plots_10 <- adonis(dist_plots_10 ~ plots_tens$status + plots_tens$cluster +
                     plots_tens$block + plots_tens$site_code, permutations = 1000)

sink("output/permanova-10-global.txt")
print(perm_plots_10)
sink()


### Hundreds

dist_plots_100 <- vegdist(matrix_hundred, "bray")
nmds_plots_100 <- metaMDS(dist_plots_100, k=2, try=100, trace=TRUE)

nmds_plots_scores_100 <- plots_hundreds %>%
  bind_cols(NMDS1 = scores(nmds_plots_100)[,1], NMDS2 = scores(nmds_plots_100)[,2])

ggplot(nmds_plots_scores_100, aes(x=NMDS1, y=(-1*NMDS2), shape=site, color=as.factor(block))) +
                              #label = big_plot)) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Block") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"))
ggsave("output/nmds-100-global.png", width = 8, height = 6)

perm_plots_100 <- adonis(dist_plots_100 ~ plots_hundreds$status + plots_hundreds$cluster +
                       plots_hundreds$block + plots_hundreds$site_code, permutations = 1000)

sink("output/permanova-100-global.txt")
print(perm_plots_100)
sink()

### Hundreds Rank
rm(species_rank)
names_list <- c("genus_species")
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
  select(genus_species, starts_with("rank_"))

sink("output/top-species-global.txt")

top_open <- rank_table %>%
  filter(rank_open <= 3) %>%
  arrange(rank_open) %>%
  select(genus_species, everything(), -rank_open)

cat("Open\n")
as.data.frame(top_open)

top_lawn <- rank_table %>%
  filter(rank_lawn == 1) %>%
  select(genus_species, everything(), -rank_lawn)

cat("\nLawn\n")
as.data.frame(top_lawn)

top_orchard <- rank_table %>%
  filter(rank_orchard == 1) %>%
  select(genus_species, everything(), -rank_orchard)

cat("\nOrchard\n")
as.data.frame(top_orchard)

top_hammock <- rank_table %>%
  filter(rank_hammock == 1) %>%
  select(genus_species, everything(), -rank_hammock)

cat("\nHammock\n")
as.data.frame(top_hammock)

top_echo <- rank_table %>%
  filter(rank_echo <= 3) %>%
  select(genus_species, everything(), -rank_echo)

cat("\nECHO\n")
as.data.frame(top_echo)
sink()
