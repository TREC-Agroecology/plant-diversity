### Create 'Community data' table for vegan package for 2020 NEON sampling data

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

surveys <- read_csv("data/NEONsurveyMar2020PlotData.csv")
status <- read_csv("data/TRECstatus.csv")
invasive <- read_csv("data/TRECinvasive.csv")

surveys_w_plots <- surveys %>%
  filter(Certainty == "ToSpp") %>%
  separate(Subplot, c("big_plot", "corner", "small_plot"), sep="\\.") %>% 
  select(block = Block, site = Land, big_plot, corner, small_plot, 
         genus_species = Code_CGM)

surveys_w_plots$block <- str_replace(surveys_w_plots$block, "B02", "B01")

pub_blocks <- data.frame(block = c("B01", "B04", "B14", "B15"), 
                         pub_block = factor(c("Old-Ag-NW", "Old-Ag-NE", "New-Ag-SW", "New-Ag-SE"), 
                                            levels = c("Old-Ag-NW", "Old-Ag-NE", "New-Ag-SW", "New-Ag-SE")))

pub_sites <- data.frame(site = c("cc", "gr"), pub_site = c("Ag", "Lawn"))

all_species <- surveys_w_plots %>%
  distinct(genus_species) %>%  # rm(, genus, species)
  arrange(genus_species)  # %>% 
#  left_join(status) %>% 
#  left_join(invasive)

tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(genus_species)

hundreds <- tens %>%
  group_by(block, site, big_plot) %>%
  distinct(genus_species)

plots_tens <- distinct(tens, block, site, big_plot, corner) %>% 
  mutate(site_code = paste(block, site, sep = ""))
plots_hundreds <- distinct(hundreds, block, site, big_plot) %>%
  mutate(site_code = paste(block, site, sep = ""))
plots_site <- distinct(surveys_w_plots, block, site) %>%
  mutate(site_code = paste(block, site, sep = ""))

## Species List

all_species_count <- surveys_w_plots %>%
  group_by(genus_species) %>%
  summarize(count=n())
write_csv(all_species_count, "output/2020-TREC-survey-list.csv")

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
  bind_cols(NMDS1 = scores(nmds_plots_10)[,1], NMDS2 = scores(nmds_plots_10)[,2]) %>% 
  left_join(pub_blocks) %>% 
  left_join(pub_sites)

ggplot(nmds_plots_scores_10, aes(x=NMDS1, y=NMDS2, shape=pub_site, color=pub_block)) +
                            #label = paste(big_plot,corner))) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Block") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))
ggsave("output/2020-nmds-10.png", width = 8, height = 6)

perm_plots_10 <- adonis(dist_plots_10 ~ plots_tens$site + plots_tens$block +
                        plots_tens$site_code, permutations = 1000)

sink("output/permanova-10.txt")
print(perm_plots_10)
sink()


### Hundreds

dist_plots_100 <- vegdist(matrix_hundred, "bray")
nmds_plots_100 <- metaMDS(dist_plots_100, k=2, try=100, trace=TRUE)

nmds_plots_scores_100 <- plots_hundreds %>%
  bind_cols(NMDS1 = scores(nmds_plots_100)[,1], NMDS2 = scores(nmds_plots_100)[,2])

ggplot(nmds_plots_scores_100, aes(x=NMDS1, y=NMDS2, shape=site, color=block)) +
  #label = big_plot)) +
  geom_point(cex=5) +
  #geom_text(color="black") +
  labs(x="NMDS1", y="NMDS2", shape="Site", color="Block") +
  theme_bw(base_size=20, base_family="Helvetica") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2"))
ggsave("output/2020-nmds-100.png", width = 8, height = 6)

perm_plots_100 <- adonis(dist_plots_100 ~ plots_hundreds$site + plots_hundreds$block +
                           plots_hundreds$site_code, permutations = 1000)

sink("output/permanova-100.txt")
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