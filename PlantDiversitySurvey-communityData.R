### Create 'Community data' table for vegan package from NEON sampling data

library(vegan)
library(tidyverse)


build_community_data <- function(species, plots, survey){
  community_data <- matrix(nrow = nrow(plots), ncol = nrow(species))
  colnames(community_data) <- all_species$genus_species
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
  mutate(genus_species = paste(tolower(str_extract(taxonID, "....")),
                               tolower(str_extract(taxonIDRemarks, "..")), 
                               sep=".")) %>%
  separate(code, c("big_plot", "corner", "small_plot"), sep="\\.")

all_species <- surveys_w_plots %>%
  filter(!is.na(taxonIDRemarks)) %>%
  distinct(genus_species) %>%
  arrange(genus_species)
  
tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(genus_species)

plots_tens <- distinct(tens, big_plot, corner)
plots_hundreds <- distinct(surveys_w_plots, block, site, big_plot)
plots_site <- distinct(surveys_w_plots, block, site)


## Shannon diversity and evenness at various scales

matrix_ten <- build_community_data(all_species, plots_tens, tens)
diversity_ten <- diversity(matrix_ten)
evenness_ten <- diversity_ten/log(specnumber(matrix_ten))

matrix_hundred <- build_community_data(all_species, plots_hundreds, tens)
diversity_hundred <- diversity(matrix_hundred)
evenness_hundred <- diversity_hundred/log(specnumber(matrix_hundred))

matrix_site <- build_community_data(all_species, plots_site, tens)
diversity_site <- diversity(matrix_site)
evenness_site <- diversity_site/log(specnumber(matrix_site))


## Non-metric Multidimensional Analysis [PCOA / metaMDS with cmdscale(), vegdist()]

dist_plots <- vegdist(matrix_ten, "bray")
nmds_plots <- metaMDS(dist_plots, k=2, trace=TRUE)
#stressplot(nmds_plots)
#ordiplot(nmds_plots, display="sites", cex=1.25)

nmds_plots_scores <- plots_tens %>%
  bind_cols(NMDS1 = scores(nmds_plots)[,1], NMDS2 = scores(nmds_plots)[,2])

ggplot(nmds_plots_scores, aes(x=NMDS1, y=NMDS2, shape=site, color=as.factor(block))) +
  geom_point(cex=3) +
  theme_bw()


## Ordination tutorial
model <- decorana(matrix_site)
shnam <- make.cepnames(colnames(matrix_site))
pl <- plot(model, dis="sp")
identify(pl, "sp", labels=shnam)

stems <- colSums(matrix_site)
plot(model, dis="sp", type="n")
sel <- orditorp(model, dis="sp", lab=shnam, priority=stems, pcol = "gray", pch="+")

plot(model, dis="sp", type="n")
ordilabel(model, dis="sp", lab=shnam, priority = stems)
