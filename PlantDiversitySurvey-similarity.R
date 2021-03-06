### This script generates jaccard dissimilarity scores at 10 m2 and 1 m2.
### Working directory is 'Source File Location'.

### TO DO:
### - Add % cover to ones
### - Redefine out group to exclude site neighbor



library(tidyverse)

## Jaccard Calculation

jaccard <- function(survey_data){
  
  plots <- distinct(tens, big_plot, corner)
  
  output_table <- data.frame()
  
  for (j in 1:nrow(plots)){
    
    ## Identify the target site records
    target <- left_join(plots[j, ], survey_data, 
                        by = c("block", "site", "big_plot", "corner"))
    target_records <- target %>%
      ungroup() %>%
      select(taxonID, taxonIDRemarks)
    
    ## Identify the plots that share a site with the target and calculate J
    in_surveys <- filter(survey_data, block == plots[j, ]$block, site == plots[j, ]$site)
    in_surveys <- anti_join(in_surveys, target,
                            by = c("block", "site", "big_plot", "corner", "taxonID", "taxonIDRemarks"))
    in_surveys_plots <- distinct(in_surveys, big_plot, corner)
    
    in_jaccard <- c()
    
    for (i in 1:nrow(in_surveys_plots)){
      other <- left_join(in_surveys_plots[i, ], survey_data,
                         by = c("block", "site", "big_plot", "corner"))
      other_records <- other %>%
        ungroup() %>%
        select(taxonID, taxonIDRemarks)
      and <- nrow(inner_join(target_records, other_records,
                             by = c("taxonID", "taxonIDRemarks")))
      or <- nrow(distinct(bind_rows(target_records, other_records)))
      jaccard <- and/or
      in_jaccard <- c(in_jaccard, jaccard)
    }
    
    ## Identify the plots at a different site than the target and calculate J
    out_surveys <- filter(survey_data, block != plots[j, ]$block)
    out_surveys <- anti_join(out_surveys, target,
                             by = c("block", "site", "big_plot", "corner", "taxonID", "taxonIDRemarks"))
    out_surveys_plots <- distinct(out_surveys, big_plot, corner)
    
    out_jaccard <- c()
    
    for (i in 1:nrow(out_surveys_plots)){
      other <- left_join(in_surveys_plots[i, ], survey_data,
                         by = c("block", "site", "big_plot", "corner"))
      other_records <- other %>%
        ungroup() %>%
        select(taxonID, taxonIDRemarks)
      and <- nrow(inner_join(target_records, other_records,
                             by = c("taxonID", "taxonIDRemarks")))
      or <- nrow(distinct(bind_rows(target_records, other_records)))
      jaccard <- and/or
      out_jaccard <- c(out_jaccard, jaccard)
    }
    
    output_row <- cbind(plots[j, ], avg_j_in = mean(in_jaccard), 
                        avg_j_out = mean(out_jaccard), 
                        j_diff = mean(in_jaccard) - mean(out_jaccard))
    output_table <- bind_rows(output_table, output_row)
  }
  
  output_table <- mutate(output_table, site_stat = paste(block, site, sep=""))
  
  jaccard_table <- output_table %>%
    group_by(block, site) %>%
    summarize(avg_jaccard_in = mean(avg_j_in),
              avg_jaccard_out = mean(avg_j_out),
              avg_diff = mean(j_diff))
  
  return(jaccard_table)
}

## Read data

surveys <- read_csv("data/PlantDiversitySurvey-surveys.csv")
surveys_w_plots <- surveys %>%
  separate(code, c("big_plot", "corner", "small_plot"), sep="\\.")


## Sort data into grid size

tens <- surveys_w_plots %>%
  filter(!is.na(corner)) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(taxonID, taxonIDRemarks)

ones <- surveys_w_plots %>%
  filter(small_plot == 1) %>%
  group_by(block, site, big_plot, corner) %>%
  distinct(taxonID, taxonIDRemarks)

## Variation across not within but between PERMANOVA (vegan:adonis)
## Plot distribution of Jaccard scores (by plot and block)
## betapart

## Jaccard Calculations

jaccard_tens <- jaccard(tens)
jaccard_ones <- jaccard(ones)


## ANOVA

avg_j_in <- aov(avg_j_in ~ block + site_stat, data=output_table)
summary(avg_j_in)

avg_j_out <- aov(avg_j_out ~ block + site_stat, data=output_table)
summary(avg_j_out)

j_diff <- aov(j_diff ~ block + site_stat, data=output_table)
summary(j_diff)

## Tukey Post Hoc

HSD.test(avg_j_in, "block")$groups
HSD.test(avg_j_out, "block")$groups
HSD.test(j_diff, "block")$groups
HSD.test(avg_j_in, "site_stat")$groups
HSD.test(avg_j_out, "site_stat")$groups
HSD.test(j_diff, "site_stat")$groups
