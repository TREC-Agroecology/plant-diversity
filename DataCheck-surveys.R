### This code cross references a record check by CM and TF.

library(tidyverse)

surveys <- read_csv("data/PlantDiversitySurvey-surveys.csv") %>%
  select(block, site, code, genus, species)

surveys_check <- read_csv("data/TREC-TFcheck.csv") %>%
  select(block, site, code, genus, species)
surveys_check$site <- str_replace(surveys_check$site, "North", "N")
surveys_check$site <- str_replace(surveys_check$site, "South", "S")
surveys_check$code <- str_replace(surveys_check$code, "31.100", "31")
surveys_check$code <- str_replace(surveys_check$code, "32.100", "32")
surveys_check$code <- str_replace(surveys_check$code, "40.100", "40")
surveys_check$code <- str_replace(surveys_check$code, "41.100", "41")

check <- anti_join(surveys_check, surveys)
check_back <- anti_join(surveys, surveys_check)

genus_check <- distinct(surveys, genus) %>%
  arrange(genus)

species_check <- distinct(surveys, genus, species) %>%
  arrange(genus, species)

write.csv(species_check, "output/PlantDiversitySampling-SpeciesList.csv")
