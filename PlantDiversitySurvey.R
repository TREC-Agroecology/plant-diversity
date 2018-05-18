### This script organizes and evaluates the NEON Plant Diversity Survey
### conducted at TREC in Homestead, FL Fall 2017 - Spring 2018

library(tidyverse)

surveys <- read_csv("data/PlantDiversitySurvey-surveys.csv")

blk1N <- filter(surveys, block == 1 & site == "N")

blk1N_1 <- blk1N[str_detect(blk1N$code, "\\.1$"),]

blk1N_10_code <- unique(blk1N$code[str_detect(blk1N$code, "\\.10")])

for (code in blk1N_10_code){
  corner <- str_replace(code, "\\.10", "")
  print(corner)
}