library(tidyverse)

### organizes and summarizes a list of species to count the number of
### species records per family
family_count <- function(species_list){
  families_table <- species_list %>%
    group_by(Family) %>%
    summarize(Count = n()) %>%
    arrange(desc(Count))
  return(families_table)
}

### plots a count of species per family in a histagram
plot_family_count <- function(families_table){
  fig <- ggplot(families_table, aes(x=Family, y=Count)) + 
    geom_bar(stat="identity", width = 0.9) +
    ylab("Species Count") +
    theme_classic(base_size=14, base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
}

### load data
list <- read_csv("data/TREClist.csv") %>%
  separate(Genus, c("Genus"), " \\[",  extra = "drop") %>%  # Removes synonyms
  separate(Species, c("Species"), " \\[",  extra = "drop")
families <- read_csv("data/TRECfamilies.csv")

### Add UF/IFAS Assessment
assessment <- read_csv("data/AssessmentList.csv") %>%
  separate(Title, c("Genus", "Species"), " ", extra = "merge")
assessment_alt <- assessment %>%
  filter(!is.na(Synonyms))

synonym_list <- data.frame(Name = c(), Assessment = c())
for (i in 1:nrow(assessment_alt)){
  synonyms <- str_split(assessment_alt$Synonyms[i], ", ")[[1]]
  for (synonym in synonyms){
    synonym_entry <- data.frame(Name = c(synonym), 
                                Assessment = c(assessment_alt$South[i]))
    synonym_list <- rbind(synonym_list, synonym_entry)
  }
}

synonym_list <- synonym_list %>%
  separate(Name, c("Genus", "Species"), " ", extra = "merge")

assessment_list <- assessment %>%
  select(Genus, Species, Assessment = South) %>%
  rbind(synonym_list) %>%
  distinct(Genus, Species, .keep_all=TRUE)

list_with_assessment <- left_join(list, assessment_list)
#write_csv(list_with_assessment, "data/TREClist.csv")

assessment_invasive <- list_with_assessment %>%
  filter(Assessment %in% c("Invasive", "High Invasion Risk", "Invasive (No Uses)", "Prohibited")) %>%
  select(Genus, Species, Assessment)
#write_csv(assessment_invasive, "data/TREClist-invasive.csv")

### All Families
list <- select(list, Family, Genus, Species, SpeciesAcro, MainCodeUSDA,
               CultivatedOnly, EstabdInFL, Native)

all_families <- family_count(list)
#write_csv(all_families, "output/TRECfamilies-all.csv")

top_families <- filter(all_families, !is.na(Family) & Count>1)
top_families <- transform(top_families, Family = reorder(Family, -Count))
plot_family_count(top_families)
#ggsave("output/TRECfamilies-all.png", width = 10, height = 4)

### Native Families
natives <- filter(list, Native == 1)
native_families <- family_count(natives)
#write_csv(native_families, "output/TRECfamilies-native.csv")

top_natives <- filter(native_families, !is.na(Family) & Count>1)
top_natives <- transform(top_natives, Family = reorder(Family, -Count))
plot_family_count(top_natives)
#ggsave("output/TRECfamilies-native.png", width = 5, height = 4)

### Non-native Established Families
establishedNonNative <- filter(list, Native == 0 & EstabdInFL == 1)
establishedNonNative_families <- family_count(establishedNonNative)
#write_csv(establishedNonNative_families, "output/TRECfamilies-establishedNonNative.csv")

top_establieshedNonNative <- filter(establishedNonNative_families, 
                                    !is.na(Family) & Count>1)
top_establieshedNonNative <- transform(top_establieshedNonNative, 
                                       Family = reorder(Family, -Count))
plot_family_count(top_establieshedNonNative)
#ggsave("output/TRECfamilies-establishedNonNative.png", width = 5, height = 4)

### Non-Established Families
nonEstablished <- filter(list, EstabdInFL == 0)
nonEstablished_families <- family_count(nonEstablished)
#write_csv(nonEstablished_families, "output/TRECfamilies-nonEstablished.csv")

top_nonEstablished <- filter(nonEstablished_families, !is.na(Family) & Count>1)
top_nonEstablished <- transform(top_nonEstablished, 
                                       Family = reorder(Family, -Count))
plot_family_count(top_nonEstablished)
#ggsave("output/TRECfamilies-nonEstablished.png", width = 5, height = 4)

### Unknown
unknown <- filter(list, Native == "?" | EstabdInFL == "?")
unknown_families <- family_count(unknown)

### Combined Table
all_families_counts <- all_families %>%
  left_join(native_families, by = "Family", suffix = c("", "_native")) %>%
  left_join(establishedNonNative_families, by = "Family",
            suffix = c("", "_estabNonNat")) %>%
  left_join(nonEstablished_families, by = "Family", 
            suffix = c("", "_nonEstab")) %>%
  left_join(unknown_families, by = "Family", 
            suffix = c("_all", "_unknown")) %>%
  filter(!is.na(Family)) %>%
  mutate_all(funs(replace(., is.na(.), 0))) %>%
  mutate(Count_sum = Count_native + Count_estabNonNat + Count_nonEstab + Count_unknown) %>%
  filter(Count_sum>1)

all_families_counts <- transform(all_families_counts, 
                                 Family = reorder(Family, -Count_sum))
all_families_counts_fig <- select(all_families_counts, -Count_all, -Count_sum)
all_families_counts_fig <- gather(all_families_counts_fig, key= "Class", value = "Count",
                              c(Count_native, Count_estabNonNat, Count_nonEstab, Count_unknown))
all_families_counts_fig$Class <- factor(all_families_counts_fig$Class,
                                        levels = c("Count_unknown", "Count_nonEstab", "Count_estabNonNat", "Count_native"))

ggplot(all_families_counts_fig, aes(x=Family, y=Count, fill=Class)) +
  geom_bar(stat="identity", width = 0.9) +
  labs(y = "Species Count", fill = "Local Establishment") +
  scale_y_continuous(breaks=seq(0,60,5)) +
  scale_fill_manual(values=c("#D3D3D3", "#00BFC4", "#F8766D", "#7CAE00"),
                    labels=c("Unknown", "Non-established", "Established Non-Natives", "Native")) +
  theme_classic(base_size=18, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        axis.title=element_text(size=24), legend.text = element_text(size=20),
        legend.title = element_text(size=24), legend.position = c(0.8, 0.8))
ggsave("output/TRECfamilies-allClass.png", width = 12, height = 10)
