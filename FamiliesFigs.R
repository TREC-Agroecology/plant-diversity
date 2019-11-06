### This script develops the families histrograms for publication.

library(tidyverse)

TREC <- read_csv("data/TRECfamilyCounts.csv") %>%
  mutate(Total = `Native` + `Non-native, established` + `Non-native, not established` + `Uncertain`) %>%
  arrange(desc(Total)) %>%
  mutate(Family = factor(Family, levels = reorder(Family, Total)))

TREC_long <- gather(TREC, class, count, -Family) %>%
  filter(class != "Total") %>%
  mutate(class = factor(class, levels = c("Uncertain", "Non-native, not established", "Non-native, established", "Native")))

ggplot(TREC_long, aes(x=Family, y=count, fill=class)) +
  geom_bar(stat="identity")+
  labs(x="Family", y="Species Count", fill="Classification") +
  scale_fill_manual(values=c("#a6cee3", "#b2df8a", "#33a02c", "#1f78b4")) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("output/TRECfamiliesCount.png")


ECHO <- read_csv("data/ECHOfamilyCounts.csv") %>%
  mutate(Total = `Native` + `Non-native, established` + `Non-native, not established` + `Uncertain`) %>%
  arrange(desc(Total)) %>%
  mutate(Family = factor(Family, levels = reorder(Family, Total)))

ECHO_long <- gather(ECHO, class, count, -Family) %>%
  filter(class != "Total") %>%
  mutate(class = factor(class, levels = c("Uncertain", "Non-native, not established", "Non-native, established", "Native")))

ggplot(ECHO_long, aes(x=Family, y=count, fill=class)) +
  geom_bar(stat="identity")+
  labs(x="Family", y="Species Count", fill="Classification") +
  scale_fill_manual(values=c("#a6cee3", "#b2df8a", "#33a02c", "#1f78b4")) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("output/ECHOfamiliesCount.png")
