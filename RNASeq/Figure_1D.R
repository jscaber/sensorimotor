library(tidyverse)

c5ips <- read_tsv("coding/featurecounts_3/gsea/c5.all.v7.0.entrez.tsv")
c5 <- subset(c5ips, padj <0.05) %>% arrange(NES) 
c5 <- c5[grep("motor|sens|neur|nerv|spina", c5$pathway, ignore.case = TRUE),]
c5$pathway <- tolower(gsub("_"," ",c5$pathway))
c5$pathway <- gsub ("gobp", "GO:BP", c5$pathway)
c5$pathway <- gsub ("hp", "HP", c5$pathway)

CairoSVG("Figure2B.svg", width = 12, height = 4)
topc5 <- top_n(c5,4,wt = -NES)
bottomc5 <- top_n(c5,5,wt = NES)
full_join(topc5,bottomc5) %>% 
  mutate(pathway = fct_reorder(pathway, dplyr::desc(NES))) %>% 
  mutate(sign = if_else(NES >= 0, 'MN', 'SN')) %>%
  ggplot(aes(NES, pathway, fill=sign))+geom_bar(stat="identity") +
  scale_fill_manual(values = c("MN" = "red", "SN" = "blue")) + xlab("Normalised Enrichment Score (NES)")+
  guides(fill=guide_legend(title="Enriched in")) + theme_cowplot() + theme(axis.title.y = element_blank())
dev.off()