library(tidyverse)
library(DESeq2)

c5ips <- read_tsv("featurecounts/featurecounts_17/gsea/neuro.tsv")

c5 <- subset(c5ips, padj <0.05) %>% arrange(NES) 
c5 <- c5[grep("motor|senso|neuron|nerv|spina|loco", c5$pathway, ignore.case = TRUE),]
c5$pathway <- tolower(gsub("_"," ",c5$pathway))
c5$pathway <- gsub ("gobp", "GO:BP", c5$pathway)
c5$pathway <- gsub ("gocc", "GO:CC", c5$pathway)
c5$pathway <- gsub ("hp", "HP", c5$pathway)

CairoSVG("Figure1D.svg", width = 10, height = 1.5)
topc5 <- top_n(c5,1,wt = -NES)
bottomc5 <- top_n(c5,3,wt = NES)
full_join(topc5,bottomc5) %>% 
  mutate(pathway = fct_reorder(pathway, dplyr::desc(NES))) %>% 
  mutate(sign = if_else(NES >= 0, 'MN', 'SN')) %>%
  ggplot(aes(NES, pathway, fill=sign))+geom_bar(stat="identity") +
  scale_fill_manual(values = c("MN" = "red", "SN" = "blue")) + xlab("Normalised Enrichment Score (NES)")+
  guides(fill=guide_legend(title="Enriched in")) + theme_cowplot() + theme(axis.title.y = element_blank(), legend.position="right")
dev.off()

