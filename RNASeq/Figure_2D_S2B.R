library(tidyverse)
library(ggpmisc)
library(ggrepel)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(goseq)
library(VennDiagram)
library(RRHO)
library(lattice)
library(gridExtra)
library(RColorBrewer)
library(Cairo)



###Load Data and Set Genes of Interest###
bulk<- read.table("featurecounts/featurecounts_21/results_full.tsv")

bulk <- as_tibble(bulk, rownames = "Ensembl")%>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange < -0 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
ips <- read.table("featurecounts/featurecounts_22/results_full.tsv")

ips <- as_tibble(ips, rownames = "Ensembl")%>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange < -0 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
sc <- read.table("featurecounts/featurecounts_23/results_full.tsv")

sc<- as_tibble(sc, rownames = "Ensembl")%>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange < -0 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
goi <- rownames(bulk[bulk$symbol %in% c("PIEZO2","AVIL","POU4F1","TRPV2","ISL1","MNX1","SCGN","NKX6-1","NEFH","CHAT", "SLC5A7", "SCN9A","SCN10A"),])


### Generate Linear corrrelation between experiments #####
bulk_vs_sc <- inner_join(bulk,sc, by="Ensembl") %>% 
  mutate(sig = case_when(padj.x < 0.05 & padj.y < 0.05 ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(padj.sum = padj.x + padj.y) %>%
  mutate(l2f.sum = log2FoldChange.x + log2FoldChange.y)
top_genes1 <- rbind(arrange(bulk_vs_sc, l2f.sum) %>% tail(10),
                    subset(bulk_vs_sc, Ensembl %in% goi),
                   arrange(bulk_vs_sc, l2f.sum) %>% head(10),
                   arrange(bulk_vs_sc, padj.sum) %>% head(10)) %>% distinct(Ensembl, .keep_all = TRUE)
bulk_vs_sc_sig <- bulk_vs_sc[bulk_vs_sc$padj.x < 0.05 & bulk_vs_sc$padj.y <0.05 & abs(bulk_vs_sc$log2FoldChange.x) > 1 &abs(bulk_vs_sc$log2FoldChange.y) > 1,]
p1 <- ggplot(bulk_vs_sc, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(aes(color = sig), size = 1/5) +theme_bw() +
  stat_poly_line() + scale_color_manual( values = c("gray50", "firebrick3")) + 
  stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                 after_stat(p.value.label),")", sep = ""))) +
  scale_x_continuous(name = "log2FC (lcm-MN vs bulk-DRG)") + scale_y_continuous(name = "log2FC (sn-MN vs sn-DRG)")
p1sig <- ggplot(bulk_vs_sc_sig, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(size = 1/5) + theme_bw()  +
  stat_poly_line() + stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                                    after_stat(p.value.label),")", sep = "")))+
  scale_x_continuous(name = "log2FC (lcm-MN vs bulk-DRG)") + scale_y_continuous(name = "log2FC (sn-MN vs sn-DRG)") +
  geom_label_repel(data = top_genes1,
                   mapping = aes(log2FoldChange.x, log2FoldChange.y, label = symbol.x),stat = "identity",
                   size = 4, max.overlaps = 100) +theme_light(base_size = 20) 

bulk_vs_ips<- inner_join(bulk,ips, by="Ensembl") %>% 
  mutate(sig = case_when(padj.x < 0.05 & padj.y < 0.05 ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(padj.sum = padj.x + padj.y) %>%
  mutate(l2f.sum = log2FoldChange.x + log2FoldChange.y)
top_genes2 <- rbind(arrange(bulk_vs_ips, l2f.sum) %>% tail(10),
                    arrange(bulk_vs_ips, l2f.sum) %>% head(10),
                    arrange(bulk_vs_ips, padj.sum) %>% head(10)) %>% distinct(Ensembl, .keep_all = TRUE)
bulk_vs_ips_sig <- bulk_vs_ips[bulk_vs_ips$padj.x < 0.05 & bulk_vs_ips$padj.y <0.05 & abs(bulk_vs_ips$log2FoldChange.x) > 1 &abs(bulk_vs_ips$log2FoldChange.y) > 1,]
p2 <-ggplot(bulk_vs_ips, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(aes(color = sig), size = 1/5) +theme_bw() +
  stat_poly_line() + scale_color_manual( values = c("gray50", "firebrick3")) + 
  stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                 after_stat(p.value.label),")", sep = ""))) +
  scale_x_continuous(name = "log2FC (lcm-MN vs bulk-DRG)") + scale_y_continuous(name = "log2FC (iPSC-MN vs iPSC-SN)")
p2sig <-ggplot(bulk_vs_ips_sig, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(size = 1/5) + theme_bw()  +
  stat_poly_line() + stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                                    after_stat(p.value.label),")", sep = "")))+
  scale_x_continuous(name = "log2FC (lcm-MN vs bulk-DRG)") + scale_y_continuous(name = "log2FC (iPSC-MN vs iPSC-SN)") +
  geom_label_repel(data = top_genes2,
                   mapping = aes(log2FoldChange.x, log2FoldChange.y, label = symbol.x),stat = "identity",
                   size = 4, max.overlaps = 100) +theme_light(base_size = 20) 
CairoSVG("FigureS2B.svg", width = 6, height = 6)
p2sig
dev.off()


sc_vs_ips<- inner_join(sc,ips, by="Ensembl") %>% 
  mutate(sig = case_when(padj.x < 0.05 & padj.y < 0.05 ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(padj.sum = padj.x + padj.y) %>%
  mutate(l2f.sum = log2FoldChange.x + log2FoldChange.y)
top_genes3 <- rbind(arrange(sc_vs_ips, l2f.sum) %>% tail(10),
                    arrange(sc_vs_ips, l2f.sum) %>% head(10),
                    arrange(sc_vs_ips, padj.sum) %>% head(10)) %>% distinct(Ensembl, .keep_all = TRUE)
sc_vs_ips_sig <- sc_vs_ips[sc_vs_ips$padj.x < 0.05 & sc_vs_ips$padj.y <0.05 & abs(sc_vs_ips$log2FoldChange.x) > 1 &abs(sc_vs_ips$log2FoldChange.y) > 1,]
p3 <- ggplot(sc_vs_ips, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(aes(color = sig), size = 1/5) +theme_bw() +
  stat_poly_line() + scale_color_manual( values = c("gray50", "firebrick3")) + 
  stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                 after_stat(p.value.label),")", sep = ""))) +
  scale_x_continuous(name = "log2FC (sn-MN vs sn-DRG)") + scale_y_continuous(name = "log2FC (iPSC-MN vs iPSC-SN)")
p3sig <- ggplot(sc_vs_ips_sig, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(size = 1/5) + theme_bw()  +
  stat_poly_line() + stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                                    after_stat(p.value.label),")", sep = "")))+
  scale_x_continuous(name = "log2FC (sn-MN vs sn-DRG)") + scale_y_continuous(name = "log2FC (iPSC-MN vs iPSC-SN)") +
  geom_label_repel(data = top_genes3,
                   mapping = aes(log2FoldChange.x, log2FoldChange.y, label = symbol.x),stat = "identity",
                   size = 4, max.overlaps = 100) +theme_light(base_size = 20) 
CairoSVG("Figure2D.svg", width = 6, height = 6)
p3sig
dev.off()