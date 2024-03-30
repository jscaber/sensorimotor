###### NORM EXPRESSION COMPARISON - MN vs SN #####
experiment<- readRDS("featurecounts/featurecounts_26/experiment_out.rds")
dds <- estimateSizeFactors(experiment)
dds <- DESeq(experiment, betaPrior = FALSE)
normcounts <- counts(dds, normalized=TRUE)

genes<- c("TUBB3", "NEFH","CHAT", "ISL1", "NKX6-1", "MNX1", "SCGN", "POU4F1")
goi <- full[full$symbol %in% genes,]
goi = rownames(goi[match(genes, goi$symbol),])
tpmtable <- makeTPMtable(goi, normcounts, colData(dds), "group")
tpmtable$track
tpmtable$contrast <- c("iPSC-MN", "iPSC-MN", "iPSC-MN", "iPSC-MN", "iPSC-MN", "iPSC-MN","iPSC-SN","iPSC-SN","iPSC-SN","iPSC-SN","iPSC-SN","iPSC-SN")
tpmtable$group2 <- c("CTR", "CTR", "CTR", "C9orf72", "C9orf72", "C9orf72","C9orf72","C9orf72","C9orf72","CTR","CTR","CTR")
p1 <- tpmtable %>% 
  gather(key = "var", value="value", -contrast, -track, -group2) %>% 
  mutate(var = factor(var, levels=unique(var))) %>%
  ggplot(aes(x = contrast, y = value, color = contrast, shape = group2)) +
  geom_point(position = position_jitter(w = 0.15, h = 0)) +
  scale_colour_manual(values=c("red","blue"))+
  facet_wrap(~ var, scales = "free",nrow = 1)+ scale_y_continuous(limits = c(0, NA)) +
  ylab("normalised counts") + guides(color = "none")+ theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank(), legend.title = element_blank()) 


CairoSVG("Figure3C.svg", width = 13, height = 3)
p1
dev.off()

###### NORM EXPRESSION COMPARISON - MN vs SN C9 #####
genes<- c("C9orf72", "TARDBP")
goi <- full[full$symbol %in% genes,]
goi = rownames(goi[match(genes, goi$symbol),])
tpmtable <- makeTPMtable(goi, normcounts, colData(dds), "group")
write_tsv(tpmtable,file = "Figure4C RNA.tsv")
