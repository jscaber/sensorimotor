library(tidyverse)
library(DESeq2)
library(pheatmap)
library(cowplot)


library(rhdf5)
library(sva)
library(ggplot2)
library(biomaRt)
library(Cairo)
library(RColorBrewer)
library(ggforce)
library(ggrepel)


### LOAD AND TRANSFORM DATA###
full <- read.table("coding/featurecounts_14/results_full.tsv")
goi <- rownames(full[full$symbol %in% c("PIEZO2","AVIL","POU4F1","TRPV2","MNX1","NKX6-1","CHAT", "SLC5A7", "SCN9A","SCN10A"),])
experiment <- readRDS("coding/featurecounts_14/experiment_out.rds")
dds <- experiment
dds <- estimateSizeFactors(dds)
vsd<- vst(dds)

### PRINCIPAL COMPONENT ANALYSIS ###
pca = prcomp(t(assay(vsd)))
variable.group <- colData(dds)[, "group"]
percentVar <- round(100 * summary(pca)$importance[2,])
scores <- data.frame(variable.group, pca$x[,1:8])
print(qplot(x=PC1, y=PC2, data=scores, colour=variable.group) +
        theme(legend.position="right") +  
        labs(colour=NULL, x=paste0("PC1 (", percentVar[1],"% of variance)"),
             y=paste0("PC2 (", percentVar[2],"% of variance)")) + theme_cowplot(font_size = 24) +
        theme(plot.title = element_text(lineheight=1, face="bold"))  + geom_point(size=4)) +
        scale_colour_manual(values=c("green3","red","blue"))
p1 <- qplot(x=PC1, y=PC2, data=scores, colour=variable.group) +
  theme(legend.position="right") +  
  labs(colour=NULL, x=paste0("PC1 (", percentVar[1],"% of variance)"),
       y=paste0("PC2 (", percentVar[2],"% of variance)")) + theme_cowplot(font_size = 20) +
  theme(plot.title = element_text(lineheight=1, face="bold"))  + geom_point(size=4) +
  scale_colour_manual(values=c("green3","red","blue"))


df <- as.data.frame(colData(dds)[,c("group")])
rownames(df) <- colData(dds)$track
colnames(df) <- "Cell Type"

# Heatmap of Top 20 Variable Genes
topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),20)
full <- read.table("salmon_mrna/results_full.tsv")
mat <- assay(vsd)[ topVarGenes, ]
rownames(mat) <- full[rownames(full) %in% rownames(mat),]$symbol
pheatmap(mat, annotation_col=df, cluster_rows=TRUE,cluster_cols = FALSE,fontsize = 14,show_colnames = FALSE)

goi <- rownames(full[full$symbol %in% c("PIEZO2","AVIL","POU4F1","TRPV2","ISL1","MNX1","SCGN","NKX6-1","NEFH","CHAT", "SLC5A7", "SCN9A","SCN10A"),])
mat <- assay(vsd)[goi, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- vsd$track
rownames(mat) <- full[rownames(full) %in% rownames(mat),]$symbol
df <- as.data.frame(colData(dds)[,c("group")])
rownames(df) <- vsd$track
colnames(df) <- c("Cell Type")
celltype <- c("green3","red","blue")
names(celltype) <-c("iPSC","iPSC-MN","iPSC-SN")
anno.colors <- list(c(celltype))
names(anno.colors) = "Cell Type"
colors <- colorRampPalette( brewer.pal(9, "Oranges") )(255)

CairoSVG("Figure1E.svg", width = 7, height = 5)
pheatmap(mat, annotation_col=df,annotation_colors = anno.colors, col=colors, scale = "row", cluster_cols = TRUE, show_colnames = FALSE, main = "Motor and Sensory Neuron Markers") 
dev.off()
p3 <-pheatmap(mat, annotation_col=df,annotation_colors = anno.colors, col=colors, scale = "row", cluster_cols = TRUE, show_colnames = FALSE, main = "Motor and Sensory Neuron Markers") 


design(dds) <- formula(~group+sex)
dds$group <- factor(dds$group, levels = c("iPSC-MN","iPSC-SN","iPSC"))
dds$group <- relevel(dds$group, ref = "iPSC-SN")
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, betaPrior=FALSE)
res <- results(dds)
resLFC <- lfcShrink(dds, coef="group_iPSC.MN_vs_iPSC.SN", type="apeglm")
DESeq2::plotMA(resLFC, ylim=c(-2,2), main="MA plot") + theme_cowplot()
res<-resLFC

res$Ensembl <- rownames(res)
full$Ensembl <- rownames(full)

res <- as_tibble(res) %>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange < 0 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

top_genes1 <- arrange(res[!is.na(res$padj),], log2FoldChange) %>% tail(10)
top_genes3 <- arrange(res[!is.na(res$padj),], log2FoldChange) %>% head(10)
top_genes4 <- arrange(res[!is.na(res$padj),], padj) %>% head(10)
top_genes2 <- res[res$Ensembl %in% goi,]

top_genes <- unique(rbind(top_genes1, top_genes2,top_genes3,top_genes4))
top_genes <- top_genes[top_genes$baseMean >300,]
top_genes <- left_join(keep = FALSE, top_genes, full, by=("Ensembl"))
top_genes


p2 <- ggplot(as.data.frame(res), aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 1/5) +
  xlab(expression("log"[2]*"FoldChange")) + 
  ylab(expression("-log"[10]*"padj")) +
  scale_color_manual( values = c("blue","gray50", "firebrick3")) +
  guides(colour = "none") + ylim(0,180) + xlim(-15,15)+
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FoldChange.x, -log(padj.x,10), label = symbol),stat = "identity",
                   size = 4, max.overlaps = 100) +theme_cowplot(font_size = 20) 


CairoSVG("Figure1BC.svg", width = 14, height = 5)
p1+ plot_spacer()+ p2+plot_layout (widths = c(1,0.3, 1.3))
dev.off()


##### DIFFERENTIAL EXPRESSION IN C9 ######
dds <- readRDS("coding/featurecounts_15/experiment_out.rds")
design(dds) <- formula(~group)
dds$group <- factor(dds$group, levels = c("C9","CTR","CRISPR", "C92"))
dds$group <- relevel(dds$group, ref = "CTR")
dds$group
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, betaPrior=FALSE)
res <- results(dds)
resLFC <- lfcShrink(dds, coef="group_C9_vs_CTR", type="normal")
CairoSVG("Figure1F.svg", width = 7, height = 6)
DESeq2::plotMA(resLFC, ylim=c(-2,2), main="MA plot of C9orf72 vs Control MNs",colSig = "red") 
dev.off()
sigres<- subset(resLFC, padj <0.05)
sigres$Ensembl <- rownames(sigres)
left_join(as_tibble(sigres), full, by="Ensembl")$symbol
