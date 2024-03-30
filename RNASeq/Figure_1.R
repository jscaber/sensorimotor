library(tidyverse)
library(DESeq2)
library(pheatmap)
library(cowplot)
library(EnsDb.Hsapiens.v86)
library(rhdf5)
library(sva)
library(ggplot2)
library(biomaRt)
library(Cairo)
library(RColorBrewer)
library(ggforce)
library(ggrepel)
library(patchwork)
library(goseq)
library(VennDiagram)
library(lattice)
library(gridExtra)


### LOAD AND TRANSFORM DATA###
full <- read.table("featurecounts/featurecounts_17/results_full.tsv")
goi <- rownames(full[full$symbol %in% c("PIEZO2","AVIL","POU4F1","MNX1","NKX6-1","CHAT", "SLC5A7", "SCN9A","SCN10A"),])
experiment <- readRDS("featurecounts/featurecounts_17/experiment_out.rds")

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



##### HEATMAP OF GENES OF INTEREST (SUPERVISED CLUSTERING) ######
goi <- rownames(full[full$symbol %in% c("PIEZO2","POU4F1","ISL1","MNX1","SCGN","NKX6-1","NEFH","CHAT", "SLC5A7","SCN10A"),])
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
CairoSVG("Figure1E.svg", width = 5, height = 3)
pheatmap(mat, annotation_col=df,annotation_colors = anno.colors, col=colors, scale = "row", cluster_cols = TRUE, show_colnames = FALSE, main = "Motor and Sensory Neuron Markers") 
dev.off()
p3 <-pheatmap(mat, annotation_col=df,annotation_colors = anno.colors, col=colors, scale = "row", cluster_cols = TRUE, show_colnames = FALSE, main = "Motor and Sensory Neuron Markers") 

##### VOLCANO PLOT MN VS SN ######
design(dds) <- formula(~group+sex)
dds$group <- factor(dds$group, levels = c("iPSC-MN","iPSC-SN","iPSC"))
dds$group <- relevel(dds$group, ref = "iPSC-SN")
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds)
res <- results(dds)
res <- results(dds, name="group_iPSC.MN_vs_iPSC.SN", cooksCutoff = FALSE)
resLFC <- lfcShrink(dds, res = res, coef="group_iPSC.MN_vs_iPSC.SN", type="apeglm")
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
  guides(colour = "none") + ylim(0,105) + xlim(-16,16)+
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FoldChange.x, -log(padj.x,10), label = symbol),stat = "identity",
                   size = 4, max.overlaps = 100, min.segment.length = 0) +theme_cowplot(font_size = 20) 

##### COMBINED VOLCANO & HEATMAP ######
CairoSVG("Figure1BC.svg", width = 12, height = 5)
p1+ plot_spacer()+ p2+plot_layout (widths = c(0.8,0.2, 1.3))
dev.off()


##### DIFFERENTIAL EXPRESSION IN C9 ######
dds <- readRDS("featurecounts/featurecounts_19/experiment_out.rds")

design(dds) <- formula(~group+sex)
dds$group <- factor(dds$group, levels = c("C9","CTR"))
dds$group <- relevel(dds$group, ref = "CTR")
dds$group
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, betaPrior=FALSE)
resLFC <- lfcShrink(dds, coef="group_C9_vs_CTR", type="normal")
CairoSVG("Figure1F.svg", width = 4.5, height = 3)
DESeq2::plotMA(resLFC, ylim=c(-2,2), main="C9orf72 vs Control Motor Neurons",colSig = "red", alpha = 0.05) 
dev.off()
resmn<- subset(resLFC, padj <0.05)
resmn$Ensembl <- rownames(resmn)


##### DIFFERENTIAL EXPRESSION IN C9 SN ######
dds <- readRDS("featurecounts/featurecounts_25/experiment_out.rds")

design(dds) <- formula(~group+sex)
dds$group <- factor(dds$group, levels = c("C9","CTR"))
dds$group <- relevel(dds$group, ref = "CTR")
dds$group
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, betaPrior=FALSE)
ressn <- results(dds)
resLFC <- lfcShrink(dds, coef="group_C9_vs_CTR", type="normal")
CairoSVG("Figure1G.svg", width = 4.5, height = 3)
DESeq2::plotMA(resLFC, ylim=c(-2,2), main="C9orf72 vs Control Sensory Neurons",colSig = "blue", alpha = 0.05) 
dev.off()
ressn<- subset(resLFC, padj <0.05)
ressn$Ensembl <- rownames(ressn)



### VENN DIAGRAM C9 MN VS SN######

all <- list(resmn, ressn)
goseq.list <- list(all)
names(goseq.list) = c("all")
tx <- transcriptsBy(EnsDb.Hsapiens.v86)
ld <- median(width(tx))

for(item in names(goseq.list)){
  de1 <- as_tibble(goseq.list[[item]][[1]])
  de2 <- as_tibble(goseq.list[[item]][[2]])
  dejoin <- inner_join(de1,de2, by="Ensembl")

  
  dataframes <- list(de1, de2)
  names(dataframes) <- c("", "")
  x <- lapply(dataframes, function(x1) as_vector(x1[,"Ensembl"]))
  str(x)
  #lapply(calculate.overlap(x), write, , append=TRUE, ncolumns=1000)
  
  if(length(unlist(x)) != 0){
    venn.diagram(x,
                 filename = paste0("Venn_",item,".png"),
                 output = TRUE ,
                 imagetype="png" ,
                 height = 700 ,
                 width = 700 ,
                 resolution = 300,
                 compression = "lzw",
                 lwd = 2,
                 lty = 'blank',
                 fill = c("red","blue"),
                 cex = 1,
                 fontface = "bold",
                 fontfamily = "sans",
                 cat.cex = 1,
                 cat.fontface = "bold",
                 cat.default.pos = "text",
                 cat.fontfamily = "sans"
    )
  }
  
}

