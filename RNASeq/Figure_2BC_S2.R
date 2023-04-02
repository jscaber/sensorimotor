library(tidyverse)
library(data.table)
library(DESeq2)
library(rhdf5)
library(ggplot2)
library(biomaRt)
library(Cairo)
library(pheatmap)
library(RColorBrewer)
library(ggforce)
library(ggrepel)
library(cowplot)
library(viridis)

full <- read.table("coding/featurecounts_1/results_full.tsv")
goi <- rownames(full[full$symbol %in% c("PIEZO2","AVIL","POU4F1","TRPV2","MNX1","NKX6-1","CHAT", "SLC5A7", "SCN9A","SCN10A"),])

### DESEQ and DATA TRANSFOPRM ###
experiment <- readRDS("coding/featurecounts_1/experiment_out.rds")
dds <- estimateSizeFactors(experiment)
dds <- DESeq(experiment, betaPrior = FALSE)
vsd<- vst(dds)


### PRINCIPAL COMPONENT ANALYSIS ###
pca = prcomp(t(assay(vsd)))
colData(dds)
variable.group <- colData(dds)[, "group"]
variable.celltype <- colData(dds)[, "celltype"]
variable.source <- colData(dds)[, "source"]
variable.subcell <- colData(dds)[, "subcell"]

percentVar <- round(100 * summary(pca)$importance[2,])
scores <- data.frame(variable.group,variable.celltype,variable.source, pca$x[,1:8])

CairoSVG("FigureS2.svg", width = 14, height = 10,pointsize = 16)
print(qplot(x=PC1, y=PC2, data=scores, colour=variable.celltype, shape=variable.source) +
        theme(legend.position="right") +  
        labs(colour="group", x=paste0("PC1 (", percentVar[1],"% of variance)"),
             y=paste0("PC2 (", percentVar[2],"% of variance)")) + geom_point(size=3, stroke = 2))+
        scale_colour_manual(values=c("violetred","slateblue","green","red","blue","black"))+
        scale_shape_manual(values=c("triangle", "plus", "circle"),) + theme_bw()+
        labs(color="Cell type", shape = "Tissue origin")  + theme_cowplot(font_size = 24)
dev.off()


### HEATMAP - HIGHEST VARIANCE (exploratory, not used) ###
df <- as.data.frame(colData(dds)[,c("group")])
rownames(df) <- colData(dds)$track
colnames(df) <- "Cell Type"
# Heatmap of Top 20 Variable Genes
topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),20)
mat <- assay(vsd)[ topVarGenes, ]
row.names(mattemp) <- mattemp$ensembl_gene_id
rownames(mat) <- full[rownames(full) %in% rownames(mat),]$symbol
pheatmap(mat, annotation_col=df, cluster_rows=TRUE,cluster_cols = FALSE,fontsize = 14,show_colnames = FALSE)


### HEATMAP - MARKER GENES ###
goi <- rownames(full[full$symbol %in% c("PIEZO2","AVIL","POU4F1","TRPV1","MNX1","SCGN","NKX6-1","NEFH","CHAT", "SLC5A7","GRIN2A", "GRIN2B"),])
mat <- assay(vsd)[goi, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- full[rownames(full) %in% rownames(mat),]$symbol
colnames(mat) <- paste(colData(dds)[,c("group")], vsd$track, sep = " - " )
df <- as.data.frame(colData(dds)[,c("celltype", "source")])
df[df$source == "adult human (pseudobulk)",]$source <- "adult human (sn-seq)"
df[df$source == "adult human (bulk)",]$source <- "adult human (bulk-seq)"
df[df$source == "in vitro (bulk)",]$source <- "in vitro (bulk-seq)"
rownames(df) <- paste(colData(dds)[,c("group")], vsd$track, sep = " - " )
colnames(df) <- c("Cell Type", "Source")
celltype <- c("green", "blue", "slateblue", "red", "violetred")
names(celltype) <- c("iPSC", "iPSC-SN", "Adult SN", "iPSC-MN", "Adult MN")
source <- c("black", "grey39" ,"white")
names(source) <- c("adult human (bulk-seq)", "adult human (sn-seq)", "in vitro (bulk-seq)")
anno.colors <- list(celltype, source)
names(anno.colors) <- c("Cell Type", "Source")
colors <- colorRampPalette( brewer.pal(9, "Oranges") )(255)
CairoSVG("Figure2C.svg", width = 7, height = 5)
pheatmap(mat, annotation_col=df,annotation_colors = anno.colors, scale = "row", 
         cluster_cols = TRUE, show_colnames = FALSE, main = "Motor and sensory neuron markers",
         color = colors)
dev.off()

###### SAMPLE DISTANCE MATRIX##########
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(colData(vsd)[,c("group")], vsd$track, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Oranges")) )(255)
df <- as.data.frame(colData(dds)[,c("celltype", "source")])
df[df$source == "adult human (pseudobulk)",]$source <- "adult human (sn-seq)"
df[df$source == "adult human (bulk)",]$source <- "adult human (bulk-seq)"
df[df$source == "in vitro (bulk)",]$source <- "in vitro (bulk-seq)"
rownames(df) <- paste(colData(dds)[,c("group")], vsd$track, sep = " - " )
colnames(df) <- c("Cell Type", "Source")
celltype <- c("green", "blue", "slateblue", "red", "violetred")
names(celltype) <- c("iPSC", "iPSC-SN", "Adult SN", "iPSC-MN", "Adult MN")
source <- c("black", "grey39" ,"white")
names(source) <- c("adult human (bulk-seq)", "adult human (sn-seq)", "in vitro (bulk-seq)")
anno.colors <- list(celltype, source)
names(anno.colors) <- c("Cell Type", "Source")

CairoSVG("Figure2B.svg", width = 7, height = 5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample distance matrix",
         col = colors,cluster_cols = TRUE, treeheight_col = 0,
         annotation_row=df,
         annotation_names_row = FALSE, annotation_colors = anno.colors, show_rownames = FALSE, )
dev.off()

