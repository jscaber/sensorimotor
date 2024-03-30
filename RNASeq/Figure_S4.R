library(tidyverse)
library(DESeq2)
library(vidger)
library(biomaRt)
library(viridis)
library(cowplot)
library(Cairo)
library(patchwork)
library(RColorBrewer)
library(pheatmap)




mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2022.archive.ensembl.org")
getmart <- function(values){
  data<- getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id", "external_gene_name", "description",'chromosome_name',
                  'start_position', 'end_position'),
    values= values,
    mart= mart,
    useCache = FALSE)
  data$description <- gsub("\t", "", data$description)
  return(data)
}

makeTPMtable  <- function(genelist, abundance, design, contrast){
  genelist.df <- getmart(genelist)
  genelist.df <- genelist.df[!duplicated(genelist.df[,1]),]
  rownames(genelist.df) <- genelist.df$ensembl_gene_id
  genelist.names <- genelist.df[genelist,]$external_gene_name
  genelist.names[is.na(genelist.names)] <- genelist[is.na(genelist.names)]
  dftemp <- as_tibble(t(abundance[genelist,]), rownames = "track")
  dftemp <- dftemp %>% rename_at(vars(all_of(genelist)), ~ genelist.names)
  dftemp$contrast <- design[dftemp$track,][,contrast]
  return(dftemp)
}

# plotTPMs function
plotTPMs <- function(dftemp, contrast_name){
  dftemp %>% 
    gather(key = "var", value="value", -contrast, -track) %>% 
    mutate(var = factor(var, levels=unique(var))) %>%
    ggplot(aes(x = contrast, y = value, color = contrast)) +
    geom_point(position = position_jitter(w = 0.15, h = 0)) +
    facet_wrap(~ var, scales = "free",nrow = 1)+
    ylab("normalised counts") + guides(color = "none")+ theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank()) 
}

full<- read.table("featurecounts/featurecounts_24/results_full.tsv")
goi <- rownames(full[full$symbol %in% c("PIEZO2","AVIL","POU4F1","TRPV2","MNX1","NKX6-1","CHAT", "SLC5A7", "SCN9A","SCN10A"),])
experiment <- readRDS("featurecounts/featurecounts_24/experiment_out.rds")
dds <- estimateSizeFactors(experiment)
dds <- DESeq(experiment, betaPrior = FALSE)
normcounts <- counts(dds, normalized=TRUE)

vsd<- vst(dds)

CairoSVG("FigureS4A.svg", width = 10, height = 10,pointsize = 2)
vsScatterMatrix(data = dds,type = "deseq", d.factor = "group", title = FALSE, grid = TRUE,)
dev.off()

Cairo("FigureS4A.png", width = 10, height = 10, units = "in",dpi = 300)
vsScatterMatrix(data = dds,type = "deseq", d.factor = "group", title = FALSE, grid = TRUE,)
dev.off()

pca = prcomp(t(assay(vsd)))
colData(dds)
variable.group <- colData(dds)[, "group"]
variable.celltype <- colData(dds)[, "celltype"]
variable.source <- colData(dds)[, "source"]
variable.subcell <- colData(dds)[, "subcell"]

percentVar <- round(100 * summary(pca)$importance[2,])
scores <- data.frame(variable.group,variable.celltype,variable.source, pca$x[,1:8])

CairoSVG("FigureS4PCA.svg", width = 14, height = 10,pointsize = 16)
print(qplot(x=PC1, y=PC2, data=scores, colour=variable.group, shape=variable.source) +
        theme(legend.position="right") +  
        labs(colour="group", x=paste0("PC1 (", percentVar[1],"% of variance)"),
             y=paste0("PC2 (", percentVar[2],"% of variance)")) + geom_point(size=3, stroke = 2))+
  scale_colour_manual(values=c("violetred","slateblue","green","red","blue","black","orange","grey"))+
  scale_shape_manual(values=c("triangle", "plus", "circle"),) + theme_bw()+
  labs(color="Cell type", shape = "Tissue origin")  + theme_cowplot(font_size = 24)
dev.off()

goi <- rownames(full[full$symbol %in% c("MNX1","NKX6-1","NEFH","CHAT", "SLC5A7","GRIN2A", "GRIN2B", "MKI67", "GFAP","S100B", "AQP4","ISL1","OLIG2","PAX6","NGFR","MAPT"),])
mat <- assay(vsd)[goi, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- full[rownames(full) %in% rownames(mat),]$symbol
colnames(mat) <- paste(colData(dds)[,c("group")], vsd$track, sep = " - " )
df <- as.data.frame(colData(dds)[,c("group")])
#df[df$source == "adult human (pseudobulk)",]$source <- "sn-seq MNs"
#df[df$source == "adult human (bulk)",]$source <- "LCM-seq MNs"
#df[df$source == "in vitro (bulk)",]$source <- "iPS-MNs"
rownames(df) <- paste(colData(dds)[,c("group")], vsd$track, sep = " - " )
colnames(df) <- c("Group")
celltype <- c("black", "grey40", "green", "darkgreen", "orange", "red", "hotpink","magenta")
names(celltype) <- c("aMN", "aMNsc", "UCL-MN-d21", "UCL-MN-d35", "Exeter-MN", "iPSC-MN", "Ulm-MN-d21", "Ulm-MN-d42")
anno.colors <- list(celltype)
names(anno.colors) <- c("Group")
colors <- colorRampPalette( brewer.pal(9, "Oranges") )(255)
CairoSVG("FigureS4HEAT.svg", width = 7, height = 5)
pheatmap(mat, annotation_col=df, scale = "row", annotation_colors = anno.colors,
         cluster_cols = TRUE, show_colnames = FALSE, main = "Motor and sensory neuron markers"
         )
dev.off()

###### SAMPLE DISTANCE MATRIX##########
colors <- colorRampPalette( rev(brewer.pal(9, "Oranges" )) )(255)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colData(vsd)[,c("group")]
colnames(sampleDistMatrix) <- NULL
CairoSVG("FigureS4B.svg", width = 7, height = 5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample distance matrix",
         col = colors,cluster_cols = TRUE, 
         annotation_names_row = FALSE, annotation_colors = anno.colors, show_rownames = TRUE, )
dev.off()



###### TPM COMPARISON - INDIVIDUAL LINES #####
genes<- c("MNX1", "ISL1", "NKX6-1","NKX6-2","CHAT", "SLC5A7", "NEFH", "GRIN2B","S100B", "MKI67","GFAP", "OLIG2")
goi <- full[full$symbol %in% genes,]
goi = rownames(goi[match(genes, goi$symbol),])
tpmtable <- makeTPMtable(goi, normcounts, colData(dds), "group")
p1 <- tpmtable %>% 
  gather(key = "var", value="value", -contrast, -track) %>% 
  mutate(var = factor(var, levels=unique(var))) %>%
  ggplot(aes(x = contrast, y = value, color = contrast)) +
  geom_point(position = position_jitter(w = 0.15, h = 0)) +
  facet_wrap(~ var, scales = "free",nrow = 3)+
  ylab("normalised counts") + guides(color = "none")+ theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank()) 
CairoSVG("FigureS4C.svg", width = 13, height = 9)
print(p1)
dev.off()

###### TPM COMPARISON - OTHER #####
experiment<- readRDS("featurecounts/featurecounts_10/experiment_out.rds")
dds <- estimateSizeFactors(experiment)
dds <- DESeq(experiment, betaPrior = FALSE)
normcounts <- counts(dds, normalized=TRUE)

genes1<- c("MNX1", "ISL1", "NKX6-1","NKX6-2")
genes2<- c("CHAT", "SLC5A7", "NEFH", "GRIN2B")
genes3<- c("S100B", "MKI67","GFAP", "OLIG2")

geneplot <- function(genes, legend=TRUE, plottitle=""){
  goi <- full[full$symbol %in% genes,]
  goi = rownames(goi[match(genes, goi$symbol),])
  tpmtable <- makeTPMtable(goi, normcounts, colData(dds), "group")
  tpmtable$track
  tpmtable$group2 <- c("Ulm.d42", "Ulm.d42", "Exeter", "Exeter", "UCL.d35", "UCL.d35","Oxford(CTR)","Oxford(CTR)","Oxford(CTR)","Oxford(C9)","Oxford(C9)","Oxford(C9)")
  tpmtable$contrast <- c("Others", "Others", "Others", "Others", "Others", "Others"," This\nStudy"," This\nStudy"," This\nStudy"," This\nStudy"," This\nStudy"," This\nStudy")
  p1 <- tpmtable %>% 
    gather(key = "var", value="value", -contrast, -track, -group2) %>% 
    mutate(var = factor(var, levels=unique(var))) %>%
    ggplot(aes(x = contrast, y = value, color = contrast, shape = group2)) +
    geom_point(position = position_jitter(w = 0.15, h = 0)) +
    facet_wrap(~ var, scales = "free",nrow = 2)+ ggtitle(plottitle)+
    ylab("normalised counts") + guides(color = "none")+ theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank(), legend.title = element_blank()) 
  if(legend == FALSE){
    p1 = p1 + theme(legend.position="none")
  }
  return(p1)
}

p1 <- geneplot(genes1, legend= FALSE, plottitle = "Transcription factors")
p2 <- geneplot(genes2, legend= FALSE, plottitle = "Functional MN genes")
p3 <- geneplot(genes3, plottitle = "Glial/other cell markers")

p1 + plot_spacer() + p2 + plot_spacer() + p3+   plot_layout(widths = c(1,0.2, 1,0.2,1))

CairoSVG("FigureS4C_alt.svg", width = 13, height = 6)
p1 + plot_spacer() + p2 + plot_spacer() + p3+   plot_layout(widths = c(1,0.2, 1,0.2,1))
dev.off()


