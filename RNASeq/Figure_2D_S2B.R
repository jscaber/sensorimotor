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
full1 <- read.table("coding/featurecounts_1/results_full.tsv")
full1 <- as_tibble(full1, rownames = "Ensembl")%>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange < -0 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
full2 <- read.table("coding/featurecounts_3/results_full.tsv")
full2 <- as_tibble(full2, rownames = "Ensembl")%>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange < -0 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
full3 <- read.table("coding/featurecounts_4/results_full.tsv")
full3<- as_tibble(full3, rownames = "Ensembl")%>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange < -0 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
goi <- rownames(full1[full1$symbol %in% c("PIEZO2","AVIL","POU4F1","TRPV2","ISL1","MNX1","SCGN","NKX6-1","NEFH","CHAT", "SLC5A7", "SCN9A","SCN10A"),])


### Generate Linear corrrelation between experiments #####
lcm_vs_sn <- inner_join(full1,full3, by="Ensembl") %>% 
  mutate(sig = case_when(padj.x < 0.05 & padj.y < 0.05 ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(padj.sum = padj.x + padj.y) %>%
  mutate(l2f.sum = log2FoldChange.x + log2FoldChange.y)
top_genes1 <- rbind(arrange(lcm_vs_sn, l2f.sum) %>% tail(10),
                    subset(lcm_vs_sn, Ensembl %in% goi),
                   arrange(lcm_vs_sn, l2f.sum) %>% head(10),
                   arrange(lcm_vs_sn, padj.sum) %>% head(10)) %>% distinct(Ensembl, .keep_all = TRUE)
lcm_vs_sn_sig <- lcm_vs_sn[lcm_vs_sn$padj.x < 0.05 & lcm_vs_sn$padj.y <0.05 & abs(lcm_vs_sn$log2FoldChange.x) > 1 &abs(lcm_vs_sn$log2FoldChange.y) > 1,]
p1 <- ggplot(lcm_vs_sn, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(aes(color = sig), size = 1/5) +theme_bw() +
  stat_poly_line() + scale_color_manual( values = c("gray50", "firebrick3")) + 
  stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                 after_stat(p.value.label),")", sep = ""))) +
  scale_x_continuous(name = "log2FC (lcm-MN vs bulk-DRG)") + scale_y_continuous(name = "log2FC (sn-MN vs sn-DRG)")
p1sig <- ggplot(lcm_vs_sn_sig, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(size = 1/5) + theme_bw()  +
  stat_poly_line() + stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                                    after_stat(p.value.label),")", sep = "")))+
  scale_x_continuous(name = "log2FC (lcm-MN vs bulk-DRG)") + scale_y_continuous(name = "log2FC (sn-MN vs sn-DRG)") +
  geom_label_repel(data = top_genes1,
                   mapping = aes(log2FoldChange.x, log2FoldChange.y, label = symbol.x),stat = "identity",
                   size = 4, max.overlaps = 100) +theme_light(base_size = 20) 

lcm_vs_drg<- inner_join(full1,full2, by="Ensembl") %>% 
  mutate(sig = case_when(padj.x < 0.05 & padj.y < 0.05 ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(padj.sum = padj.x + padj.y) %>%
  mutate(l2f.sum = log2FoldChange.x + log2FoldChange.y)
top_genes2 <- rbind(arrange(lcm_vs_drg, l2f.sum) %>% tail(10),
                    arrange(lcm_vs_drg, l2f.sum) %>% head(10),
                    arrange(lcm_vs_drg, padj.sum) %>% head(10)) %>% distinct(Ensembl, .keep_all = TRUE)
lcm_vs_drg_sig <- lcm_vs_drg[lcm_vs_drg$padj.x < 0.05 & lcm_vs_drg$padj.y <0.05 & abs(lcm_vs_sn$log2FoldChange.x) > 1 &abs(lcm_vs_sn$log2FoldChange.y) > 1,]
p2 <-ggplot(lcm_vs_drg, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(aes(color = sig), size = 1/5) +theme_bw() +
  stat_poly_line() + scale_color_manual( values = c("gray50", "firebrick3")) + 
  stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                 after_stat(p.value.label),")", sep = ""))) +
  scale_x_continuous(name = "log2FC (lcm-MN vs bulk-DRG)") + scale_y_continuous(name = "log2FC (iPSC-MN vs iPSC-SN)")
p2sig <-ggplot(lcm_vs_drg_sig, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(size = 1/5) + theme_bw()  +
  stat_poly_line() + stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                                    after_stat(p.value.label),")", sep = "")))+
  scale_x_continuous(name = "log2FC (lcm-MN vs bulk-DRG)") + scale_y_continuous(name = "log2FC (iPSC-MN vs iPSC-SN)") +
  geom_label_repel(data = top_genes2,
                   mapping = aes(log2FoldChange.x, log2FoldChange.y, label = symbol.x),stat = "identity",
                   size = 4, max.overlaps = 100) +theme_light(base_size = 20) 

sn_vs_drg<- inner_join(full3,full2, by="Ensembl") %>% 
  mutate(sig = case_when(padj.x < 0.05 & padj.y < 0.05 ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(padj.sum = padj.x + padj.y) %>%
  mutate(l2f.sum = log2FoldChange.x + log2FoldChange.y)
top_genes3 <- rbind(arrange(sn_vs_drg, l2f.sum) %>% tail(10),
                    arrange(sn_vs_drg, l2f.sum) %>% head(10),
                    arrange(sn_vs_drg, padj.sum) %>% head(10)) %>% distinct(Ensembl, .keep_all = TRUE)
sn_vs_drg_sig <- sn_vs_drg[sn_vs_drg$padj.x < 0.05 & sn_vs_drg$padj.y <0.05 & abs(lcm_vs_sn$log2FoldChange.x) > 1 &abs(lcm_vs_sn$log2FoldChange.y) > 1,]
p3 <- ggplot(sn_vs_drg, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(aes(color = sig), size = 1/5) +theme_bw() +
  stat_poly_line() + scale_color_manual( values = c("gray50", "firebrick3")) + 
  stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                 after_stat(p.value.label),")", sep = ""))) +
  scale_x_continuous(name = "log2FC (sn-MN vs sn-DRG)") + scale_y_continuous(name = "log2FC (iPSC-MN vs iPSC-SN)")
p3sig <- ggplot(sn_vs_drg_sig, aes(log2FoldChange.x,log2FoldChange.y)) + geom_point(size = 1/5) + theme_bw()  +
  stat_poly_line() + stat_poly_eq(aes(label = paste("atop(",after_stat(rr.label),",",
                                                    after_stat(p.value.label),")", sep = "")))+
  scale_x_continuous(name = "log2FC (sn-MN vs sn-DRG)") + scale_y_continuous(name = "log2FC (iPSC-MN vs iPSC-SN)") +
  geom_label_repel(data = top_genes3,
                   mapping = aes(log2FoldChange.x, log2FoldChange.y, label = symbol.x),stat = "identity",
                   size = 4, max.overlaps = 100) +theme_light(base_size = 20) 
p1 + plot_spacer()+  p2 +plot_spacer()+ p3+ plot_layout(widths = c(3, 1 ,3,1,3))

CairoSVG("Figure2D.svg", width = 20, height = 6)
p1sig +plot_spacer()+ p2sig +plot_spacer()+ p3sig + plot_layout(widths = c(3, 1 ,3,1,3))
dev.off()


###Investigate overlap with VENN DIAGRAMS ######

sigfull1up <- subset(full1, padj < 0.05 &log2FoldChange > 2)
sigfull1down <- subset(full1, padj < 0.05 &log2FoldChange < -2)
sigfull2up <- subset(full2, padj < 0.05 &log2FoldChange > 2)
sigfull2down <- subset(full2, padj < 0.05 &log2FoldChange < -2)
sigfull3up <- subset(full3, padj < 0.05 &log2FoldChange > 2)
sigfull3down <- subset(full3, padj < 0.05 &log2FoldChange < -2)

up_lfcips <- list(sigfull1up, sigfull2up)
up_lfcsn <- list(sigfull1up, sigfull3up)
up_snips <- list(sigfull3up, sigfull2up)
down_lfcips <- list(sigfull1down, sigfull2down)
down_lfcsn <- list(sigfull1down, sigfull3down)
down_snips <- list(sigfull3down, sigfull2down)
goseq.list <- list(down_snips, down_lfcsn, down_lfcips, up_snips, up_lfcsn, up_lfcips)
names(goseq.list) = c("down_snips", "down_lfcsn", "down_lfcips", "up_snips", "up_lfcsn", "up_lfcips")

tx <- transcriptsBy(EnsDb.Hsapiens.v86)
ld <- median(width(tx))

for(item in names(goseq.list)){
  de1 <- goseq.list[[item]][[1]]
  de2 <- goseq.list[[item]][[2]]
  dejoin <- inner_join(de1,de2, by="Ensembl")
  res.name= as.name(item)
  de.df <- full1
  de.genes <- as.integer(de.df$Ensembl %in% dejoin$Ensembl)
  names(de.genes) = de.df$Ensembl
  temp <- ld[match(de.df$Ensembl, names(ld))]
  pwf=nullp(de.genes,bias.data=temp)
  all = goseq(pwf,"hg38","ensGene", method="Hypergeometric")
  sigall <- all
  names(sigall) <- c("category","pvalue","underrepresented_pvalue","numberDE", "numberTOT", "term", "ontology")
  sigall$pvalue <- p.adjust(sigall$pvalue, method="BH")
  sigall$percent <- sigall$numberDE/sigall$numberTOT
  sigall <- sigall[which(sigall[,2] < 0.05),]
  cats <- sigall$category
  write.table(sigall[,c("category", "pvalue")], file=paste0("GO_",res.name,".tsv"), quote=FALSE, sep='\t', row.names = FALSE)
  write.table(sigall, paste0("GO_annotated_",res.name,".tsv"), quote=FALSE, sep='\t', row.names = FALSE)
  write.table(all, paste0("GO_complete_",res.name,".tsv"), quote=FALSE, sep='\t', row.names = FALSE)
  
  dataframes <- list(de1, de2)
  names(dataframes) <- c("", "")
  x <- lapply(dataframes, function(x1) as_vector(x1[,"Ensembl"]))
  str(x)
  #lapply(calculate.overlap(x), write, , append=TRUE, ncolumns=1000)
  
  if(length(unlist(x)) != 0){
    venn.diagram(x,
                 filename = paste0("Venn_",res.name,".png"),
                 output = TRUE ,
                 imagetype="png" ,
                 height = 700 ,
                 width = 700 ,
                 resolution = 300,
                 compression = "lzw",
                 lwd = 2,
                 lty = 'blank',
                 fill = c("blue","darkblue"),
                 cex = 0.4,
                 fontface = "bold",
                 fontfamily = "sans",
                 cat.cex = 0.4,
                 cat.fontface = "bold",
                 cat.default.pos = "text",
                 cat.fontfamily = "sans"
    )
  }
  
}


#### VISUALIZE AND QUANTIFY OVERLAP WITH RRHO ######
RRHO.lcm_vs_sn <- RRHO(as.data.frame(full1[,c("Ensembl","log2FoldChange")]),
                       as.data.frame(full3[,c("Ensembl","log2FoldChange")]),
                       BY=TRUE, alternative='enrichment', stepsize = 100, labels = c("lcm-MN vs bulk-DRG", "sn-MN vs sn-DRG"), plots = TRUE, outputdir =  "RRHO/")
RRHO.sn_vs_ips <- RRHO(as.data.frame(full3[,c("Ensembl","log2FoldChange")]),
                     as.data.frame(full2[,c("Ensembl","log2FoldChange")]),
                     BY=TRUE, alternative='enrichment', stepsize = 100, labels = c("lcm-MN vs bulk-DRG", "iPSC-MN vs iPSC-SN"), plots = TRUE, outputdir =  "RRHO/")
RRHO.lcm_vs_ips <- RRHO(as.data.frame(full1[,c("Ensembl","log2FoldChange")]),
                        as.data.frame(full2[,c("Ensembl","log2FoldChange")]),
                        BY=TRUE, alternative='enrichment', stepsize = 100, labels = c("sn-MN vs bulk-DRG", "iPSC-MN vs iPSC-SN"), plots = TRUE, outputdir =  "RRHO/")

my.min <-0
my.max <- 1000
my.at <- seq(1, 1000, 9)



q1 <- levelplot(RRHO.lcm_vs_ips$hypermat, col.regions = rev(colorRampPalette(brewer.pal(9, "Spectral"))(100)),
                at=seq(0, 600, length.out=99),xlab='lcm-MN ~ bulk-DRG (Rank)', ylab='iPSMN ~ iPSSN (Rank)')
q2 <- levelplot(RRHO.lcm_vs_sn$hypermat, col.regions = rev(colorRampPalette(brewer.pal(9, "Spectral"))(100)),
                at=seq(0, 600, length.out=99),xlab='lcm-MN ~ bulk-DRG (Rank)', ylab='sn-MN ~ sn-DRG (Rank)')
q3<- levelplot(RRHO.sn_vs_ips$hypermat, col.regions = rev(colorRampPalette(brewer.pal(9, "Spectral"))(100)), 
               at=seq(0, 600, length.out=99),xlab='sn-MN ~ sn-DRG (Rank)', ylab='iPSMN ~ iPSSN (Rank)')

CairoSVG("Figure2E.svg", width = 12.5, height = 4)
grid.arrange(q2, q1, q3, ncol =3)
dev.off()

pval.testing.lcm_vs_sn <- pvalRRHO(RRHO.lcm_vs_sn, 500)
pval.testing.lcm_vs_ips <- pvalRRHO(RRHO.lcm_vs_ips, 500)
pval.testing.sn_vs_ips <- pvalRRHO(RRHO.sn_vs_ips, 500)

xs<- seq(0, 10, length=100)
pval.testing.lcm_vs_ips$pval
r1 <- plot(Vectorize(pval.testing.lcm_vs_ips$FUN.ecdf)(xs)~xs,
           xlab='-log(pvalue)', ylab='ECDF', type='S')
pval.testing.lcm_vs_sn$pval
r2 <- plot(Vectorize(pval.testing.lcm_vs_sn$FUN.ecdf)(xs)~xs,
     xlab='-log(pvalue)', ylab='ECDF', type='S')
pval.testing.sn_vs_ips$pval
r3 <- plot(Vectorize(pval.testing.sn_vs_ips$FUN.ecdf)(xs)~xs,
     xlab='-log(pvalue)', ylab='ECDF', type='S')
save.image("RRHO.RData")
lattice::levelplot(RRHO.lcm_vs_sn$hypermat.by)


full1$metric <- sign(full1$log2FoldChange) * full1$padj
full2$metric <- sign(full2$log2FoldChange) * full2$padj
full3$metric <- sign(full3$log2FoldChange) * full3$padj

pRRHO.lcm_vs_sn <- RRHO(as.data.frame(full1[,c("Ensembl","metric")]),
                       as.data.frame(full3[,c("Ensembl","metric")]),
                       BY=TRUE, alternative='enrichment', stepsize = 100, labels = c("lcm-MN vs bulk-DRG", "sc-MN vs bulk-DRG"), plots = TRUE, outputdir =  "zheng2/")
pRRHO.sn_vs_ips <- RRHO(as.data.frame(full3[,c("Ensembl","metric")]),
                       as.data.frame(full2[,c("Ensembl","metric")]),
                       BY=TRUE, alternative='enrichment', stepsize = 100, labels = c("lcm-MN vs bulk-DRG", "iPSC-MN vs iPSC-SN"), plots = TRUE, outputdir =  "zheng2/")
pRRHO.lcm_vs_ips <- RRHO(as.data.frame(full1[,c("Ensembl","metric")]),
                        as.data.frame(full2[,c("Ensembl","metric")]),
                        BY=TRUE, alternative='enrichment', stepsize = 100, labels = c("sn-MN vs bulk-DRG", "iPSC-MN vs iPSC-SN"), plots = TRUE, outputdir =  "zheng2/")

pq1 <- levelplot(pRRHO.lcm_vs_ips$hypermat, col.regions = rev(colorRampPalette(brewer.pal(9, "Spectral"))(100)), xlab='lcm-MN vs bulk-DRG(Rank)',at=seq(0, 1000, length.out=99), ylab='iPSC-MN vs iPSC-SN (Rank)')
pq2 <- levelplot(pRRHO.lcm_vs_sn$hypermat, col.regions = rev(colorRampPalette(brewer.pal(9, "Spectral"))(100)), xlab='lcm-MN vs bulk-DRG(Rank)',at=seq(0, 1000, length.out=99), ylab='sc-MN vs bulk-DRG (Rank)')
pq3<- levelplot(pRRHO.sn_vs_ips$hypermat, col.regions = rev(colorRampPalette(brewer.pal(9, "Spectral"))(100)), xlab='sc-MN vs bulk-DRG (Rank)',at=seq(0, 1000, length.out=99), ylab='iPSC-MN vs iPSC-SN (Rank)')

grid.arrange(pq2, pq1, pq3, ncol =3)

ppval.testing.lcm_vs_sn <- pvalRRHO(pRRHO.lcm_vs_sn, 500)
ppval.testing.lcm_vs_ips <- pvalRRHO(pRRHO.lcm_vs_ips, 500)
ppval.testing.sn_vs_ips <- pvalRRHO(pRRHO.sn_vs_ips, 500)

xs<- seq(0, 10, length=100)
ppval.testing.lcm_vs_ips$pval
pr1 <- plot(Vectorize(ppval.testing.lcm_vs_ips$FUN.ecdf)(xs)~xs,
           xlab='-log(pvalue)', ylab='ECDF', type='S')
ppval.testing.lcm_vs_sn$pval
pr2 <- plot(Vectorize(ppval.testing.lcm_vs_sn$FUN.ecdf)(xs)~xs,
           xlab='-log(pvalue)', ylab='ECDF', type='S')
ppval.testing.sn_vs_ips$pval
pr3 <- plot(Vectorize(ppval.testing.sn_vs_ips$FUN.ecdf)(xs)~xs,
           xlab='-log(pvalue)', ylab='ECDF', type='S')

save.image("RRHO.RData")

lattice::levelplot(RRHO.lcm_vs_sn$hypermat.by)
