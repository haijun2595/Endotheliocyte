######################################
##############################
library(DESeq2)
library(limma)
library(ggpubr)
library(pheatmap)
library(limma)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(tidyverse)
library(stringr)
library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(bigmemory)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(ggpubr)
library(ggalluvial)
library(reshape2)
diff <-Endo_merge_after
table(diff$group)
table(diff$subcelltype)
table(diff$subcelltype)

diff$subcelltype <- as.factor(diff$subcelltype)
Idents(diff) <- "subcelltype"
new.cluster.ids <- c("Venous_ECs", "Capillary_ECs", "LUM_ECs", "Tip_like_ECs", "PRRX1_ECs", "Lymphatic_ECs", "Proli_ECs")
names(new.cluster.ids) <- levels(diff)
diff <- RenameIdents(diff, new.cluster.ids)
diff@meta.data$subcelltype <- diff@active.ident
Idents(diff) <- "subcelltype"
Endo_MCAM <-subset(diff,idents = c("Tip_like_ECs"))
table(Endo_MCAM$group)
table(Endo_MCAM$subcelltype)
Idents(Endo_MCAM) <- "subcelltype"
Endo_MCAM@meta.data$subcelltype <- Endo_MCAM@active.ident
Idents(Endo_MCAM) <- "group"
Endo_LN <-subset(Endo_MCAM,idents = c("lung","NT"))

Idents(Endo_LN) <- "group"
table(Endo_LN$group)
Endo_LN@meta.data$group <- Endo_LN@active.ident
saveRDS(Endo_LN,file = "Endo_LN.rds")
table(Endo_LN$subcelltype)
table(Endo_LN$group)
LN <- FindMarkers(Endo_LN, ident.1 = "lung", ident.2 = "NT", 
                  slot = "data", 
                  logfc.threshold = 0.25, min.pct = 0.1, 
                  test.use = "wilcox")
#####
resOrdered <- LN[order(LN$p_val_adj),]
resOrdered=as.data.frame(resOrdered)
write.csv(LN, file = "Endo_Lung_vs_NT.csv")

head(resOrdered)
resOrdered$logP <- -log10(resOrdered$p_val_adj)
resOrdered$Group <-"Not-significant"
resOrdered$Group[which((resOrdered$p_val_adj <(0.05))&(resOrdered$avg_log2FC >(0.25)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$p_val_adj <(0.05))&(resOrdered$avg_log2FC <(-0.25)))] <-"Down-regulated"
table(resOrdered$Group)
resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$p_val_adj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)##
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene
#plot
p <-ggscatter(resOrdered,x = "avg_log2FC",y = "logP",
              color = "Group",
              palette = c('#3C5866',"#BBBBBB",'#AB3282'),
              size = 1,
              label = resOrdered$Label,
              font.label = 8,
              repel = T,
              xlab = "avg_log2FC",
              ylab = "-log10(p_val_adj)")  + geom_hline(yintercept = 1.3, linetype = "dashed")+
  geom_vline(xintercept = c(-0.25, 0.25),linetype = "dashed")+ylim(0,85)
p
ggsave("Endo_Lung_vs_NT.pdf",p,width = 6,height = 6)
dev.off()
#######################

if(F){
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library("org.Hs.eg.db")
  library(clusterProfiler)
  Idents(Endo_merge_after) <- "subcelltype" 
  Endo_tiplike <- subset(Endo_merge_after,ident=("Tip_like_ECs"))
  cells_sub <- subset(Endo_tiplike@meta.data, 
                      group %in% c("lung","NT"))
  scRNA_sub <- subset(Endo_tiplike, 
                      cells=row.names(cells_sub))
  
  
  Idents(scRNA_sub) <- "group" 
  scRNA_sub$group <- scRNA_sub@active.ident
  table(scRNA_sub$group)
  table(scRNA_sub$subcelltype)
  Idents(scRNA_sub) <- "subcelltype" 
  scRNA_sub$subcelltype <- scRNA_sub@active.ident
  sub.markers <- FindMarkers(scRNA_sub,group.by = 'group',
                             ident.1 = 'lung',ident.2 = 'NT',logfc.threshold = 0.01)
  sub.markers.sig <- subset(sub.markers, p_val_adj<0.05 & abs(avg_log2FC) >0.01)
    mydata <- data.frame(Gene=rownames(sub.markers.sig),logFC=sub.markers.sig$avg_log2FC)%>%
    arrange(desc(logFC))
  gsea_input <- mydata$logFC
  names(gsea_input) <- mydata$Gene
  geneset <- read.gmt("c5.go.v2023.1.Hs.symbols.gmt")  
  # Run GSEA
  gg <- GSEA(gsea_input, TERM2GENE=geneset,verbose=F,
             pvalueCutoff=0.1, pAdjustMethod = "BH")
  sortgg<- gg[order(gg$NES, decreasing = T),]
  sortgg<- sortgg[sortgg$p.adjust <0.05,]
  write.csv(sortgg,file = "sortgg_lung_VS_PT.csv")
  # plot
  library(enrichplot)
  geneset_plot <- c("GOBP_VASCULATURE_DEVELOPMENT",
                    "GOBP_ENDOTHELIAL_CELL_MIGRATION")
  mycol <- pal_nejm()(8)
  p<-gseaplot2(gg,
               geneSetID = geneset_plot,
               color = mycol[c(1:2)],
               title = 'Tip__Met Lung Vs PC',
               rel_heights = c(2.5, 0.5, 1.0,1.5,2.0),
               pvalue_table = T)
  p
  ggsave(filename = "Tip_Lung_VS_PC.pdf",p,width = 7.5,height = 5)
}




############
 library(ggplot2)
 library(Seurat)
 library(pheatmap)
 library(scRNAtoolVis)
 library(magrittr)
 library(dplyr)
 library(scRNAtoolVis)
 #BiocManager::install("cols4all")
 library(Seurat)
 library(tidyverse)
 library(dplyr)
 library(reshape2)
 library(S4Vectors)
 library(tibble)
 library(SingleCellExperiment)
 library(DESeq2)
 library(ggplot2)
 library(cowplot)
 library(ggrepel)
 library(pheatmap)
 library(Seurat)
 library(dplyr)
 library(pheatmap)
 library(cols4all)
 
 library(ComplexHeatmap)
 library(circlize)
 ############################################33
 table(Endo_merge_after$subcelltype)
 Idents(Endo_merge_after) <-"group"
 
 p<-VlnPlot(Endo_merge_after,features = "MCAM",pt.size = 0.1)
 p
 ggsave(p,file="MCAM_Endo.pdf",width = 6,height = 4)
 

 

 table(Endo_merge_after$subcelltype)
 Idents(Endo_merge_after) <-"subcelltype"
 markers <- FindAllMarkers(Endo_merge_after,
                           logfc.threshold = 0.25,
                           min.pct = 0.25,
                           only.pos = T)
 
 markers_filtered <- markers[!grepl("^RPS|^RPL", markers$gene), ]
 head(markers_filtered)
 
 table(Endo_merge_after$subcelltype)
 head(Endo_merge_after)

 
 
 sig_markers <- markers_filtered %>%
   group_by(cluster)%>%
   top_n(n = 6, wt = avg_log2FC)
 head(sig_markers)
 
 
 genes <- unique(sig_markers$gene)
 
 aver_dt <- AverageExpression(Endo_merge_after,
                              features = genes,
                              group.by = c('subcelltype'),  # Include both subcelltype and group
                              # group.by = 'subcelltype' 
                              slot = 'data')
 
 aver_dt <- as.data.frame(aver_dt$RNA)
 aver_dt[1:3,1:3]
 
 ####################################################
 library(RColorBrewer)
 #mycol <- colorRampPalette(c("#5E3C99", "white", "#E66101"))(50)
 #mycol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(50)
 #mycol <- colorRampPalette(c("#00A8C6", "white", "#F18D05"))(50)
 
 

 mycol <- colorRampPalette(c("#0da9ce", "white", "#e74a32"))(50)
 

 cell_anno <- data.frame(Subcelltype = factor(colnames(aver_dt), levels = unique(colnames(aver_dt))))
 rownames(cell_anno) <- colnames(aver_dt)
 
 anno_col <- list(
   Subcelltype = c("Venous_ECs" = "#384998",
                   "Capillary_ECs" = "#4988A2", 
                   "LUM_ECs" = "#A65697", 
                   "Tip_like_ECs" = "#E9BF1F",
                   "PRRX1_ECs"="#2C702C",
                   "Lymphatic_ECs"="#DD9A63",
                   "Proli_ECs"="#E67A8A" )
 )
 
 
 pdf(file="genepattern.pdf",width = 6,height = 6)
 pheatmap(as.matrix(aver_dt),
          scale = "row",
          cluster_rows = FALSE,  
          cluster_cols = FALSE, 
          annotation_col = cell_anno, 
          annotation_colors = anno_col,
          color = mycol,  
          border_color = 'white',
          angle_col = "45") 
 dev.off()
 table(Endo_merge_after$subcelltype)


