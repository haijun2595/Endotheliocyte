library(scRNAtoolVis)
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
table(Endo_merge_after$subcelltype)
table(Endo_merge_after$subcelltype)
new.cluster.ids <- c("Venous_ECs", "Capillary_ECs", "LUM_ECs", "Tip_like_ECs", "PRRX1_ECs", "Lymphatic_ECs", "Proli_ECs")
names(new.cluster.ids) <- levels(Endo_merge_after)
Endo_merge_after <- RenameIdents(Endo_merge_after, new.cluster.ids)
Endo_merge_after@meta.data$subcelltype <- Endo_merge_after@active.ident
view(Endo_merge_after@meta.data)
table(Endo_merge_after$subcelltype)
table(Endo_merge_after$orig.ident)
table(Endo_merge_after$orig.ident)
Idents(Endo_merge_after) <- "orig.ident"
Endo_merge_after@meta.data$group = ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                            c("lungmeta10","lungmeta17"),"lung",
                                          ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                                   c("NT01","NT03","NT07","NT08"),"NT","PT"))

table(Endo_merge_after$group)
Idents(Endo_merge_after) <- "group" 
###########################
############
library(ggplot2)
library(Seurat)
library(pheatmap)
library(scRNAtoolVis)
library(magrittr)
library(dplyr)
library(scRNAtoolVis)
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
############################################
table(Endo_merge_after$subcelltype)
Idents(Endo_merge_after) <-"group"

p<-VlnPlot(Endo_merge_after,features = "MCAM",pt.size = 0.1)
p
ggsave(p,file="MCAM_Endo.pdf",width = 6,height = 4)

Idents(Endo_merge_after) <-"subcelltype"

Tip_like_ECs <- subset(Endo_merge_after,ident=c("Tip_like_ECs"))
Idents(Tip_like_ECs) <-"group"
p<-VlnPlot(Tip_like_ECs,features = "MCAM",pt.size = 0)
ggsave(p,file="MCAM_TIP.pdf",width = 4,height = 4)

table(Endo_merge_after$subcelltype)
Idents(Endo_merge_after) <-"subcelltype"
markers <- FindAllMarkers(Endo_merge_after,
                          logfc.threshold = 0.25,
                          min.pct = 0.25,
                          only.pos = T)
table(Endo_merge_after$subcelltype)
head(Endo_merge_after)

sig_markers <- markers %>%
  group_by(cluster)%>%
  top_n(n = 5, wt = avg_log2FC)
head(sig_markers)


genes <- unique(sig_markers$gene)

aver_dt <- AverageExpression(Endo_merge_after,
                             features = genes,
                             group.by = c('subcelltype'),  # Include both subcelltype and group
                             # group.by = 'subcelltype' 
                             slot = 'data')

aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt[1:3,1:3]

##########################################################
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
         cluster_cols = FALSE,  #
         annotation_col = cell_anno,  #
         annotation_colors = anno_col,  #
         color = mycol,  # 
         border_color = 'white',
         angle_col = "45") 
dev.off()
table(Endo_merge_after$subcelltype)
###################################
table(Endo_merge_after$subcelltype)
Idents(Endo_merge_after) <- "subcelltype"
Tipcell <- subset(Endo_merge_after,ident=c("Tip_like_ECs"))
table(Tipcell$subcelltype)
Tipcell$subcelltype <- Tipcell@active.ident
table(Tipcell$group)
ibrary(ClusterGVis)
library(org.Hs.eg.db)
library(ComplexHeatmap)
Idents(Tipcell) <- "group"
markers.all <- Seurat::FindAllMarkers(Tipcell,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)
markers <- markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

head(pbmc.markers)
st.data <- prepareDataFromscRNA(object = Tipcell,
                                diffData = markers,
                                showAverage = TRUE)
str(st.data)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
head(enrich)
markGenes = unique(markers$gene)[sample(1:length(unique(markers$gene)),40,
                                        replace = F)]
pdf(file = "Gene.pdf",height = 5,width=8)
visCluster(object = st.data,
           plot.type = "line")
dev.off()

my_colors <- c(
  "#FF0000", "darkblue","#008080","#800080", "#FF7F50","#0000FF", "#6495ED","#4682B4", "#FF00FF", "#00FFFF",
  "#800000", "#808000", "#008000",  "#808080", 
  "#C0C0C0", "#FFA500", "#A52A2A", "#DEB887", "#5F9EA0", "#7FFF00",
  "#D2691E",   "#DC143C", "#00FA9A"
)
repeated_colors <- rep(my_colors, each = 5)
go.col <- my_colors
go.col <- repeated_colors[1:15]  
groups <- unique(Tipcell$group)
subcelltype_colors <- data.frame(
  group = groups,  
  color = c("#5DCADB", 
            "#EBC9DF",
            "#F292BC") 
)
subcelltype_color_list <- setNames(subcelltype_colors$color, subcelltype_colors$group)
set.seed(5201314)
pdf('sc2.pdf',height = 8,width = 12,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:7),
           go.col = go.col,
           add.bar = F,
           sample.col = subcelltype_color_list)
dev.off()
#####################################################
rm(list=ls())
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(tidyverse)
library(stringr)
library(harmony)
library(ggplot2)
library(ggpubr)
library(AUCell) 
library(clusterProfiler)
###################################
scRNA <- Endo_merge_after
gene_list <- read.table("Angiogenesis.txt", header = T)
###########################AddModuleScore
DefaultAssay(scRNA) <- "RNA"

scRNA <- AddModuleScore(scRNA,
                        features = gene_list, 
                        ctrl = 100, 
                        name = "AddModuleScore")
colnames(scRNA@meta.data)[11] <- 'Angioenesis_Score' 
p1 <- VlnPlot(scRNA, 
              features = 'Angioenesis_Score', 
              pt.size = 0, 
              adjust = 2,
              group.by = "subcelltype") +
  scale_fill_manual(values = c("#384998","#4988A2","#A65697", 
                               "#E9BF1F","#2C702C","#DD9A63","#E67A8A")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
p1


avg_scores <- tapply(scRNA@meta.data$Angioenesis_Score, scRNA@meta.data$subcelltype, mean)
avg_scores_df <- data.frame(subcelltype = names(avg_scores), avg_score = avg_scores)

p2 <- p1 + geom_text(data = avg_scores_df, 
                     aes(x = subcelltype, y = avg_score, label = round(avg_score, 4)),
                     inherit.aes = FALSE, 
                     position = position_nudge(y = 0.0), 
                     size = 4)
p2
ggsave(filename = "Angiogenesis_Score.pdf",p2,width=6,height=6)
dev.off()
####################AUcell
expr_matrix <- scRNA@assays$RNA@counts
geneSets <- gene_list
cells_rankings <- AUCell_buildRankings(expr_matrix)
aucell_res <- AUCell_calcAUC(geneSets, cells_rankings)  
auc_scores <- aucell_res@assays@data$AUC
auc_scores_vector <- as.vector(auc_scores)
scRNA[["Angioenesis_Score"]] <- auc_scores_vector

p3 <- VlnPlot(scRNA, 
              features = 'Angioenesis_Score', 
              pt.size = 0, 
              adjust = 2,
              group.by = "subcelltype")+scale_fill_manual(values = c(
                "#384998","#4988A2","#A65697","#E9BF1F","#2C702C","#DD9A63","#E67A8A"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
p3
dev.off()
average_auc <- tapply(scRNA@meta.data$Angioenesis_Score, scRNA@meta.data$subcelltype, mean)
avg_scores_df <- data.frame(subcelltype = names(average_auc), avg_score = average_auc)

p4 <- p3 + geom_text(data = avg_scores_df, 
                     aes(x = subcelltype, y = avg_score, label = round(avg_score, 4)),
                     inherit.aes = FALSE, 
                     position = position_nudge(y = 0.0), 
                     size = 4)

p4 <- p4 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4
ggsave(filename = "Angioenesis_AUC_Score.pdf",p4,width = 8,height = 5)

