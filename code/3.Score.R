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

