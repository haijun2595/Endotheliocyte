rm(list = ls())
library(Seurat)
library(tidyverse)
library(patchwork)
table(OS.combined_harmony$orig.ident)
OS.combined_harmony@meta.data$group = ifelse(OS.combined_harmony@meta.data$orig.ident %in%
                                      c("lungmeta10","lungmeta17"),"lung",
                                    ifelse(OS.combined_harmony@meta.data$orig.ident %in%
                                             c("NT01","NT03","NT07","NT08"),"NT","PT"))

Idents(OS.combined_harmony) <- "celltype"
DimPlot(OS.combined_harmony,split.by = 'group',raster=FALSE)
view(OS.combined_harmony@meta.data)
OS.combined_harmony$group <- as.factor(OS.combined_harmony$group)  
table(OS.combined_harmony$group)

######
Idents(OS.combined_harmony) <- "group"
Tumor <- subset(OS.combined_harmony,idents = c("PT","lung"))

table(Tumor$celltype)
Idents(Tumor) <- "celltype"
scRNA <- subset(Tumor,idents = c("Osteoblasts","Myeloid cells","Proliferative cells","T/NK cells","Pericytes","Osteoclasts"))
table(Endo_merge_after$orig.ident)
Endo_merge_after@meta.data$group = ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                               c("lungmeta10","lungmeta17"),"lung",
                                             ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                                      c("NT01","NT03","NT07","NT08"),"NT","PT"))

Idents(Endo_merge_after) <- "subcelltype" 
view(Endo_merge_after@meta.data)
Endo_merge_after$group <- as.factor(Endo_merge_after$group)  ###
table(Endo_merge_after$group)

Idents(Endo_merge_after) <- "group"
Endo_Tumor <- subset(Endo_merge_after,idents = c("PT","lung"))
table(Endo_Tumor$orig.ident)
table(Endo_Tumor$subcelltype)
table(Endo_Tumor$group)
table(Endo_Tumor$subcelltype) 
Idents(Endo_Tumor) <- "subcelltype"
MCAM_Tumor <- subset(Endo_Tumor,idents = c("S0_Endo_MCAM"))
table(MCAM_Tumor$subcelltype) 
table(MCAM_Tumor$group) 
view(MCAM_Tumor@meta.data)

MCAM_Tumor@meta.data$subcelltype <- gsub("S0_Endo_MCAM", "tip_like_ECs", MCAM_Tumor@meta.data$subcelltype)
table(MCAM_Tumor$subcelltype) 
Idents(MCAM_Tumor) <-"subcelltype"
newLabels<- c("tip_like_ECs"="tip_like_ECs")
names(newLabels)=levels(MCAM_Tumor)
MCAM_Tumor=RenameIdents(MCAM_Tumor, newLabels)
MCAM_Tumor@meta.data$celltype1 <- MCAM_Tumor@active.ident
View(MCAM_Tumor@meta.data)
table(MCAM_Tumor$subcelltype) 
table(MCAM_Tumor$celltype1) 
rm(Endo_Tumor)
scRNA@meta.data$celltype <- as.factor(scRNA@meta.data$celltype)
table(scRNA$celltype)
Idents(scRNA) <-"celltype"
scRNA@meta.data$celltype <- scRNA@active.ident
newLabels<- c("Osteoblasts"="Osteoblasts",
              "Myeloid cells"="Myeloid cells",
              "Proliferative cells"="Proliferative cells",
              "T/NK cells"="T/NK cells",
              "Pericytes"="Pericytes",
              "Osteoclasts"="Osteoclasts") 
names(newLabels)=levels(scRNA)
scRNA=RenameIdents(scRNA, newLabels)
scRNA@meta.data$celltype1 <- scRNA@active.ident
table(scRNA$celltype1)
View(scRNA@meta.data)




Tumor_MCAM_celltype<- merge(scRNA, y = list(MCAM_Tumor), add.cell.ids = c("scRNA","MCAM_Tumor"), project = "merge")

saveRDS(Tumor_MCAM_celltype, "Tumor_MCAM_celltype.rds")  ###
View(Tumor_MCAM_celltype@meta.data)
table(Tumor_MCAM_celltype$orig.ident)
table(Tumor_MCAM_celltype$group)


geneID_10x <- read.table('features.tsv', header=F, sep='\t')
geneID_10x<-geneID_10x[,-3]
names(geneID_10x) <- c('Ensembl','Gene')

Tumor_MCAM_celltype <- readRDS("Tumor_MCAM_celltype.rds")
library(Seurat)
table(Tumor_MCAM_celltype@meta.data$celltype1)

sp1<-Tumor_MCAM_celltype
rm(Tumor_MCAM_celltype)
sp1[["celltype1"]]<-sp1@active.ident
sp1_counts <- data.frame(sp1@assays$RNA@data, check.names = F,stringsAsFactors = F)
sp1_counts <- data.frame(Gene=rownames(sp1_counts), sp1_counts,check.names = F,stringsAsFactors = F)
sp1_counts <- inner_join(geneID_10x,sp1_counts)
sp1_counts <- sp1_counts[,-which(names(sp1_counts)=="Gene")]
sp1_meta <- data.frame(Cell=rownames(sp1@meta.data), cell_type=sp1@meta.data$celltype1,stringsAsFactors = F)
write.table(sp1_counts, "sp1_counts.txt", row.names=F, sep='\t')
write.table(sp1_meta, "sp1_meta.txt", row.names=F, sep='\t')



conda activate cellphone
cd /path
cellphonedb method statistical_analysis sp1_meta.txt sp1_counts.txt --threads 43

cellphonedb plot heatmap_plot sp1_meta.txt \
--pvalues-path out/pvalues.txt \
--output-path out

library(ktplots)
library(devtools)
#library(remotes)

devtools::install_github("sqjin/CellChat")
#remotes::install_github('zktuong/ktplots', dependencies = TRUE)
count_net <- read.delim("count_network.txt", check.names = FALSE)
inter_net<-read.delim("interaction_count.txt",check.names=FALSE)
pvalues <- read.delim("pvalues.txt", check.names = FALSE)
means <- read.delim("means.txt", check.names = FALSE)
sig.means <- read.delim("significant_means.txt", check.names = FALSE)
library(CellChat)

deconvoluted <- read.delim("deconvoluted.txt", check.names = FALSE)

pvals <- read.delim("pvalues.txt", check.names = FALSE)
scRNA <- Tumor_MCAM_celltype
table(scRNA$celltype1)
means <- read.delim(paste0("means.txt"), check.names = FALSE)

plot_cpdb(cell_type1 = 'tip_like_ECs', cell_type2 = "", scdata = scRNA,       
          idents = 'celltype1', means = means, pvals = pvals,       
          gene.family = 'costimulatory',#highlight = "blue",       
          keep_significant_only=T) +  
  theme(axis.text  = element_text(size = 10, color = 'black'))

plot_cpdb(
  cell_type1 = 'tip_like_ECs', 
  cell_type2 = '', 
  scdata = scRNA, 
  #idents = 'celltype1', 
  means = means,  
  pvals = pvals,  
  celltype_key = 'celltype1',  
 # gene.family = 'costimulatory', 
  highlight_col = 'blue',  
  keep_significant_only = T
) + 
  theme(axis.text = element_text(size = 10, color = 'black'))




################
p <- plot_cpdb(
  scdata = scRNA,  
  cell_type1 = "tip_like_ECs",
  cell_type2 = "Osteoblast",
  means = means,
  pvals = pvals,
  celltype_key = "celltype1",
  col_option = "#F9A46B",
  gene_family = NULL,  
  genes = c("CD46", "LTBR","VEGFA","KDR","MIF","EFNB1"),  
  highlight_col = "blue",  
  keep_significant_only = TRUE  
) + 
  theme(
    axis.text = element_text(size = 10, color = 'black'),  
    panel.grid.major = element_line(color = "gray80", size = 0.1),  
    panel.grid.minor = element_line(color = "gray90", size = 0.1), 
    panel.background = element_blank(),  
    plot.background = element_blank()  
  )  
  #coord_flip()
p
















