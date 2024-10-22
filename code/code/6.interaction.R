rm(list = ls())
library(Seurat)
library(tidyverse)
library(patchwork)
#############################
setwd("path")
geneID_10x <- read.table('features.tsv', header=F, sep='\t')
geneID_10x<-geneID_10x[,-3]
names(geneID_10x) <- c('Ensembl','Gene')
NT_OB_Endo <- readRDS("/NT_OB_Endo.rds")
sp1<-NT_OB_Endo
rm(NT_OB_Endo)
sp1[["celltype1"]]<-sp1@active.ident
sp1_counts <- data.frame(sp1@assays$RNA@data, check.names = F,stringsAsFactors = F)
sp1_counts <- data.frame(Gene=rownames(sp1_counts), sp1_counts,check.names = F,stringsAsFactors = F)
sp1_counts <- inner_join(geneID_10x,sp1_counts)
sp1_counts <- sp1_counts[,-which(names(sp1_counts)=="Gene")]
sp1_meta <- data.frame(Cell=rownames(sp1@meta.data), cell_type=sp1@meta.data$celltype1,stringsAsFactors = F)
write.table(sp1_counts, "sp1_counts.txt", row.names=F, sep='\t')
write.table(sp1_meta, "sp1_meta.txt", row.names=F, sep='\t')
###################################
conda activate cellphone
cd /path/NT
cellphonedb method statistical_analysis sp1_meta.txt sp1_counts.txt --threads 43
cellphonedb plot heatmap_plot sp1_meta.txt \
--pvalues-path out/pvalues.txt \
--output-path out


conda activate cellphone
cd /path/NT/out
cellphonedb plot dot_plot --means-path means.txt --pvalues-path pvalues.txt \
--rows row.txt --columns col.txt \   
--output-path outdotplot/ --output-name EC.pdf

####################################
library(Seurat)
library(tidyverse)
library(patchwork)
PT_OB_Endo <- readRDS("path/PT_OB_Endo.rds")
table(PT_OB_Endo$subcelltype)
table(PT_OB_Endo$celltype)
table(PT_OB_Endo$celltype1)
table(PT_OB_Endo$group)
table(PT_OB_Endo$orig.ident)
##############################################
###############################
Idents(PT_OB_Endo) <- "orig.ident" 
PT_OB_Endo$Secgroup = ifelse(PT_OB_Endo$orig.ident %in%
                               c("PT02","PT03","PT04","PT07"),"Lung_Metastasis","No_Lung_Metastasis")
Idents(PT_OB_Endo) <- "Secgroup" 
PT_OB_Endo$Secgroup <- as.factor(PT_OB_Endo@active.ident)
table(PT_OB_Endo$Secgroup)
table(PT_OB_Endo$subcelltype)
table(PT_OB_Endo$orig.ident)
table(PT_OB_Endo$celltype1)
PT_OB_Endo$celltype1 <-as.factor(PT_OB_Endo$celltype1)

Idents(PT_OB_Endo) <- "celltype1" 
table(PT_OB_Endo$celltype1)
new.cluster.ids <- c("Osteoblast","Venous_ECs", "Capillary_ECs", "LUM_ECs", "Tip_like_ECs", "PRRX1_ECs", "Lymphatic_ECs", "Proli_ECs")
names(new.cluster.ids) <- levels(PT_OB_Endo)
PT_OB_Endo <- RenameIdents(PT_OB_Endo, new.cluster.ids)
PT_OB_Endo@meta.data$celltype1 <- PT_OB_Endo@active.ident
view(PT_OB_Endo@meta.data)
table(PT_OB_Endo$celltype1)
table(PT_OB_Endo$Secgroup)
################################## 
#######
#############################
Idents(PT_OB_Endo) <- "celltype1" 
table(PT_OB_Endo$celltype1)
Osteoblast_PT<-subset(PT_OB_Endo,idents = c("Osteoblast"))
rm(PT_OB_Endo)
table(Osteoblast_PT$Secgroup)
#############################
Osteoblast_PT <- NormalizeData(Osteoblast_PT, normalization.method = "LogNormalize", scale.factor = 10000)
Osteoblast_PT <- FindVariableFeatures(Osteoblast_PT, selection.method = "vst", nfeatures = 2000)
table(Osteoblast_PT$Secgroup)
#####
all.genes <- rownames(Osteoblast_PT)
Osteoblast_PT <- ScaleData(Osteoblast_PT, features = all.genes)
Osteoblast_PT <- RunPCA(Osteoblast_PT, features = VariableFeatures(object = Osteoblast_PT))
VizDimLoadings(object = Osteoblast_PT, dims = 1:4, reduction = "pca",nfeatures = 20)
ElbowPlot(Osteoblast_PT)
library("harmony")
Osteoblast_PT_harmony <- Osteoblast_PT %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
view(Osteoblast_PT_harmony@meta.data)
Osteoblast_PT_harmony <- Osteoblast_PT_harmony %>%
  RunPCA(npcs = 15, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  RunUMAP(reduction="harmony", dims=1:15)

resolution <- c(0.01,0.03,0.05, 0.07,0.1, 0.15, 0.2, 0.25, 0.3)#  
scRNA<-Osteoblast_PT_harmony 
rm(Osteoblast_PT_harmony)
scRNA <-FindClusters(scRNA,resolution = 0.07)
new.cluster.ids <- c("OB_C0","OB_C1","OB_C2","OB_C3","OB_C4","OB_C5")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA@meta.data$subgroup <- scRNA@active.ident

#view(scRNA@meta.data)
table(scRNA$subgroup)
#################
table(scRNA$Secgroup)
table(scRNA$subgroup)
Idents(scRNA) <- "subgroup"
C1 <- subset(scRNA,idents = c("OB_C1"))
save(C1,file="C1.Rdata")
####################################
##############
rm(scRNA_harmony)
Idents(PT_OB_Endo) <- "celltype1" 
table(PT_OB_Endo$celltype1)
Endo_PT<-subset(PT_OB_Endo,idents = c("Venous_ECs","Capillary_ECs","LUM_ECs","Tip_like_ECs","PRRX1_ECs","Lymphatic_ECs","Proli_ECs"))
rm(Osteoblast_PT_harmony)
table(Endo_PT$celltype1)
view(Endo_PT@meta.data)
view(C1@meta.data)
rm(expr_matrix)
table(Endo_PT$celltype1)
Idents(Endo_PT) <-"celltype1"
Endo_PT$celltype1 <- Endo_PT@active.ident
newLabels<- c("Venous_ECs"="Venous_ECs",  
              "Capillary_ECs"="Capillary_ECs",
              "LUM_ECs"="LUM_ECs",
              "Tip_like_ECs"="Tip_like_ECs",
              "PRRX1_ECs"="PRRX1_ECs",
              "Lymphatic_ECs"="Lymphatic_ECs",
              "Proli_ECs"="Proli_ECs")
names(newLabels)=levels(Endo_PT)
Endo_PT=RenameIdents(Endo_PT, newLabels)
Endo_PT@meta.data$cellchatuse <- Endo_PT@active.ident
View(Endo_PT@meta.data)
table(Endo_PT$cellchatuse)

C1@meta.data$subgroup <- as.factor(C1@meta.data$subgroup)
table(C1$subgroup)
Idents(C1) <-"subgroup"
C1@meta.data$subgroup <- C1@active.ident

newLabels<- c("OB_C1"="OB_C1") #
names(newLabels)=levels(C1)
C1=RenameIdents(C1, newLabels)
C1@meta.data$cellchatuse <- C1@active.ident
View(C1@meta.data)
table(C1$cellchatuse)

C1_OB_Endo<- merge(C1, y = list(Endo_PT), add.cell.ids = c("C1","Endo_PT"), project = "merge")
#saveRDS(OB_Endo, "OB_Endo.rds")  
table(C1_OB_Endo$cellchatuse)
rm(data.input)
#############################
C1_OB_Endo <- NormalizeData(C1_OB_Endo, normalization.method = "LogNormalize", scale.factor = 10000)
C1_OB_Endo <- FindVariableFeatures(C1_OB_Endo, selection.method = "vst", nfeatures = 2000)
#####PCA
all.genes <- rownames(C1_OB_Endo)
C1_OB_Endo <- ScaleData(C1_OB_Endo, features = all.genes)
C1_OB_Endo <- RunPCA(C1_OB_Endo, features = VariableFeatures(object = C1_OB_Endo))
VizDimLoadings(object = C1_OB_Endo, dims = 1:4, reduction = "pca",nfeatures = 20)
ElbowPlot(C1_OB_Endo)

library("harmony")
C1_OB_Endo_harmony <- C1_OB_Endo %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
#load("Endo_harmony1.Rds")
view(C1_OB_Endo_harmony@meta.data)
#metadata<-Endo_harmony@meta.data
#dim(Endo_harmony)
#table(Endo_harmony)
#Idents(Endo_harmony)<-"orig.ident"
C1_OB_Endo_harmony <- C1_OB_Endo_harmony %>%
  RunPCA(npcs = 15, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  RunUMAP(reduction="harmony", dims=1:15)
Idents(C1_OB_Endo_harmony) <- "Secgroup"
table(C1_OB_Endo_harmony$Secgroup)
table(C1_OB_Endo_harmony$orig.ident)
table(C1_OB_Endo_harmony$cellchatuse)
Idents(C1_OB_Endo_harmony) <- "Secgroup"
Metastasis_EC_C1_OB<-subset(C1_OB_Endo_harmony,idents = c("Metastasis"))
table(Metastasis_EC_C1_OB$Secgroup)
table(Metastasis_EC_C1_OB$celltype1)
table(Metastasis_EC_C1_OB$orig.ident)
Idents(Metastasis_EC_C1_OB) <- "cellchatuse"
Metastasis_EC_C1_OB$cellchatuse <- Metastasis_EC_C1_OB@active.ident
#############################
rm(identity)
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(reticulate)
library(ComplexHeatmap)
library(CellChat)
view(Metastasis_EC_OS@meta.data)
scRNA_harmony <- Metastasis_EC_C1_OB
data.input <- GetAssayData(scRNA_harmony, assay="RNA",  slot = "data")
identity <- subset(scRNA_harmony@meta.data, select = "cellchatuse")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "cellchatuse")

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
unique(CellChatDB$interaction$annotation)
##############
cellchat@DB <- CellChatDB # 
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 6)   
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
table(Metastasis_EC_C1_OB$cellchatuse)
table(Metastasis_EC_C1_OB$Secgroup)
view(cellchat@DB)
##################################
df.net.1 <- subsetCommunication(cellchat,slot.name = "netP")
df.net.2 <- subsetCommunication(cellchat )
set.seed(1213)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways

pdf("Cellchat.pdf", width = 8, height = 6)
par(mfrow = c(1,1), xpd=TRUE)  
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
dev.off()
interaction_count <- cellchat@net$count
print(interaction_count)
