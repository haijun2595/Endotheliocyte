library(devtools)
library(Seurat)
library(SPOTlight)
library(tidyverse)
library(patchwork)
library(hdf5r)
library(ggpubr)
library(ggplot2)
MLN07_ST <- Load10X_Spatial(data.dir ="/path/", 
                            filename = "filtered_feature_bc_matrix.h5",
                            assay = "Spatial",
                            slice ="OS")
dir.create("QC")
dir.create("Data")
dir.create("Cluster")

MLN07_ST[["percent.mt"]] <- PercentageFeatureSet(MLN07_ST, pattern = "^mt-")
MLN07_ST[["percent.rb"]] <- PercentageFeatureSet(MLN07_ST, pattern = "^Rp[ls]")
p1 <- VlnPlot(MLN07_ST, features = "nCount_Spatial") + NoLegend() + theme(axis.text.x = element_blank())
p2 <- SpatialFeaturePlot(MLN07_ST, features = "nCount_Spatial") + theme(legend.position = "right")
p <- p1|p2
ggsave("QC/nCount_Spatial_before.pdf", p, width = 10, height = 6)
p1 <- VlnPlot(MLN07_ST, features = "nFeature_Spatial") + NoLegend() + theme(axis.text.x = element_blank())
p2 <- SpatialFeaturePlot(MLN07_ST, features = "nFeature_Spatial") + theme(legend.position = "right")
p <- p1|p2
ggsave("QC/nFeature_Spatial_before.pdf", p, width = 10, height = 6)
p1 <- VlnPlot(MLN07_ST, features = "percent.mt") + NoLegend() + theme(axis.text.x = element_blank())
p2 <- SpatialFeaturePlot(MLN07_ST, features = "percent.mt") + theme(legend.position = "right")
p <- p1|p2
ggsave("QC/percent_mt_before.pdf", p, width = 10, height = 6)
p1 <- VlnPlot(MLN07_ST, features = "percent.rb") + NoLegend() + theme(axis.text.x = element_blank())
p2 <- SpatialFeaturePlot(MLN07_ST, features = "percent.rb") + theme(legend.position = "right")
p <- p1|p2
ggsave("QC/percent_rb_before.pdf", p, width = 10, height = 6) 
## 质控
minCount = 1500
minFeature = 500
maxMT = 25
MLN07_ST <- subset(MLN07_ST, nCount_Spatial>minCount&nFeature_Spatial>minFeature&percent.mt<maxMT)
p1 <- VlnPlot(MLN07_ST, features = "nCount_Spatial") + NoLegend() + theme(axis.text.x = element_blank())
p2 <- SpatialFeaturePlot(MLN07_ST, features = "nCount_Spatial") + theme(legend.position = "right")
p <- p1|p2
ggsave("QC/nCount_Spatial_after.pdf", p, width = 10, height = 6)
p1 <- VlnPlot(MLN07_ST, features = "nFeature_Spatial") + NoLegend() + theme(axis.text.x = element_blank())
p2 <- SpatialFeaturePlot(MLN07_ST, features = "nFeature_Spatial") + theme(legend.position = "right")
p <- p1|p2
ggsave("QC/nFeature_Spatial_after.pdf", p, width = 10, height = 6)
p1 <- VlnPlot(MLN07_ST, features = "percent.mt") + NoLegend() + theme(axis.text.x = element_blank())
p2 <- SpatialFeaturePlot(MLN07_ST, features = "percent.mt") + theme(legend.position = "right")
p <- p1|p2
ggsave("QC/percent_mt_after.pdf", p, width = 10, height = 6)
p1 <- VlnPlot(MLN07_ST, features = "percent.rb") + NoLegend() + theme(axis.text.x = element_blank())
p2 <- SpatialFeaturePlot(MLN07_ST, features = "percent.rb") + theme(legend.position = "right")
p <- p1|p2
ggsave("QC/percent_rb_after.pdf", p, width = 10, height = 6) 
MLN07_ST <- SCTransform(MLN07_ST, assay = "Spatial") 
save(MLN07_ST,file = "Data/MLN07_ST_SCT.RData")
MLN07_ST <- RunPCA(MLN07_ST, assay = "SCT", verbose = FALSE)
ElbowPlot(MLN07_ST, ndims = 50)
MLN07_ST <- FindNeighbors(MLN07_ST, reduction = "pca", dims = 1:10)
MLN07_ST <- FindClusters(MLN07_ST, resolution = 0.4, verbose = FALSE)
MLN07_ST <- RunUMAP(MLN07_ST, reduction = "pca", dims = 1:10)
colors <- c("#f33b99","#E7481B","#1b79af","#179b73","#8da0cb","#bbe173","#970030","#3d4a78","#f6e36d","#ad5f2c")
p1 <- DimPlot(MLN07_ST, reduction = "umap", label = TRUE,cols = colors)
p1
ggsave("Cluster/cluster_0.4.pdf", p1, width = 10, height = 8) 
p2 <- SpatialDimPlot(MLN07_ST, label = TRUE, label.size = 3, group.by = "seurat_clusters", 
                     cols = c("0" = "#f33b99",
                              "1" = "#E7481B",
                              "2" = "#1b79af",
                              "3" =  "#179b73",
                              "4" =  "#8da0cb",
                              "5" = "#bbe173",
                              "6" = "#1b79af",
                              "7" = "#bbe173"))
p2
ggsave("Cluster/louvain_0.4.pdf", p2, width = 12, height = 8) 
p <- p1 + p2
ggsave("Cluster/cluster_louvain_0.4.pdf", p, width = 15, height = 7)  
save(MLN07_ST,file = "Data/MLN07_ST_cluster_0.4.RData")
scRNA <- readRDS("/MLN07.rds")
Idents(scRNA) <- "celltype"
table(scRNA$celltype)
scRNA <- SCTransform(scRNA, verbose = TRUE)

sc.marker <-FindAllMarkers(scRNA,assay = "SCT",slot = "data",only.pos = T,min.pct = 0.25) 
#sc.marker <- FindAllMarkers(scRNA, assay = "RNA", slot = "data", only.pos = TRUE, min.pct = 0.25)


write.csv(sc.marker, file = "Data/marker.csv",row.names = F)
save(sc.marker,file = "Data/sc_marker.RData")
stRNA <- MLN07_ST
st.marker <-FindAllMarkers(stRNA,only.pos = T,min.pct = 0.25) 
write.csv(st.marker, file = "Data/MLN07.csv",row.names = F)
save(st.marker,file = "Data/st_marker.RData")
all.markers = sc.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top100_OB = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>%
  filter(cluster=="Osteoblast")
score_list=list(gene=top100_OB$gene)
stRNA <-AddModuleScore(stRNA,features = score_list,name = "Osteoblast")  
table(scRNA$celltype)
all.markers = sc.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top100_Mye = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>%
  filter(cluster=="Myeloid")

score_list=list(gene=top100_Mye$gene)
stRNA <-AddModuleScore(stRNA,features = score_list,name = "Myeloid_cell") 
all.markers = sc.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top100_T = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>%
  filter(cluster=="T cell")

score_list=list(gene=top100_T$gene)
stRNA <-AddModuleScore(stRNA,features = score_list,name = "T_cell")  
all.markers = sc.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top100_B = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>%
  filter(cluster=="B cell")

score_list=list(gene=top100_B$gene)
stRNA <-AddModuleScore(stRNA,features = score_list,name = "B_cell") 
all.markers = sc.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top100_Pla = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>%
  filter(cluster=="Plasma cell")
score_list=list(gene=top100_Pla$gene)
stRNA <-AddModuleScore(stRNA,features = score_list,name = "Plasma_cell") 
all.markers = sc.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top100_EC = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>%
  filter(cluster=="Endothelial cell")

score_list=list(gene=top100_EC$gene)
stRNA <-AddModuleScore(stRNA,features = score_list,name = "Endothelial_cell")  
stRNA_df <- as.data.frame(stRNA@meta.data)
average_scores <- stRNA_df %>%
  group_by(SCT_snn_res.0.4) %>%
  summarise(
    Osteoblastic_cells1_avg = mean(Osteoblast1, na.rm = TRUE),
    Myeloid_cells1_avg = mean(Myeloid_cell1, na.rm = TRUE),
    T_cells1_avg = mean(T_cell1, na.rm = TRUE),
    #Fibroblasts1_avg = mean(Fibroblasts1, na.rm = TRUE),
    B_cells1_avg = mean(B_cell1, na.rm = TRUE),
    #TNK_cells1_avg = mean(TNK_cells1, na.rm = TRUE),
    Plasma_cells1_avg = mean(Plasma_cell1, na.rm = TRUE),
    Endothelial_cells1_avg = mean(Endothelial_cell1, na.rm = TRUE)
  )
save(stRNA,file = "Data/MLN07_ST_0.4_addmo.RData")
library(ggplot2)
metadata <-stRNA@meta.data
my_colors <- c("0" = "#f33b99",
               "1" = "#E7481B",
               "2" = "#1b79af",
               "3" =  "#179b73",
               "4" =  "#8da0cb",
               "5" = "#bbe173",
               "6" = "#1b79af",
               "7" = "#bbe173")

p <- ggplot(metadata, aes(x = as.factor(seurat_clusters), y = Osteoblast1, fill = seurat_clusters)) +
  geom_violin(trim = FALSE) +
  labs(x = "seurat_clusters_ST", y = "Osteoblastic_cells_SC_score") +
  theme_minimal(base_size = 14) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")+
  scale_fill_manual(values = my_colors)
p
library(ggplot2)
ggsave("Cluster/Osteoblastic_cells_ST.pdf", plot = p, width = 7, height = 4, limitsize = FALSE)
p <- ggplot(metadata, aes(x = as.factor(seurat_clusters), y = Myeloid_cell1, fill = as.factor(seurat_clusters))) +
  geom_violin(trim = FALSE) +
  labs(x = "seurat_clusters_ST", y = "Myeloid_cells_SC_score") +
  theme_minimal(base_size = 14) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
p
ggsave("Cluster/Myeloid_cells1_ST.pdf", plot = p, width = 6, height = 4, limitsize = FALSE)

###################
p <- ggplot(metadata, aes(x = as.factor(seurat_clusters), y = Endothelial_cell1, fill = seurat_clusters)) +
  geom_violin(trim = FALSE) +
  labs(x = "seurat_clusters_ST", y = "Endothelial_cell_SC_score") +
  theme_minimal(base_size = 14) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")+
  scale_fill_manual(values = my_colors)
p

ggsave("Cluster/Endothelial_cell1_ST..pdf", plot = p, width = 7, height = 4, limitsize = FALSE)

########################################

p <- ggplot(metadata, aes(x = as.factor(seurat_clusters), y = T_cell1, fill = as.factor(seurat_clusters))) +
  geom_violin(trim = FALSE) +
  labs(x = "seurat_clusters_ST", y = "T_cell_SC_score") +
  theme_minimal(base_size = 14) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
p
dev.off()
ggsave("Cluster/T_cell_ST.pdf", plot = p, width = 6, height = 4, limitsize = FALSE)

###############################
p <- ggplot(metadata, aes(x = as.factor(seurat_clusters), y = B_cell1, fill = as.factor(seurat_clusters))) +
  geom_violin(trim = FALSE) +
  labs(x = "seurat_clusters_ST", y = "B_cell1_SC_score") +
  theme_minimal(base_size = 14) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
p
ggsave("Cluster/B_cell_ST.pdf", plot = p, width = 6, height = 4, limitsize = FALSE)

#################################################
p <- ggplot(metadata, aes(x = as.factor(seurat_clusters), y = Plasma_cell1, fill = as.factor(seurat_clusters))) +
  geom_violin(trim = FALSE) +
  labs(x = "seurat_clusters_ST", y = "Plasma_cells_SC_score") +
  theme_minimal(base_size = 14) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
p
ggsave("Cluster/Plasma_cells_ST.pdf", plot = p, width = 6, height = 4, limitsize = FALSE)


all.markers = st.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top50 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster=="0")

score_list=list(gene=top50$gene)

scRNA <-AddModuleScore(scRNA,features = score_list,name = "ST_cluster0")


all.markers = st.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top50 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster=="1")

score_list=list(gene=top50$gene)

scRNA <-AddModuleScore(scRNA,features = score_list,name = "ST_cluster1")

all.markers = st.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top50 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster=="2")

score_list=list(gene=top50$gene)

scRNA <-AddModuleScore(scRNA,features = score_list,name = "ST_cluster2")

all.markers = st.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top50 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster=="3")

score_list=list(gene=top50$gene)

scRNA <-AddModuleScore(scRNA,features = score_list,name = "ST_cluster3")

all.markers = st.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top50 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster=="4")

score_list=list(gene=top50$gene)

scRNA <-AddModuleScore(scRNA,features = score_list,name = "ST_cluster4")

all.markers = st.marker %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top50 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  filter(cluster=="5")

score_list=list(gene=top50$gene)

scRNA <-AddModuleScore(scRNA,features = score_list,name = "ST_cluster5")

scRNA_df <- as.data.frame(scRNA@meta.data)
average_scores_SC <- scRNA_df %>%
  group_by(celltype) %>%
  summarise(
    ST_cluster01_avg = mean(ST_cluster01, na.rm = TRUE),
    ST_cluster11_avg = mean(ST_cluster11, na.rm = TRUE),
    ST_cluster21_avg = mean(ST_cluster21, na.rm = TRUE),
    ST_cluster31_avg = mean(ST_cluster31, na.rm = TRUE),
    ST_cluster41_avg = mean(ST_cluster41, na.rm = TRUE)
  )

cell.type <- c("Cluster0","Cluster1","Cluster2","Osteoblast","Endothelial cell","Cluster5")
Idents(MLN07_ST) <- "seurat_clusters"
names(cell.type) <- levels(MLN07_ST)
MLN07_ST <- RenameIdents(MLN07_ST, cell.type)
MLN07_ST$celltype <- Idents(MLN07_ST)
table(MLN07_ST$celltype)
table(MLN07_ST$seurat_clusters)

colors <- c("#f33b99","#E7481B","#1b79af","#179b73","#8da0cb","#bbe173","#970030","#3d4a78","#f6e36d","#ad5f2c")
p1 <- DimPlot(MLN07_ST, reduction = "umap", label = TRUE,cols = colors)
p1
ggsave("Cluster/cluster_0.4_.pdf", p1, width = 12, height = 8) 

p2 <- SpatialDimPlot(MLN07_ST, label = TRUE, label.size = 3, group.by = "celltype", 
                     cols = c("Cluster0" = "#f33b99",
                              "Cluster1" = "#E7481B",
                              "Cluster2" = "#1b79af",
                              "Osteoblast" =  "#179b73",
                              "Endothelial cell" =  "#8da0cb",
                              "Cluster5" = "#bbe173"))
p2
ggsave("Cluster/louvain_0.4.pdf", p2, width = 12, height = 8) 
ggsave("Cluster/louvain_0.4.png", p2, width = 12, height = 8)
p <- p1 + p2
ggsave("Cluster/cluster_louvain_0.4.pdf", p, width = 15, height = 7)  
ggsave("Cluster/cluster_louvain_0.4.png", p, width = 15, height = 7)

p1<-SpatialFeaturePlot(MLN07_ST, features = "MCAM")
p1
ggsave("Cluster/MCAM_ST.png", plot = p1, width = 10, height = 8, dpi = 300)
ggsave("Cluster/MCAM_ST.pdf", plot = p1, width = 10, height = 8, dpi = 1000)
