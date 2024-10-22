rm(list=ls())
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
setwd("/path")
#####Merging Data
PT01.data <- Read10X(data.dir = "path/PT01")
PT02.data <- Read10X(data.dir = "path/PT02")
PT03.data <- Read10X(data.dir = "path/PT03")
PT04.data <- Read10X(data.dir = "path/PT04")
PT05.data <- Read10X(data.dir = "path/PT05")
PT06.data <- Read10X(data.dir = "path/PT06")
PT07.data <- Read10X(data.dir = "path/PT07")
PT08.data <- Read10X(data.dir = "path/PT08")
NT01.data <- Read10X(data.dir = "path/NT01")
NT03.data <- Read10X(data.dir = "path/NT03")
NT07.data <- Read10X(data.dir = "path/NT07")
NT08.data <- Read10X(data.dir = "path/NT08")
lungmeta10.data <- Read10X(data.dir = "path/lungmeta10")
lungmeta17.data <- Read10X(data.dir = "path/lungmeta17")
BC2.data <- Read10X(data.dir = "path/BC2")
BC3.data <- Read10X(data.dir = "path/BC3")
BC5.data <- Read10X(data.dir = "path/BC5")
BC6.data <- Read10X(data.dir = "path/BC6")
BC16.data <- Read10X(data.dir = "path/BC16")


PT01 <- CreateSeuratObject(counts = PT01.data, project = "PT01", min.cells = 3, min.features = 200)
PT02 <- CreateSeuratObject(counts = PT02.data, project = "PT02", min.cells = 3, min.features = 200)
PT03 <- CreateSeuratObject(counts = PT03.data, project = "PT03", min.cells = 3, min.features = 200)
PT04 <- CreateSeuratObject(counts = PT04.data, project = "PT04", min.cells = 3, min.features = 200)
PT05 <- CreateSeuratObject(counts = PT05.data, project = "PT05", min.cells = 3, min.features = 200)
PT06 <- CreateSeuratObject(counts = PT06.data, project = "PT06", min.cells = 3, min.features = 200)
PT07 <- CreateSeuratObject(counts = PT07.data, project = "PT07", min.cells = 3, min.features = 200)
PT08 <- CreateSeuratObject(counts = PT08.data, project = "PT08", min.cells = 3, min.features = 200)
NT01 <- CreateSeuratObject(counts = NT01.data, project = "NT01", min.cells = 3, min.features = 200)
NT03 <- CreateSeuratObject(counts = NT03.data, project = "NT03", min.cells = 3, min.features = 200)
NT07 <- CreateSeuratObject(counts = NT07.data, project = "NT07", min.cells = 3, min.features = 200)
NT08 <- CreateSeuratObject(counts = NT08.data, project = "NT08", min.cells = 3, min.features = 200)
lungmeta10 <- CreateSeuratObject(counts = lungmeta10.data, project = "lungmeta10", min.cells = 3, min.features = 200)
lungmeta17 <- CreateSeuratObject(counts = lungmeta17.data, project = "lungmeta17", min.cells = 3, min.features = 200)
BC2 <- CreateSeuratObject(counts = BC2.data, project = "BC2", min.cells = 3, min.features = 200)
BC3 <- CreateSeuratObject(counts = BC3.data, project = "BC3", min.cells = 3, min.features = 200)
BC5 <- CreateSeuratObject(counts = BC5.data, project = "BC5", min.cells = 3, min.features = 200)
BC6 <- CreateSeuratObject(counts = BC6.data, project = "BC6", min.cells = 3, min.features = 200)
BC16 <- CreateSeuratObject(counts = BC16.data, project = "BC16", min.cells = 3, min.features = 200)


OS.combined <- merge(PT01,  y = list(PT02,PT03,PT04,PT05,PT06,PT07,PT08,NT01,NT03,NT07,NT08,lungmeta10,lungmeta17,BC2,BC3,BC5,BC6,BC16), add.cell.ids = c("PT01","PT02","PT03","PT04","PT05","PT06","PT07","PT08","NT01","NT03","NT07","NT08","lungmeta10","lungmeta17","BC2","BC3","BC5","BC6","BC16"), project = "Endothelial")
saveRDS(OS.combined, "OS.combined.rds")
view(OS.combined@meta.data)


#####Quality Control
#load("combined.Rdata")
table(OS.combined$orig.ident)
OS.combined[["percent.mt"]] <- PercentageFeatureSet(OS.combined, pattern = "^MT-")
pdf(file = "QC1.pdf",height = 8,width = 15)
VlnPlot(OS.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
plot1 <- FeatureScatter(OS.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE)
plot2 <- FeatureScatter(OS.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE)
plot1 + plot2
dev.off()

pdf(file="QC2.pdf",height = 8,width = 15)
VlnPlot(object = OS.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

OS.combined<- subset(OS.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)
OS.combined <- NormalizeData(OS.combined, normalization.method = "LogNormalize", scale.factor = 10000)
OS.combined <- FindVariableFeatures(OS.combined, selection.method = "vst", nfeatures = 2000)
dim(OS.combined)


pdf(file = "afterQC1.pdf")
p<-VlnPlot(OS.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
ggsave("afterQC1.pdf", p, width = 10, height = 6)


plot1 <- FeatureScatter(OS.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE)
plot2 <- FeatureScatter(OS.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE)
p<-plot1 + plot2
ggsave("afterQC2.pdf", p, width = 10, height = 6)



# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(OS.combined), 10)

# plot variable features with and without labels
pdf(file = "QC3.pdf")
plot1 <- VariableFeaturePlot(OS.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
p<-plot1 + plot2
ggsave("afterQC3.pdf", p, width = 10, height = 6)


all.genes <- rownames(OS.combined)
OS.combined <- ScaleData(OS.combined, features = all.genes)
OS.combined <- RunPCA(OS.combined, features = VariableFeatures(object = OS.combined))
VizDimLoadings(object = OS.combined, dims = 1:4, reduction = "pca",nfeatures = 20)
pdf(file = "PCA.pdf")
DimPlot(object = OS.combined, reduction = "pca",raster=FALSE)
dev.off()
print(OS.combined[["pca"]], dims = 1:5, nfeatures = 5)

pdf(file = "ElbowPlot.pdf")
ElbowPlot(OS.combined)
dev.off()


library("harmony")
OS.combined_harmony <- OS.combined %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
OS.combined_harmony <- OS.combined_harmony %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  RunUMAP(reduction="harmony", dims=1:15)
saveRDS(OS.combined_harmony, file = "OS.combined_harmony.rds")    

#######################################################
OS.combined_harmony <-FindClusters(OS.combined_harmony,resolution = 0.03)
view(OS.combined_harmony@meta.data)

pdf(file = "UMAP.pdf")    
DimPlot(OS.combined_harmony, label = T,raster=FALSE)
dev.off()
pdf(file = "UMAP2.pdf")    
DimPlot(OS.combined_harmony, group.by = "orig.ident",raster=FALSE)
dev.off()

pdf(file = "tSNE.pdf")
DimPlot(OS.combined_harmony, label = T, reduction = "tsne",raster=FALSE)
dev.off()
pdf(file = "tSNE2.pdf")    
DimPlot(OS.combined_harmony, group.by = "orig.ident", reduction = "tsne",raster=FALSE)
dev.off()

OS.combined.markers <- FindAllMarkers(OS.combined_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(OS.combined.markers, file = "OS_combined_markers.csv")


new.cluster.ids <- c("Osteoblast","Myeloid cell","Proliferative cell","T/NK cell","Pericyte","Osteoclast","Endothelial cell")
names(new.cluster.ids) <- levels(OS.combined_harmony)
OS.combined_harmony <- RenameIdents(OS.combined_harmony, new.cluster.ids)
OS.combined_harmony@meta.data$celltype <- OS.combined_harmony@active.ident


p <- DimPlot(OS.combined_harmony,reduction="tsne",label=T,raster=FALSE)
saveRDS(OS.combined_harmony, file = "data.rds")    

###########################
cell.group <-OS.combined_harmony
rm(OS.combined_harmony)
DimPlot(cell.group,raster=FALSE)
table(cell.group$orig.ident)
cell.group@meta.data$group = ifelse(cell.group@meta.data$orig.ident %in%
                                      c("lungmeta10","lungmeta17"),"lung",
                                    ifelse(cell.group@meta.data$orig.ident %in%
                                             c("NT01","NT03","NT07","NT08"),"NT","PT"))

Idents(cell.group) <- "celltype" 
DimPlot(cell.group,split.by = 'group',raster=FALSE)
view(cell.group@meta.data)
cell.group$group <- as.factor(cell.group$group)  
table(cell.group$group)
table(cell.group$celltype)

##############################################################
new_color <- c("lungmeta10"="#EFE2AA",
               "lungmeta17"="#83B4EF",
               "NT01" = "#FFA240" ,
               "NT03" ="#EED0E0", 
               "NT07"= "#EBAEA9",
               "NT08"= "#95A6DA",
               "PT01"= "#BFA6C9",
               "PT02"= "#F5E0BA",
               "PT03"= "#AED0DF",
               "PT04"= "#89B780",
               "PT05"= "#F5D8D0",
               "PT06"= "#CB95BB",
               "PT07"= "#AAD0AC",
               "PT08"= '#FF6F91',
               "BC2"=  "#00A287",
               "BC3"=  "#6899D3",
               "BC5"=  '#FF8066',
               "BC6"=  "#8ECFF8",
               "BC16"= '#C9B9A0')


group_color<- c("PT"="#66C2A5",
                "NT"="#F3AE63",
                "lung" ="#8DA0CB")

Idents(cell.group) <- "orig.ident"
table(cell.group$orig.ident)
cell.group$orig.ident <- cell.group@active.ident
new_order <- c("lungmeta10","lungmeta17","NT01", "NT03","NT07","NT08","PT01","PT02","PT03","PT04","PT05","PT06","PT07","PT08","BC2","BC3","BC5","BC6","BC16")
cell.group$orig.ident <- factor(cell.group$orig.ident, levels = new_order)



p <- DimPlot(cell.group, group.by = "orig.ident", reduction = "tsne",label = F, raster = FALSE, cols = new_color)
p
ggsave("tsne1.pdf", plot = p, width = 10, height = 8)

p <- DimPlot(cell.group, group.by = "group", reduction = "tsne",label = F, raster = FALSE, cols = group_color)
p
ggsave("tsne2.pdf", plot = p, width = 10, height = 8)

####################################
#############################
colors <- c("T/NK cell" = "#f57c6e",
            "Osteoblast" =  "#F2B56F",
            "Myeloid cell" = "#FAE69E",
            "Endothelial cell" = "#84C3B7",
            "Proliferative cell"="#88D8DB",
            "Pericyte"="#F2A7DA",
            "Osteoclast"="#BE7FBC")
Idents(cell.group) <- "celltype"

p <- DimPlot(cell.group, split.by = "group", reduction = "tsne",group.by = "celltype",label = F,raster=FALSE,cols = colors)
ggsave("tSNE4.pdf",plot = p,width = 14, height = 7)

###############################################

Idents(cell.group) <- "celltype"
table(cell.group$celltype)
markers <- list(
  Osteoblast=c("ALPL","RUNX2","CLEC11A","IBSP","CDH11"),
  Myeloid=c("CD14", "APOE","S100A9","FCGR3A","LYZ"),
  Proliferative=c("MKI67","CDK1","TOP2A","PCNA"),
  T_NK=c("NKG7","TRAC","CD3D","GZMA","GZMK"),
  Pericyte=c("ACTA2","RGS5","TAGLN","MYL9"),
  Osteoclast=c("MMP9","ACP5","CTSK"),
  Endothelial=c("PECAM1","EGFL7","VWF","CDH5","CD34","CLDN5"))

unique_markers <- unique(unlist(markers))
############################################
p <- DotPlot(cell.group, features = unique_markers) +
  scale_colour_gradientn(
    colours = c("#3D95AB", "#AABC7E", "#DCAF17", "#E13611"),
    values = scales::rescale(c(0, 1, 2, 3))  
  ) +
  coord_flip() + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 10),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  xlab("Markers") + 
  ylab("Cell Types") 

print(p)
ggsave("dot_gene.pdf", plot = p, width = 10, height = 10)

######################
pdf("Sample.pdf", width=8.2, height=6)

doughnut <- function (x, labels = names(x), edges = 200, outer.radius = 0.8,
                      inner.radius=0.6, clockwise = FALSE,
                      init.angle = if (clockwise) 90 else 0, density = NULL,
                      angle = 45, col = NULL, border = FALSE, lty = NULL,
                      main = NULL, ...)
{
  if (!is.numeric(x) || any(is.na(x) | x < 0))
    stop("'x' values must be positive.")
  if (is.null(labels))
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L])
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col))
    col <- if (is.null(density))
      palette()
  else par("fg")
  col <- rep(col, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  twopi <- if (clockwise)
    -2 * pi
  else 2 * pi
  t2xy <- function(t, radius) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p),
         y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),
              outer.radius)
    polygon(c(P$x, 0), c(P$y, 0), density = density[i],
            angle = angle[i], border = border[i],
            col = col[i], lty = lty[i])
    Pout <- t2xy(mean(x[i + 0:1]), outer.radius)
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * Pout$x, c(1, 1.05) * Pout$y)
      text(1.1 * Pout$x, 1.1 * Pout$y, labels[i],
           xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0),
           ...)
    }      
    Pin <- t2xy(seq.int(0, 1, length.out = n*nx),
                inner.radius)
    polygon(Pin$x, Pin$y, density = density[i],
            angle = angle[i], border = border[i],
            col = "white", lty = lty[i])
  }
  
  title(main = main, ...)
  invisible(NULL)
}

df <- table(cell.group$celltype) %>% as.data.frame()
labs <- paste0(df$Var1," (", round(df$Freq/sum(df$Freq)*100,2), "%)")

p.circle <- doughnut(
  df$Freq,
  labels=labs, 
  init.angle=50,    
  col = c("T/NK cell" = "#f57c6e",
          "Osteoblast" =  "#F2B56F",
          "Myeloid cell" = "#FAE69E",
          "Endothelial cell" = "#84C3B7",
          "Proliferative cell"="#88D8DB",
          "Pericyte"="#F2A7DA",
          "Osteoclast"="#BE7FBC") , 
  border="white",    
  inner.radius= 0.4,
  cex = 1.5)       

dev.off()
#######################

#############################
library(ggplot2)
scRNA <- cell.group
view(scRNA@meta.data)
table(scRNA$orig.ident)

cellnum <- table(scRNA$orig.ident, scRNA$celltype)

cell.prop<-as.data.frame(prop.table(cellnum,1)) 
colnames(cell.prop)<-c("Sample","Celltype","Proportion")  
new_order <- c("lungmeta10", "lungmeta17", "NT01", "NT03", "NT07", "NT08", "PT01", "PT02", "PT03", "PT04", "PT05", "PT06", "PT07", "PT08", "BC2", "BC3", "BC5", "BC6", "BC16")
cell.prop$Sample <- factor(cell.prop$Sample, levels = new_order)
###############
colors <- c("T/NK cell" = "#f57c6e",
            "Osteoblast" =  "#F2B56F",
            "Myeloid cell" = "#FAE69E",
            "Endothelial cell" = "#84C3B7",
            "Proliferative cell"="#88D8DB",
            "Pericyte"="#F2A7DA",
            "Osteoclast"="#BE7FBC")
ggplot(cell.prop, aes(x = Sample, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1), 
    axis.text.y = element_text(size = rel(1))  
  )
ggsave("Sample.pdf", width = 8, height = 6)

OS.combined_harmony<-cell.group
save(OS.combined_harmony, file = "data_defined.Rdata")   
#############


#####################################################
rm(OS.combined)
Endo <- subset(OS.combined_harmony,idents = c("Endothelial cells"))###
view(Endo@meta.data)
table(Endo$orig.ident)
#####PCA
all.genes <- rownames(Endo)
Endo <- ScaleData(Endo, features = all.genes)
Endo <- RunPCA(Endo, features = VariableFeatures(object = Endo))
VizDimLoadings(object = Endo, dims = 1:4, reduction = "pca",nfeatures = 20)
pdf(file = "PCA.pdf")
DimPlot(object = Endo, reduction = "pca")
dev.off()
print(Endo[["pca"]], dims = 1:5, nfeatures = 5)
pdf(file = "ElbowPlot.pdf")
ElbowPlot(Endo)
dev.off()
library("harmony")
Endo_harmony <- Endo %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
#load("Endo_harmony1.Rds")
view(Endo_harmony@meta.data)
Endo_harmony <- Endo_harmony %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  RunUMAP(reduction="harmony", dims=1:15)
saveRDS(Endo_harmony, file = "Endo_harmony.rds")

Endo<-Endo_harmony ###
##################
Endo <-FindClusters(Endo,resolution = 0.03)
new.cluster.ids <- c("S0","S1","S2","S3")
names(new.cluster.ids) <- levels(Endo)
Endo <- RenameIdents(Endo, new.cluster.ids)
Endo@meta.data$subcelltype <- Endo@active.ident
view(Endo@meta.data)
table(Endo$subcelltype)
pdf(file = "UMAP.pdf")    
DimPlot(Endo, label = T)
dev.off()
pdf(file = "UMAP2.pdf")    
DimPlot(Endo, group.by = "orig.ident")
dev.off()
pdf(file = "tSNE.pdf")
DimPlot(Endo, label = T, reduction = "tsne")
dev.off()
pdf(file = "tSNE2.pdf")    
DimPlot(Endo, group.by = "orig.ident", reduction = "tsne")
dev.off()

Endo.markers <- FindAllMarkers(Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Endo.markers, file = "Endo_markers.csv")
table(Endo$subcelltype)
new.cluster.ids <- c("S0","S1_Endo_PRRX1","S2_Endo_TFF3","S3_Endo_CENPF")
names(new.cluster.ids) <- levels(Endo)
Endo <- RenameIdents(Endo, new.cluster.ids)
Endo@meta.data$subcelltype <- Endo@active.ident
view(Endo@meta.data)
table(Endo$subcelltype)
view(Endo@meta.data) 
saveRDS(Endo, file = "Endo_define_first.rds")
table(Endo$subcelltype)
#####################################
S0_Endo <- subset(Endo,idents = c("S0"))#########
view(S0_Endo@meta.data)
table(S0_Endo$orig.ident)
Idents(S0_Endo) <- "sample" 
table(S0_Endo$orig.ident)

#####PCA
all.genes <- rownames(S0_Endo)
S0_Endo <- ScaleData(S0_Endo, features = all.genes)
S0_Endo <- RunPCA(S0_Endo, features = VariableFeatures(object = Endo))
VizDimLoadings(object = S0_Endo, dims = 1:4, reduction = "pca",nfeatures = 20)
pdf(file = "S0_EndoPCA.pdf")
DimPlot(object = S0_Endo, reduction = "pca")
dev.off()
print(S0_Endo[["pca"]], dims = 1:5, nfeatures = 5)

pdf(file = "S0_EndoElbowPlot.pdf")
ElbowPlot(S0_Endo)
dev.off()
library("harmony")
S0_Endo_harmony <- S0_Endo %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
view(S0_Endo_harmony@meta.data)

S0_Endo_harmony <- S0_Endo_harmony %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  RunUMAP(reduction="harmony", dims=1:15)
saveRDS(S0_Endo_harmony, file = "S0_Endo_harmony_1.rds")

#################
S0_group <-FindClusters(S0_Endo_harmony,resolution = 0.1)
new.cluster.ids <- c("S0-0","S0-1","S0-2","S0-3","S0-4")
names(new.cluster.ids) <- levels(S0_group)
S0_group <- RenameIdents(S0_group, new.cluster.ids)
S0_group@meta.data$subcelltype <- S0_group@active.ident
view(S0_group@meta.data)
table(S0_group$subcelltype)
saveRDS(S0_group, file = "S0_group_.rds")
S0_group.markers <- FindAllMarkers(S0_group, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(S0_group.markers, file = "S0_group.markers.csv")
new.cluster.ids <- c("S0_Endo_MCAM","S0_Endo_ACKR1","S0_Endo_LUM","S0_Endo_CA4","S0_Endo_SDC1")
names(new.cluster.ids) <- levels(S0_group)
S0_group <- RenameIdents(S0_group, new.cluster.ids)
S0_group@meta.data$subcelltype <- S0_group@active.ident
view(S0_group@meta.data)
table(S0_group$subcelltype)
view(S0_group@meta.data)  
saveRDS(S0_group, file = "S0_group_defined.rds")
table(S0_group$subcelltype)
#####remove the low quality cell
Idents(S0_group) <- "subcelltype"
S0_group_Dele <- subset(S0_group,idents = c("S0_Endo_MCAM","S0_Endo_ACKR1","S0_Endo_LUM","S0_Endo_CA4"))
saveRDS(S0_group_Dele, file = "S0_group_Dele.rds")
############Running Merge
load("S0_group_Dele.rds")

#S0_group <- S0_group_zhushi.rds
view(Endo@meta.data)
table(Endo$subcelltype)
table(S0_group$subcelltype)
Endo_S1 <- subset(Endo,idents = c("S1_Endo_PRRX1"))
Endo_S2 <- subset(Endo,idents = c("S2_Endo_TFF3"))
Endo_S3 <- subset(Endo,idents = c("S3_Endo_CENPF"))
view(S0_group)
Endo_merge_afte <- merge(S0_group_Dele, y = list(Endo_S1,Endo_S2,Endo_S3), add.cell.ids = c("S0_group_Dele","Endo_S1","Endo_S2","Endo_S3"), project = "Endothelial_merge")
saveRDS(Endo_merge_afte, "Endo_merge_after.rds")

################################################
table(Endo_merge$subcelltype)
Endo_2 <- Endo_merge_after
Idents(Endo_2) <- "subcelltype"
Endo_2$subcelltype <- Endo_2@active.ident
library(dplyr)
##################################
table(Endo_2$subcelltype)
new.cluster.ids <- c("Venous_ECs", "Capillary_ECs", "LUM_ECs", "Tip_like_ECs", "PRRX1_ECs", "Lymphatic_ECs", "Proli_ECs")
names(new.cluster.ids) <- levels(Endo_2)
Endo_2 <- RenameIdents(Endo_2, new.cluster.ids)
Endo_2@meta.data$subcelltype <- Endo_2@active.ident
view(Endo_2@meta.data)
table(Endo_2$subcelltype)
view(Endo_2@meta.data)  
##########################
library(Seurat)
Endo_2 <- NormalizeData(Endo_2)
Endo_2 <- FindVariableFeatures(Endo_2, selection.method = "vst", nfeatures = 2000)
all.genes <- VariableFeatures(Endo_2)
Endo_2 <- ScaleData(Endo_2, features = all.genes)
Endo_2 <- RunPCA(Endo_2, features = all.genes)
Endo_2 <- RunUMAP(Endo_2,reduction = "pca", dims = 1:20)
Endo_2 <- RunTSNE(Endo_2,reduction = "pca",dims = 1:20)
######################
table(Endo_2$orig.ident)
Endo_2@meta.data$group = ifelse(Endo_2@meta.data$orig.ident %in%
                                  c("lungmeta10","lungmeta17"),"lung",
                                ifelse(Endo_2@meta.data$orig.ident %in%
                                         c("NT01","NT03","NT07","NT08"),"NT","PT"))

Idents(Endo_2) <- "subcelltype"
table(Endo_2$subcelltype)
Endo_2$group <- as.factor(Endo_2$group)
table(Endo_2$group)
#########################################################################
colors <- c("Venous_ECs" = "#384998",
            "Capillary_ECs" = "#4988A2", 
            "LUM_ECs" = "#A65697", 
            "Tip_like_ECs" = "#E9BF1F",
            "PRRX1_ECs"="#2C702C",
            "Lymphatic_ECs"="#DD9A63",
            "Proli_ECs"="#E67A8A")

p <- DimPlot(Endo_2, label = F,reduction = "umap",group.by = "subcelltype",raster=FALSE,cols = colors)
p
ggsave("Merge_UMAP.pdf",plot = p,width = 10, height = 8)

p <- DimPlot(Endo_2, reduction = "umap", split.by = "group", group.by = "subcelltype",label = F,raster = FALSE,cols = colors)
ggsave("Merge_UMAP2.pdf",plot = p,width = 12, height = 8)

p <- DimPlot(Endo_2, label = FALSE, reduction = "tsne", group.by = "subcelltype", raster = FALSE, cols = colors) +
  theme(plot.title = element_text(size = 20),  
        legend.text = element_text(size = 18)) 
p <- p + labs(title = "Subcelltype")
print(p)
ggsave("Merge_tSNE.pdf",plot = p,width = 10, height = 8)

p <- DimPlot(Endo_2, reduction = "tsne",split.by = "group", group.by = "subcelltype",label = F,raster=FALSE,cols = colors)
ggsave("Merge_tSNE2.pdf",plot = p,width = 10, height = 8)
###########
###########
table(Endo_2$subcelltype)
markers <- list(
  Venous_ECs=c("ACKR1","SELP","VCAM1"),
  Capillary_ECs=c("CA4","CD36","RGCC"),
  LUM_ECs=c("SPP1","LUM"),
  Tip_like_ECs=c("CXCR4","ESM1","NID2","PGF","LXN","ANGPT2"),
  PRRX1_Endo=c("RGS5","ACTA2","PRRX1"),
  Lymphatic_ECs=c("LYVE1","PDPN","PROX1"),
  Proli_ECs=c("MKI67","CDK1","TOP2A","PCNA")
)

unique_markers <- unique(unlist(markers))
#################################3
p <- DotPlot(Endo_2, features = unique_markers) +
  scale_colour_gradientn(
    colours = c("#3D95AB", "#AABC7E", "#DCAF17", "#E13611"), 
    values = scales::rescale(c(0, 1, 2, 3)) 
  ) +
  coord_flip() + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 10),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  xlab("Markers") + 
  ylab("Cell Types") 

# Print the plot
print(p)
ggsave("EC_markers.pdf", plot = p, width = 10, height = 9)

###############
Endo_merge_after <- Endo_2
Endo_merge_after@meta.data$group = ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                  c("lungmeta10","lungmeta17"),"lung",
                                ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                         c("NT01","NT03","NT07","NT08"),"NT","PT"))
############################
library(ggplot2)
table(Endo_merge_after$group)
table(Endo_merge_after$subcelltype)
cellnum <- table(Endo_merge_after$group,Endo_merge_after$subcelltype)
cell.prop<-as.data.frame(prop.table(cellnum,1)) 
colnames(cell.prop)<-c("Group","Subcelltype","Proportion")  
new_order <- c("PT","NT","lung")
cell.prop$Group <- factor(cell.prop$Group,levels = new_order)

###############
colors <- c("Venous_ECs" = "#384998",
            "Capillary_ECs" = "#4988A2", 
            "LUM_ECs" = "#A65697", 
            "Tip_like_ECs" = "#E9BF1F",
            "PRRX1_ECs"="#2C702C",
            "Lymphatic_ECs"="#DD9A63",
            "Proli_ECs"="#E67A8A")


p <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Subcelltype, 
                           stratum = Subcelltype, alluvium = Subcelltype)) + 
  geom_col(position = 'stack', width = 0.8) +  
  geom_stratum(width = 0.8, color = 'white') +  
  scale_fill_manual(values = colors) + 
  theme_bw() + 
  theme(
    axis.text = element_text(colour = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid = element_blank() 
  scale_y_continuous(expand = c(0, 0)) + 
  coord_flip()  
p
ggsave("col.pdf", width = 6, height = 3)

print(p)

##########################

table(Endo_merge_after$subcelltype)

####
colors <- c("Venous_ECs" = "#384998",
            "Capillary_ECs" = "#4988A2", 
            "LUM_ECs" = "#A65697", 
            "Tip_like_ECs" = "#E9BF1F",
            "PRRX1_ECs"="#2C702C",
            "Lymphatic_ECs"="#DD9A63",
            "Proli_ECs"="#E67A8A")
Idents(Endo_merge_after) <- "subcelltype"

my_comparisons <- list(c("Tip_like_ECs", "PRRX1_ECs"),
                       c("Tip_like_ECs", "Proli_ECs"),
                       c("Tip_like_ECs", "Venous_ECs"))
################################

veous_data <- FetchData(Endo_merge_after, vars = "MCAM", cells = WhichCells(Endo_merge_after, idents = "Venous_ECs"))
capillary_data <- FetchData(Endo_merge_after, vars = "MCAM", cells = WhichCells(Endo_merge_after, idents = "Capillary_ECs"))
tip_like_data <- FetchData(Endo_merge_after, vars = "MCAM", cells = WhichCells(Endo_merge_after, idents = "Tip_like_ECs"))
prrx1_data <- FetchData(Endo_merge_after, vars = "MCAM", cells = WhichCells(Endo_merge_after, idents = "PRRX1_ECs"))
proli_data <- FetchData(Endo_merge_after, vars = "MCAM", cells = WhichCells(Endo_merge_after, idents = "Proli_ECs"))

 shapiro.test(venous_data$MCAM)    
 shapiro.test(capillary_data$MCAM) 
 shapiro.test(tip_like_data$MCAM)  
 shapiro.test(prrx1_data$MCAM)    
 shapiro.test(proli_data$MCAM)  


 my_comparisons <- list(c("Tip_like_ECs", "PRRX1_ECs"),
                        c("Tip_like_ECs", "Proli_ECs"),
                        c("Tip_like_ECs", "Venous_ECs"))

VlnPlot(Endo_merge_after, features = "MCAM", pt.size = 0,cols = colors) &
  stat_compare_means(method = "wilcox.test", hide.ns = F,         
                     comparisons = my_comparisons,              
                     label = "p.signif",                 
                     bracket.size = 0.8,             
                     tip.length = 0,               
                     size = 6) &
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave("MCAM1_wilcox.testt.pdf",width = 10, height = 8)

###############################
Idents(Endo_merge_after) <- "orig.ident"
Endo_merge_after@meta.data$group = ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                            c("lungmeta10","lungmeta17"),"lung",
                                          ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                                   c("NT01","NT03","NT07","NT08"),"NT","PT"))

table(Endo_merge_after$group)
Idents(Endo_merge_after) <- "group" #
Endo_before_after_DEseq <-subset(Endo_merge_after,idents = c("PT"))
table(Endo_before_after_DEseq$subcelltype)
table(Endo_before_after_DEseq$group)
table(Endo_before_after_DEseq$orig.ident)

Idents(Endo_before_after_DEseq) <- "orig.ident" 
Endo_before_after_DEseq$group = ifelse(Endo_before_after_DEseq$orig.ident %in%
                                         c("PT01","PT02","PT03","PT04","PT05","PT06"),"before","after")
Idents(Endo_before_after_DEseq) <- "group" 
Endo_before_after_DEseq$group <- as.factor(Endo_before_after_DEseq@active.ident) 

table(Endo_before_after_DEseq$group)
table(Endo_before_after_DEseq$subcelltype)

Idents(Endo_before_after_DEseq) <- "subcelltype" 
table(Endo_before_after_DEseq$orig.ident)
table(Endo_before_after_DEseq$group)
table(Endo_before_after_DEseq$subcelltype)
Endo_before_after_DEseq) <- "group" 
colors <- c("before" = "#179b73",
            "after" ="#d48aff")
Idents(Endo_before_after_DEseq) <- "group"

############
################################
before_data <- FetchData(Endo_before_after_DEseq, vars = "MCAM", cells = WhichCells(Endo_before_after_DEseq, idents = "before"))
after_data <- FetchData(Endo_before_after_DEseq, vars = "MCAM", cells = WhichCells(Endo_before_after_DEseq, idents = "after"))


 shapiro.test(before_data$MCAM)    
 shapiro.test(after_data$MCAM)
my_comparisons <- list(c("after", "before"))

#################################
###############
colors <- c("Venous_ECs" = "#384998",
            "Capillary_ECs" = "#4988A2", 
            "LUM_ECs" = "#A65697", 
            "Tip_like_ECs" = "#E9BF1F",
            "PRRX1_ECs"="#2C702C",
            "Lymphatic_ECs" = "#DD9A63",
            "Proli_ECs" = "#E67A8A")
table(Endo_before_after_DEseq$subcelltype)
table(Endo_before_after_DEseq$group)
Idents(Endo_before_after_DEseq) <- "subcelltype"

p <- DimPlot(Endo_before_after_DEseq, label = FALSE, reduction = "tsne", group.by = "subcelltype",
              raster = FALSE) +
  theme(plot.title = element_text(size = 23),      
        legend.text = element_text(size = 18),      
        strip.text = element_text(size = 25),      
        axis.line = element_line(size = 1.5),        
        axis.text = element_text(size = 12),        
        axis.title = element_text(size = 14))       
print(p)
dev.off()
ggsave("tSNE2.pdf",width = 12, height = 9)
##################################################

table(Endo$subcelltype)
table(Endo$group)

Idents(Endo) <- "group"
PT_group<-subset(Endo,idents = c("PT"))
Idents(PT_group) <- "orig.ident" 
PT_group$Secgroup = ifelse(PT_group$orig.ident %in%
                             c("PT02","PT03","PT04"),"Lung_Metastasis","No_Lung_Metastasis")
Idents(PT_group) <- "Secgroup" 
PT_group$Secgroup <- as.factor(PT_group@active.ident)  
table(PT_group$Secgroup)
table(PT_group$subcelltype)
table(PT_group$orig.ident)

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
library(ggalluvial)
library(reshape2)

table(PT_group$subcelltype)
#############
colors <- c("Venous_ECs" = "#384998",
            "Capillary_ECs" = "#4988A2", 
            "LUM_ECs" = "#A65697", 
            "Tip_like_ECs" = "#E9BF1F",
            "PRRX1_ECs"="#2C702C",
            "Lymphatic_ECs"="#DD9A63",
            "Proli_ECs"="#E67A8A")
library(ggplot2)
table(PT_group$Secgroup)
table(PT_group$subcelltype)
cellnum <- table(PT_group$Secgroup,PT_group$subcelltype)
cell.prop<-as.data.frame(prop.table(cellnum,1)) 

colnames(cell.prop)<-c("Group","Subcelltype","Proportion") 
#new_order <- c("PT","NT","lung")
#cell.prop$Group <- factor(cell.prop$Group,levels = new_order)

colors <- c("Venous_ECs" = "#384998",
            "Capillary_ECs" = "#4988A2", 
            "LUM_ECs" = "#A65697", 
            "Tip_like_ECs" = "#E9BF1F",
            "PRRX1_ECs"="#2C702C",
            "Lymphatic_ECs"="#DD9A63",
            "Proli_ECs"="#E67A8A")


p <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Subcelltype, 
                           stratum = Subcelltype, alluvium = Subcelltype)) + 
  geom_col(position = 'stack', width = 0.8) + 
  geom_stratum(width = 0.8, color = 'white') + 
  #geom_alluvium(alpha = 0.4, width = 0.8, color = 'white', 
  #             linewidth = 1, curve_type = "linear") +  
  scale_fill_manual(values = colors) +
  theme_bw() + 
  theme(
    axis.text = element_text(colour = 'black'), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid = element_blank()  
  ) + 
  xlab('') +  
  ylab('') + 
  scale_y_continuous(expand = c(0, 0)) + 
  coord_flip() 
p
ggsave("col2.pdf", width = 6, height = 2)








