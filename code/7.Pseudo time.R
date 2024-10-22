library(scater)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(patchwork)
scRNA <- readRDS("Endo_merge_after.rds")

table(scRNA$subcelltype)
new.cluster.ids <- c("Venous_ECs", "Capillary_ECs", "LUM_ECs", "Tip_like_ECs", "PRRX1_ECs", "Lymphatic_ECs", "Proli_ECs")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA@meta.data$subcelltype <- scRNA@active.ident

Idents(scRNA) <- "subcelltype" 
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
#####PCA
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
ElbowPlot(scRNA)
library("harmony")
scRNA <- scRNA %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

scRNA <- scRNA %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  RunUMAP(reduction="harmony", dims=1:15)
saveRDS(scRNA, file = "scRNA.rds")
Idents(scRNA) <- "orig.ident"
scRNA@meta.data$group = ifelse(scRNA@meta.data$orig.ident %in%
                                 c("lungmeta10","lungmeta17"),"lung",
                               ifelse(scRNA@meta.data$orig.ident %in%
                                        c("NT01","NT03","NT07","NT08"),"NT","PT"))

scRNA$group <- as.factor(scRNA$group)
table(scRNA$orig.ident)
table(scRNA$group)
table(scRNA$subcelltype)
library(slingshot)
library(tradeSeq)
library(BiocParallel)
library(RColorBrewer)
sce <- scRNA
sce <- as.SingleCellExperiment(scRNA, assay = "RNA") # SCT
#run slingshot
sce_slingshot1 <- slingshot(sce, 
                            reducedDim = 'UMAP',
                            clusterLabels = sce$subcelltype) #subcelltype
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_slingshot1$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce_slingshot1)$UMAP, col = plotcol, pch=16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(sce_slingshot1), lwd=2, col='black')

library(CytoTRACE)
library(ggplot2)
library(ggpubr)
table(is.na(scRNA$subcelltype))
phenotype = as.character(scRNA$subcelltype)
names(phenotype) <- rownames(scRNA@meta.data)
mat <- as.matrix(scRNA@assays$RNA@counts)
results <- CytoTRACE(mat = mat,ncores = 1)
emb <- scRNA@reductions[["umap"]]@cell.embeddings

plotCytoTRACE(results,phenotype = phenotype,
              emb = emb, 
              outputDir = "cytotrace")
slingsce<-SlingshotDataSet(sce_slingshot1)
pseudotimeED <- slingPseudotime(slingsce, na = FALSE)
cellWeightsED <- slingCurveWeights(slingsce)
counts<-sce_slingshot1@assays@data@listData$counts
###################
########################
# Run Monocle2
rm(list = ls())
rm(fData)
if(F){
  library(monocle)
  library(Seurat)
  library(dplyr)
  #load("Endo_merge_after.rds")
  View(Endo_merge_after@meta.data)
  table(Endo_merge_after$subcelltype)
  new.cluster.ids <- c("Venous_ECs", "Capillary_ECs", "LUM_ECs", "Tip_like_ECs", "PRRX1_ECs", "Lymphatic_ECs", "Proli_ECs")
  names(new.cluster.ids) <- levels(Endo_merge_after)
  Endo_merge_after <- RenameIdents(Endo_merge_after, new.cluster.ids)
  Endo_merge_after@meta.data$subcelltype <- Endo_merge_after@active.ident
  view(Endo_merge_after@meta.data)
  table(Endo_merge_after$subcelltype)
  view(Endo_merge_after@meta.data)  
    scRNAsub <- Endo_merge_after  
  ###Extract data, phenotype data, and feature data from the SeuratObject
  data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data =scRNAsub@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  ###Construct monocle cds
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  my_cds<-monocle_cds
  my_cds
  slotNames(my_cds)
  my_cds <- estimateSizeFactors(my_cds)
  my_cds <- estimateDispersions(my_cds)
  my_cds <- detectGenes(my_cds, min_expr = 0.1)
  print(head(fData(my_cds)))

  expressed_genes <- row.names(subset(fData(my_cds),
                                      num_cells_expressed >= 10))
  head(fData(my_cds))
  summary(fData(my_cds)$num_cells_expressed)
  #keep only genes expressed in greater than 5% of cells
  fData(my_cds)$use_for_ordering <- fData(my_cds)$num_cells_expressed > 0.05 * ncol(my_cds)
  table(fData(my_cds)$use_for_ordering)
  
  
  #PCA
  plot_pc_variance_explained(my_cds, return_all = F)
  my_cds <- reduceDimension(my_cds,max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',residualModelFormulaStr = "~orig.ident+ num_genes_expressed",verbose = TRUE)
  my_cds <- clusterCells(my_cds, verbose = F)
  head(pData(my_cds))
  plot_cell_clusters(my_cds, color_by = 'as.factor(Cluster)')
  clustering_DEG_genes <- differentialGeneTest(my_cds,fullModelFormulaStr = '~Cluster',cores = 6)
  dim(clustering_DEG_genes)
  clustering_DEG_genes %>% arrange(qval) %>% head()
  my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:400]
  my_cds <- setOrderingFilter(my_cds, ordering_genes = my_ordering_genes)
  my_cds <- reduceDimension(my_cds, method = 'DDRTree')
  my_cds <- orderCells(my_cds)
  saveRDS(my_cds,file = "my_cds.rds")
  
  Idents(Endo_merge_after) <- "orig.ident"
  Endo_merge_after@meta.data$group = ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                              c("lungmeta10","lungmeta17"),"lung",
                                            ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                                     c("NT01","NT03","NT07","NT08"),"NT","PT"))
  
  Endo_merge_after$group <- as.factor(Endo_merge_after$group)  ##
  Idents(Endo_merge_after) <- "group" 
  my_cds$group <- Endo_merge_after@active.ident
  
  table(my_cds$group)

  table(my_cds$subcelltype)
  my_cds$subcelltype <- as.factor(my_cds$subcelltype)

  new.cluster.ids <- c("Venous_ECs", "Capillary_ECs", "LUM_ECs", "Tip_like_ECs", "PRRX1_ECs", "Lymphatic_ECs", "Proli_ECs")
  
  levels(my_cds$subcelltype) <- new.cluster.ids
  table(my_cds$subcelltype)
  
  p <- plot_cell_trajectory(my_cds, color_by = "Pseudotime")
  p
  ggsave("Pseudotime.pdf",p,width = 8,height = 8)

  p <- plot_cell_trajectory(my_cds, color_by = "subcelltype") +scale_color_manual(values = c("#384998",
                                                                                             "#4988A2", 
                                                                                             "#A65697", 
                                                                                             "#E9BF1F",
                                                                                             "#2C702C",
                                                                                             "#DD9A63",
                                                                                             "#E67A8A"))+
    facet_wrap(~subcelltype, nrow = 3)
  p
  ggsave("subcelltype_single.pdf",p,width = 8,height = 8)
  
##########################
BEAM_res <- BEAM(my_cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)

top20<-BEAM_res[1:20,]

plot_genes_branched_heatmap(my_cds[row.names(subset(top20,
                                                    qval < 3.026364e-35)),], #
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 40,
                            use_gene_short_name = T,
                            show_rownames = T)

genes_set <- c("MCAM") 
p<-plot_genes_in_pseudotime(my_cds[genes_set, ], 
                            color_by = "subcelltype",     
                            nrow = 1,   
                            ncol = NULL) 
p <- p + 
  theme(legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12),  
        plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12),    
        axis.title = element_text(size = 14)) + 
  scale_color_manual(values = c("#384998",
                                "#4988A2", 
                                "#A65697", 
                                "#E9BF1F",
                                "#2C702C",
                                "#DD9A63",
                                "#E67A8A"))  

print(p)
ggsave(filename = "MCAM.pdf", plot = p, width = 6, height = 4)
