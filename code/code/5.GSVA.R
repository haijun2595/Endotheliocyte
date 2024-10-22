#####################

library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(GSEA)
library(dplyr)
library(limma)
library(Seurat)
library(pheatmap)
library(stringr)
library(ggplot2)
library(psych)
library(Matrix)
library(ggpubr)
library(SingleR)
rm(list=ls())

setwd("path")
table(Endo_merge_after$subcelltype)
table(Endo_merge_after$subcelltype)
new.cluster.ids <- c("Venous_ECs", "Capillary_ECs", "LUM_ECs", "Tip_like_ECs", "PRRX1_ECs", "Lymphatic_ECs", "Proli_ECs")
names(new.cluster.ids) <- levels(Endo_merge_after)
Endo_merge_after <- RenameIdents(Endo_merge_after, new.cluster.ids)
Endo_merge_after@meta.data$subcelltype <- Endo_merge_after@active.ident
Idents(Endo_merge_after) <- "orig.ident"
Endo_merge_after@meta.data$group = ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                            c("lungmeta10","lungmeta17"),"lung",
                                          ifelse(Endo_merge_after@meta.data$orig.ident %in%
                                                   c("NT01","NT03","NT07","NT08"),"NT","PT"))

Endo_merge_after$group <- as.factor(Endo_merge_after$group)  
table(Endo_merge_after$group)
PDX_analysis<-Endo_merge_after
PDX_analysis@meta.data$seurat_clusters <- as.factor(PDX_analysis@meta.data$seurat_clusters);######
str(PDX_analysis@meta.data$seurat_clusters);
Idents(PDX_analysis)="subcelltype"
###################GSVA
gsvatestdata <- PDX_analysis
View(gsvatestdata@meta.data)
Idents(gsvatestdata) <- "subcelltype"
table(gsvatestdata$subcelltype)
Average <- AverageExpression(gsvatestdata, assays = NULL, features = NULL, return.seurat = FALSE,
                             add.ident = NULL, slot = "data", use.scale = FALSE, use.counts = FALSE, verbose = TRUE)

data <- as.matrix(Average$RNA)
####################################################
genesets2 <- read.gmt("h.all.v2023.1.Hs.symbols.gmt")
genesets2 = split(genesets2$gene, genesets2$term)
GSVAresult <- gsva(data, genesets2, min.sz=10, max.sz=Inf, tau=1, method="gsva", kcdf="Poisson",
                   mx.diff=TRUE, abs.ranking=FALSE, verbose=TRUE, parallel.sz=10)
write.csv(GSVAresult, file = "OB_hall.csv")
t <- t(scale(t(GSVAresult)))
write.csv(t, file = "OB_hall_t.csv")
getwd()
###
t1 <- read.csv("OB_hall_t_.csv", header = T, row.names = 1)
t1 <- as.matrix(t1)
p <- pheatmap::pheatmap(t1, 
                        cluster_rows = F, 
                        cluster_cols = F,
                        show_colnames = T,
                        fontsize_row = 8, 
                        fontsize_col = 8, 
                        border_color = "white",
                        color = colorRampPalette(c("#4876B7", "#E9F6FA", "#EA6A3B"))(30),
                        angle_col = 45
)
ggsave("Hall.pdf",p,width = 6,height = 6)

########################################

library(Seurat)
library(GSVA)
library(GSEABase)
library(dplyr)
library(limma)
Idents(Endo_merge_after) <- "group"
before_after<-subset(Endo_merge_after,idents = c("PT"))
Idents(before_after) <- "orig.ident" 
before_after$group = ifelse(before_after$orig.ident %in%
                              c("PT01","PT02","PT03","PT04","PT05","PT06"),"before","after")
Idents(before_after) <- "group" 
before_after$group <- as.factor(before_after@active.ident)  
table(before_after$group)
table(before_after$subcelltype)

scRNA <- before_after
table(scRNA$subcelltype)
Idents(scRNA)="subcelltype"
Tip_like_ECs=subset(scRNA,ident="Tip_like_ECs")

table(Tip_like_ECs$subcelltype)
Tip_like_ECs$subcelltype <- Tip_like_ECs@active.ident
table(Tip_like_ECs$group)

#######################
###########
sce2 <- Tip_like_ECs
Idents(sce2) <- "Anno" 
expr2 <- as.matrix(sce2@assays$RNA@data)
GeneSet <- getGmt("h.all.v2023.1.Hs.symbols.gmt")

gsva.kegg2 <- gsva(expr2, 
                   gset.idx.list = GeneSet, 
                   kcdf="Gaussian",
                   method = "gsva",
                   parallel.sz=20)

group_list2 <- sce2@meta.data[,c("group")]
design <- model.matrix(~0+factor(group_list2))
colnames(design)=levels(factor(group_list2))
rownames(design)=colnames(expr)

contrast.matrix<-makeContrasts(contrasts = "after-before",levels = design)
fit <- lmFit(gsva.kegg2,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff2 <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(diff2)
diff2$group <- ifelse(diff2$logFC > 0 & diff2$P.Value < 0.05 ,"up" ,
                      ifelse(diff2$logFC < 0 & diff2$P.Value < 0.05 ,"down","noSig")
)

head(diff2)
diff2 <- diff2 %>% 
  mutate(hjust2 = ifelse(t>0,1,0)) %>% 
  mutate(nudge_y = ifelse(t>0,-0.1,0.1)) %>% 
  filter(group != "noSig") %>% 
  arrange(t) %>% 
  rownames_to_column("ID")
diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))
interest_pathways <- c( "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                        "HALLMARK_KRAS_SIGNALING_DN",
                        "HALLMARK_ANGIOGENESIS", 
                        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                        "HALLMARK_MARK_OXIDATIVE_PHOSPHORYLATION",
                        "HALLMARK_APICAL_SURFACE",
                        "HALLMARK_HEDGEHOG_SIGNALING",
                        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                        "HALLMARK_APOPTOSIS",
                        "HALLMARK_APICAL_JUNCTION",
                        "HALLMARK_NOTCH_SIGNALING",
                        "HALLMARK_HYPOXIA",
                        "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                        "HALLMARK_MYC_TARGETS_V1"
)
diff2_interest <- diff2 %>% filter(ID %in% interest_pathways)
p <- ggplot(diff2_interest, aes(ID, t, fill = group)) + 
  geom_bar(stat = 'identity', alpha = 0.7) + 
  scale_fill_manual(breaks = c("down", "up"), 
                    values = c("#008020", "#08519C")) +
  geom_text(data = diff2_interest, aes(label = diff2_interest$ID, 
                                       y = diff2_interest$nudge_y),
            nudge_x = 0, nudge_y = 0, hjust = diff2_interest$hjust,
            size = 3) + 
  labs(x = "HALLMARKERS pathways", 
       y = paste0("t value of GSVA score"),
       title = "Tip_like_ECs Post-chemotherapy VS Pre-chemotherapy") +
  scale_y_continuous(limits = c(-limt, limt)) +
  coord_flip() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5, size = 12),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right")

print(p)
ggsave(p, file="GSVA2.pdf", width = 7, height = 5)



