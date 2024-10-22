#############
library(Seurat)
data <- read.table("TARGETcount(1)2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1) 
patients <- read.table("TARGET_Clinical2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)   
patient_ids <- patients[[1]] 

transposed_data <- t(data)
transposed_df <- as.data.frame(transposed_data)
specific_patient_data <- transposed_df[rownames(transposed_df) %in% patient_ids, ]
print(specific_patient_data)
scRNA <- readRDS("Endo_merge_after.rds")
Idents(scRNA) <- "subcelltype"
markers_MCAM <- FindMarkers(scRNA, ident.1 = "Tip_like_ECs", only.pos = TRUE)
significant_markers_MCAM <- markers_MCAM[markers_MCAM$p_val < 0.05, ]
sorted_genes <- significant_markers_MCAM[order(-significant_markers_MCAM$avg_log2FC), ]
print(head(sorted_genes))
top_1471_genes <- rownames(head(sorted_genes,30))###############
library(GSVA)
gene_list <- list(top1471=top_1471_genes)
expr_matrix <- as.matrix(specific_patient_data)
expr_matrix_transposed <- t(expr_matrix)
ssgsea_scores <- gsva(expr_matrix_transposed, gene_list, method = "ssgsea")
print(ssgsea_scores)
df_ssgsea <- as.data.frame(t(ssgsea_scores))
df_ssgsea$ID <- rownames(df_ssgsea)
merged_data <- merge(df_ssgsea, patients[, c("ID", "Effect")], 
                     by = "ID", all.x = TRUE)
head(merged_data)
rownames(merged_data) <- merged_data$ID
merged_data$ID <- NULL  
colnames(merged_data)[colnames(merged_data) == "top1471"] <- "score"
head(merged_data)
merged_data$Effect <- ifelse(merged_data$Effect %in% c("0"), "poor", "good")
library(ggplot2)

######
shapiro_test_poor <- shapiro.test(merged_data$score[merged_data$Effect == "poor"])
shapiro_test_good <- shapiro.test(merged_data$score[merged_data$Effect == "good"])
cat("Shapiro test p-value for poor:", shapiro_test_poor$p.value, "\n") 
cat("Shapiro test p-value for good:", shapiro_test_good$p.value, "\n")

#######
var_test <- var.test(score ~ Effect, data = merged_data)
cat("Var test p-value:", var_test$p.value, "\n")

#######
t_test_result <- t.test(score ~ Effect, data = merged_data, var.equal = T)  
cat("T-test p-value:", t_test_result$p.value, "\n")

#######
wilcox_test_result <- wilcox.test(score ~ Effect, data = merged_data)
cat("Mann-Whitney U test p-value:", wilcox_test_result$p.value, "\n")

plot <- ggplot(merged_data, aes(x = Effect, y = score, fill = Effect)) +
  geom_boxplot(width = 0.7) +
  geom_text(aes(label = sprintf("p = %.3f", t_test_result$p.value), y = max(score) * 0.95, x = 1.5),
            size = 5, fontface = "plain") + 
  scale_fill_manual(values = c("pink", "#EFC000FF")) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_blank(), 
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"), 
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12) 
  )
print(plot)
ggsave("plot.pdf", plot, width = 8, height = 5)

################################################
library(survival)
library(survminer)
inputFile="survival.txt"   ####TARGET
rt=read.table(inputFile,header=TRUE,sep="\t",check.names=FALSE)
rt=rt[,c("futime","fustat","MCAM")]
colnames(rt)=c("futime","fustat","var")
res.cut=surv_cutpoint(rt,time="futime",event="fustat",variables=c("var"))
res.cut
res.cat=surv_categorize(res.cut)

fit=survfit(Surv(futime,fustat)~var,data=res.cat)

diff=survdiff(Surv(futime,fustat)~var,data=res.cat)
pValue=1-pchisq(diff$chisq,df=1) 

if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.3f",pValue))
}
surPlot = ggsurvplot(fit,
                     data = res.cat,
                     conf.int = TRUE,  
                     pval = pValue, 
                     legend.title = "MCAM", 
                     legend.labs = c("High", "Low"), 
                     xlab = "Time (years)",
                     ylab = "Survival Probability",  
                     palette = c("#FF7F0E", "#65c2a4"), 
                     risk.table = TRUE,  
                     risk.table.title = "Number at risk (MCAM)", 
                     risk.table.col = "strata",  
                     risk.table.height = 0.25, 
                     risk.table.y.text.col = TRUE,  
                     risk.table.y.text = TRUE,  
                     break.time.by = 1) 

pdf("survival_MCAM_TARGET.pdf",width = 6,height = 6)
print(surPlot)
dev.off()



