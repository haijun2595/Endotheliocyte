# Endotheliocyte
The ***endotheliocyte*** used for building **Single-cell and spatial transcriptomics reveals the key role of MCAM+ tip-like endothelial cells in osteosarcoma metastasis**.

## Table of Contents
1. System requirements
2. Installation guide
3. Codes
5. Contact

## System requirements
The following are the version numbers of the software or algorithms used in this study.

AUCell 1.12.0

Seurat 4.3.0.1

CellphoneDB 3.0.0

survival 3.2-10

Ubuntu 18.04

R 4.0.5 and 4.1.0 

CytoTRACE (v0.3.3)  

Monocle2 (v2.22.0)  

GSVA(v1.42.0 )      

GSEA(v1.2.0)      

Seurat(v4.3.0.1)  

CellChat (1.6.1)   

ClusterGVis (v0.1.1)    

survival (v3.5-5)    

survminer (v0.4.9)   

pRRophetic (v0.5)

## Installation guide
Python libraries can be installed in a shell environment using the "pip install" command. 

	pip install "library_name"

R packages can be installed in the R environment using the "install.packege()" or "BiocManager::install()" commands.

	install.packege("packege_name")
	if(!"BiocManager"%in%installed.packages()){ 
	install.packages("BiocManager")}
 	if(!"packege_name"%in%installed.packages()){ 
	BiocManager::install("packege_name")}
	if (!"devtools" %in% installed.packages()) {
  	install.packages("devtools")}
   	devtools::install_github("packege_name")
## Codes
Specific descriptions of the codes can be found in the corresponding documents.


"1.Quality control and subgroup.R" is used for quality control and cell cluster identification in single-cell data.

"2.DEG.R" primarily performs differential gene expression analysis.

"3.Score.R" is used to score functional gene sets within cell subpopulations.

"4.bulk.R" analyzes bulk data.

"5.GSVA.R" conducts functional enrichment pathway analysis for cell subpopulations.

"6.Interaction.R" is mainly applied to analyze cell-cell communication between subpopulations.

"7.Pseudo time.R" analyzes the differentiation trajectories of cell subpopulations.

"8.ST.R" is primarily used for spatial transcriptomics analysis of single-cell samples.

## Contact
If you have any questions or feedback, feel free to contact me at tanghaijun@sr.gxmu.edu.cn.
