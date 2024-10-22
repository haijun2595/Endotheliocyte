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
	CellphoneDB 5.0.0
	survival 3.2-10
 	Ubuntu 18.04
	R 4.0.5 and 4.1.0 
        CellChat 1.6.1
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
### 1. Download data
"1_SRRraw data download.txt" is used to download the single-cell SRR raw file and decompress it into a "fastq.gz" file, and finally execute the SARTsolo command to get the processed file.
	 
"2_FASTQ raw data download.txt" is used to download the single-cell "fastq.gz" and execute the SARTsolo command to get the processed file.
### 2. Quality control and subgroup
"1_pancancer_data_python_process.txt" is used to import the processed data into the python environment and perform quality control to eliminate low quality cells and genes and finally perform dimensionality reduction clustering on the data.

"2_pan-cancer huge group analysis.txt" was used to define overall cellular subpopulations and plot subpopulation umap plots, cell scale plots, and other plots

### 3. NKT cell
"1_NKT cell extraction code (marker).txt" was used to extract NKT cell subclusters.
"2_T_celL_DEGs.txt" was used to plot the umap of T-cell subpopulations, cell scale plots, and other graphs.

"3_AddModule.txt" was used to score T cell subpopulations (e.g., TEFF and TEX scores).

"4_python_DNT_marker_heatmap.txt" was used to draw the DNT subcluster Marker heat map.

"5_CD4+CD8+ DNT correlation analysis.txt" was used for correlation analysis between T cell subclusters.

"6_DNT cell verification.txt" was used to reference external data to verify DNT presence.

"7_Heat map of cell ratio and RO.txt" was used to draw heat maps of DNT in different tissue sources with RO enrichment analysis.

"8_palantir.txt" was used for pseudotiming analysis of T cell subclusters.

"9_machinelearning.txt" was used to construct an immunotherapy prognostic model for the T-cell subcluster combined with the TIGER database.


### 5. Myeloid cell
"1_Extracellular matrix scoring.txt" was used to score extracellular matrix remodeling for subclusters of cells between samples.

"2_M1M2score.txt" was used to assess the M1/M2 polarization phenotype of myeloid cell subpopulations by scores.

"3_scFEA.txt" was used for the analysis of cellular metabolic reprogramming.

"4_Prognostic analysis.txt" was used for prognostic analysis of cell subpopulations in the joint the TCGA database.

### 6. B cell
"1_SCENIC.txt" was used to analyze the expression levels of transcription factors in cellular subclusters among different sample groups.

"2_monocle3.txt" was used to characterize the developmental trajectory of B cells.

"3_Gene_Trajectory.txt" was used to analyze gene trajectory activity in B cells.

"4_cellphoneDB.txt" was used to analyze ligand-receptor pairs between different cell types.

### 7. Generic analysis
"1_OR.txt" was used to show the dominance ratio of cell subpopulation distribution in each tissue.

"2_GSVA.txt" was used to explore enrichment pathways in cellular subclusters.
### 8. Spatial transcriptome
"TESLA.txt" was used to verify the presence of specific cells in the spatial transcriptome.
## Contact
If you have any questions or feedback, feel free to contact me at tanghaijun@sr.gxmu.edu.cn.
