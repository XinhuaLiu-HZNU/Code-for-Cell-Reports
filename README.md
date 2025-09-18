# Code-for-Cell-Reports
Code for RNA-seq, scRNA-seq and CUT&amp;Tag analysis reside in my work in Cell Reports
Integrated NGS Analysis Pipelines
This repository contains automated pipelines for processing and analyzing next-generation sequencing (NGS) data from three major modalities: CUT&Tag, bulk RNA-seq, and single-cell RNA-seq (scRNA-seq). All pipelines are configured for the mouse genome (mm10).

Overview
Pipeline	Description	Key Inputs	Key Outputs
CUT&Tag	Processes paired-end sequencing data for histone modification or transcription factor binding analysis.	Raw FastQ files	BAM files, BED fragments, BigWig coverage tracks
Bulk RNA-seq	Processes paired-end RNA-seq data for transcriptome-wide gene expression quantification.	Raw FastQ files	BAM files, gene count tables
scRNA-seq	Analyzes single-cell RNA-seq data from CellRanger output for clustering, annotation, and differential expression.	CellRanger filtered_feature_bc_matrix directories	Annotated Seurat object, cluster markers, differential expression results
1. CUT&Tag Analysis Pipeline
Description
This pipeline processes paired-end CUT&Tag data, performing quality control, alignment, filtering, and generation of genome-wide coverage tracks.

Key Steps:
Quality Control & Adapter Trimming: fastp

Alignment: bowtie2 (very-sensitive mode) to mm10.

BAM Processing: Sorting, indexing, and filtering for MAPQ ≥ 5.

Fragment Analysis: Extracts valid paired-end fragments (≤ 2000 bp).

Coverage Track Generation: Creates normalized BigWig files for visualization.

Usage:
Create a text file (fqfiles.txt) listing your sample names (one per line).

Ensure reference paths (e.g., mm10index, mm10.chrsize) are correct.

Run the main bash script: bash CUT&Tag_Mouse.sh

2. Bulk RNA-seq Analysis Pipeline
Description
This pipeline processes standard bulk RNA-seq data for gene expression quantification using the STAR aligner.

Key Steps:
Quality Control & Trimming: fastp

Alignment & Quantification: STAR to mm10.

BAM Processing: Convert SAM to BAM, sort by name for HTSeq.

Gene Counting: htseq-count using an NCBI RefSeq GTF file.

Usage:
Create a text file (fqfiles.txt) listing your sample names.

Ensure the paths to the STAR genome index and GTF file are correct.

Run the main bash script: bash RNAseqAna_PE_Mouse.sh

3. Single-Cell RNA-seq Analysis Pipeline
Description
This R-based pipeline performs end-to-end analysis of scRNA-seq data, starting from a count matrix (e.g., CellRanger output) through clustering, annotation, and differential expression.

Key Steps:
Data Loading & Integration: Reads multiple samples and integrates them using Seurat.

Quality Control: Filters low-quality cells based on gene counts and mitochondrial percentage.

Clustering: Performs PCA, UMAP, and graph-based clustering. Optimal resolution is determined using clustree.

Cell Annotation: Annotates cell types using SingleR with the HumanPrimaryCellAtlasData reference.

Differential Expression & Pathway Analysis: Finds markers for clusters and performs GSEA using fgsea on KEGG pathways.

Usage:
Set the scdir variable to the path containing your CellRanger output folders.

The script is designed to handle multiple samples. Sample names are automatically parsed from the directory structure.

Run the R script step-by-step in an R environment with all required packages installed. Critical parameters (like the clustering resolution integrated_snn_res.1.4) should be adjusted based on the clustree output for your dataset.

Required R Packages:
r
Seurat, SeuratObject, dplyr, reshape2, ggplot2, patchwork, SingleR, clustree, celldex, pheatmap, clusterProfiler, EnsDb.Hsapiens.v86, fgsea, msigdbr
Installation & Dependencies
Software Requirements:
FastP: For read trimming.

Bowtie2: For CUT&Tag alignment.

STAR: For RNA-seq alignment.

Samtools: For BAM file processing.

Bedtools: For BED file operations.

deepTools (bamCoverage): For BigWig generation.

HTSeq: For gene counting (bulk RNA-seq).

R (v4.0+) with the packages listed above.

Reference Genomes:
All pipelines are configured for the mouse mm10 genome. Ensure you have the correct indices and annotation files:

Bowtie2 index

STAR index

Chromosome sizes file (mm10.chrsize)

GTF annotation file (e.g., mm10.ncbiRefSeq.gtf)

License
This project is licensed under the MIT License.

Citation
If you use this pipeline in your research, please consider citing this repository and the relevant primary software (Seurat, STAR, Bowtie2, etc.).

For detailed instructions and parameters, please refer to the comments within each script.
