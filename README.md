# Bulk-RNAseq-Workflow

At CIG-Core, we perform data analysis for Bulk RNA-seq data on regular basis. As a result we needed a workflow that
allowed us to repeatedly analyze different datasets with minimum changes and new inputs. Hence, we developed the following
workflow that preprocess fastq files to get abundance and counts data followed by downstream analysis. Here, we use the
following tools:

1 FastQC - To check quality of Fastq files
2. MultiQC - To prepare QC reports for quality report of raw, trimmed and aligned data
3. Trimmomatic - To trim sequenced data as and when needed (For trimming adapter sequences, over-represented sequences etc.)
4. Kallisto - Alignment (pseudoalignment with bootstraps)
5. Sleuth - Read in the counts and abundance data and perform DE analysis at transcript level.

The folder comprises of two main R Markdown scripts that can be used to "preprocess.Rmd" and "downstream_analysis.Rmd".
The "preprocess.Rmd" calls on different inhouse R functions that creates SLURM scripts and submits them to HPC cluster to 
perform preprocessing steps like QC, trimming and alignment and generates reports.
The "downstream_analysis.Rmd" performs DE Analysis using Sleuth followed by Over Representation Analysis (ORA) and Gene Set 
Enrichment Analysis (GSEA) using clusterProfiler. An html report is generated upon knitting the file and also various images, 
excel files and csv files are generated and saveed at respective output folder.
