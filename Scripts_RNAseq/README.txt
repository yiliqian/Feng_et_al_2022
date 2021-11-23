Commands and python scripts for processing RNA sequencing data from raw fastq files to feature counts, RPKM, and statistical calculations

Description:

Each file contains general sample code to process a single RNA sequencing fastq file. "Commands.txt" is the only file provided that is not original code

- Commands.txt sample commands, listed in order of execution, for using FastQC tool, packages from BBTools suite including BBDuk, BBSplit, BBMap and a package from the SubRead suite named featureCounts.

- computeRPKM.py, inputs .count files and aligns GeneID with "product name" from a reference file and computes read per kilobase million (RPKM) for each geneID

- CountsAndRPKM.py, inputs .rpkm files from previous script and compiles the counts and RPKM values by geneID for every sample into a single .csv for later use with DESeq

- DESeq.R, R script using DESeq2 Bioconductor library to normalize read counts across samples and quantify differential gene expression between samples and control
