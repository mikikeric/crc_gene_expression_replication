# ===========================================================================
# Gene Expression Analysis of Colorectal Cancer Progression
# Replication of: Kim SK et al. (2014), "A nineteen gene-based
#                 risk score classifier predicts prognosis of
#                 colorectal cancer patients", Molecular Oncology
# DOI:10.1016/j.molonc.2014.06.016 | PMID: 25049118
# GitHun repo: https://github.com/mikikeric/crc_gene_expression_replication.git
# Script: 01_Setup.R - project setup and package installation
# ===========================================================================

# Confirm working directory
getwd()

# Optional step - saves a copy of everything that shows up in my console into a
# text file, so I have backup record of what happened when I ran this.
library(TeachingDemos)
txtStart("01_paper_replication/output/Rsave_01_Setup.txt")

# Installing the packages I'll need later on:
# - ggplot2: for making plots like my volcano plot
# - pheatmap: for making my heatmap
# - survival & survminer: for the survival analysis later
install.packages(c("ggplot2", "pheatmap", "survival", "survminer"))

# DESeq2 checks gene by gene whether activity differs between the groups
# (normal vs tumour vs metastasis), or if it's just noise.
#
# It comes from Bioconductor, not R's normal package store, so I need
# BiocManager first to install it.
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("DESeq2")

# Loading everything now, just need to double check they all are actually
# installed properly and there's no errors.
library(ggplot2)
library(pheatmap)
library(survival)
library(survminer)
library(DESeq2)

# Stop logging the console - closes off the log file
txtStop()
