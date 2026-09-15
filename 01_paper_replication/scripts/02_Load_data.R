# ===========================================================================
# Gene Expression Analysis of Colorectal Cancer Progression
# Replication of: Kim SK et al. (2014), "A nineteen gene-based
#                 risk score classifier predicts prognosis of
#                 colorectal cancer patients", Molecular Oncology
# DOI:10.1016/j.molonc.2014.06.016 | PMID: 25049118
# GitHub repo: https://github.com/mikikeric/crc_gene_expression_replication.git
# Script: 02_Load_Data.R - downloading and checking the dataset
# ===========================================================================

# Confirm working directory
getwd()

# Loading the packages I had already installed in 01_Setup.R, restarted my R
# session and need to load them again before I can use them.
library(ggplot2)
library(pheatmap)
library(survival)
library(survminer)
library(DESeq2)

# Optional step - saves copy of everything that shows up in my console into
# a text file, so I have a record of what happened.
library(TeachingDemos)
txtStart("01_paper_replication/output/Rsave_02_Load_Data.txt")

# This package lets me pull the dataset straight from Expression Atlas (a
# database that hosts already-processed public data), instead of downlaoding
# files from a website manually.
if (!require("ExpressionAtlas"))
  BiocManager::install("ExpressionAtlas")
library(ExpressionAtlas)

# Downloading the actual dataset from the paper.
# GSE50760 is this dataset's ID number in the GEO database - confirmed in the
# paper's Methods section.
#
# Dataset is 54 RNA-seq samples (normal colon, primary tumour,
# liver metastasis) from 18 colorectal cancer patients.
#
# Paper cites this dataset as GSE50760 (GEO format) but ExpressionAtlas needs
# the ArrayExpress format ID instead - same dataset, different ID system.
# So changed to E-GEOD-50760 to fix this.
crc_data = getAtlasData("E-GEOD-50760")

# Quick look at what got downloaded, just to see its structure.
crc_data

# crc-data has one item inside it, "list of length 1".
# Pulling that single item out so I can work with it directly instead of
# having to dig into the list every time.
crc_experiment = crc_data[["E-GEOD-50760"]]

# Looking at what's inside this now
crc_experiment

# There was one more item to obtain within it called "rnaseq", which should be
# the actual data object I need
crc_rnaseq = crc_experiment[["rnaseq"]]

# Looking inside this now
crc_rnaseq

# Extracting the sample information table - this tells me which sample belongs
# to which group (normal/tumour/metastasis)
# colData() pulls out the sample info
# as.data.frame() turns it into a normal table
sample_info = as.data.frame(colData(crc_rnaseq))

# Looking to see what the groups actually look like
# head() shows just the first few rows, so I'm not printing all 54
head(sample_info)

# Checking specifically what categories exist in the biopsy_site column, since
# this is likely what tells me normal vs tumour vs metastasis
# table() counts how many samples fall into each category
# the $ makes sure it gives you only the biopsy_site column as a list of values.
table(sample_info$biopsy_site)

# Extracting actual counts table (gene expression numbers) from
# inside crc_rnaseq.
# assay() pulls out a specific data table stored inside the object, in this
# case the one named "counts"
counts_matrix = assay(crc_rnaseq, "counts")

# Making sure the sample names in my counts table match the sample names in
# my sample info table in the same order - matters because DESeq2 assumes they
# line up correctly.
#
# colnames() = column names of the counts table (sample IDs)
# rownames() = row names of the sample info table (also sample IDs)
# == is for checking if each pair matches
# all() checks that every pair matches
all(colnames(counts_matrix) == rownames(sample_info))

# Stopping the console log for this script
txtStop()
