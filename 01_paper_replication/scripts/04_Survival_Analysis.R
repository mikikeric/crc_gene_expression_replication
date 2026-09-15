# ===========================================================================
# Gene Expression Analysis of Colorectal Cancer Progression
# Replication of: Kim SK et al. (2014), "A nineteen gene-based
#                 risk score classifier predicts prognosis of
#                 colorectal cancer patients", Molecular Oncology
# DOI:10.1016/j.molonc.2014.06.016 | PMID: 25049118
# GitHub repo: https://github.com/mikikeric/crc_gene_expression_replication.git
# Script: 04_Survival_Analysis.R - checking for survival data and running a
#                                  simplified survival analysis
#
# Note: this script uses objects created in 02_Load_Data.R and
#       03_Differential_Expression.R (sample_info, significant_genes, dds)
#       - run those scripts first in same session
# ===========================================================================

# Confirm working directory
getwd()

library(ggplot2)
library(survival)
library(survminer)
library(DESeq2)

# Optional step - saves a copy of everything that shows up in my console into
# a text file, so I have a record of what happened.
library(TeachingDemos)
txtStart("01_paper_replication/output/Rsave_04_Survival_Analysis.txt")

# Checking what columns exist in my sample info to see if there's any patient
# survival/outcome data I can use for my survival analysis step
names(sample_info)

# Result: no survival or follow-up columns are present in this dataset -
# only organism, individual, organism_part, disease, disease_staging and
# biopsy_site.
#
# This is consistent with the paper's methodology. The AMC cohort (this
# dataset) was used exclusively for the RNA-seq differential expression
# analysis comparing normal, primary tumour and metastatic tissue. The paper's 
# survival/prognosis analysis was developed and validated using 4 separate 
# external cohorts (CIT, AUS, UAPC and HI), comprising over 1000 additional 
# patients - a scale of validation beyond the scope of this replication project.
#
# Conclusion: a Kaplan-Meier survival analysis cannot be performed on this
# dataset as it was not designed to support one

txtStop()
