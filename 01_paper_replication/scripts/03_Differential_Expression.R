# ===========================================================================
# Gene Expression Analysis of Colorectal Cancer Progression
# Replication of: Kim SK et al. (2014), "A nineteen gene-based
#                 risk score classifier predicts prognosis of
#                 colorectal cancer patients", Molecular Oncology
# DOI:10.1016/j.molonc.2014.06.016 | PMID: 25049118
# GitHub repo: https://github.com/mikikeric/crc_gene_expression_replication.git
# Script: 03_Differential_Expression.R - running DESeq2, finding significant
#                                        genes and visualising results
#
# Note: this script uses objects created in 02_Load_Data.R (sample_info,
#       counts_matrix) - run that script first in same session.
# ===========================================================================

# Confirm working directory
getwd()

library(ggplot2)
library(pheatmap)
library(DESeq2)
library(org.Hs.eg.db)

# Optional step - saves a copy of everything that shows up in my console into
# a text file, so I have a record of what happened.
library(TeachingDemos)
txtStart("01_paper_replication/output/Rsave_03_Differential_Expression.txt")

# Setting "normal" as the reference group to compare everything else against,
# since that's the natural baseline instead of it being DESeq's default choice
# which is alphabetical)
sample_info$biopsy_site = relevel(factor(sample_info$biopsy_site),
                                  ref = "normal")

# Building the DESeq2 dataset by telling it my counts, sample info and that I
# want to compare groups based on biopsy_site (normal vs primary tumour vs
# liver metastasis).
#
# DESeqDataSetFromMatrix() builds special object DESeq2 needs to work with
# design = ~ biopsy_site tells DESeq2 to group samples by biopsy_site and
# compare gene expression between those groups.
dds = DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_info, 
  design = ~ biopsy_site)

# Performing the actual differential expression analysis
# DESeq() does all the stastitical work - this function runs whole comparison
dds = DESeq(dds)

# Checking what comparisons are available noe
# This lists the specific group comparisons DESeq2 has calculated 
# (e.g. tumour vs normal)
resultsNames(dds)

# Extracting the actual results table for primary tumour vs normal comparison
# results() extracts a specific comparison's stats (fold-change, p-value)
# from dds
tumour_vs_normal = results(dds,
                           name = "biopsy_site_primary.tumor_vs_normal")

# Taking a look at first few rows
head(tumour_vs_normal)

# Summary of the results - how many gene went up, down, how many were
# significant
summary(tumour_vs_normal)

# Turning results into a normal table (data frame) so filtering and
# sorting is more easy
#
# as.data.frame() converts DESeq2's special results object into ordinary table
tumour_results_df = as.data.frame(tumour_vs_normal)

# Filtering down to genes that are actually significant - using a stricter
# threshold closer to what the paper used (p < 0.001)
# Also requiring at least a 2-fold difference in expression (matches paper's
# criteria)
#
# subset() filters rows based on a condition
# pvalue < 0.001 matches the paper's significance threshold
# abs(log2FoldChange) >= 1 means at least a 2-fold change (up or down) since
# abs() ignores the +/- direction
significant_genes = subset(tumour_results_df,
                           pvalue < 0.001 & abs(log2FoldChange) >= 1)

# Sorting so the most significant genes appear first
# order() sorts by a column
# smallest p-value (most significant) first
significant_genes = significant_genes[order(significant_genes$pvalue), ]

# How many genes made the cut?
# nrow() counts how many rows (genes) are in this filtered table
nrow(significant_genes)

# Taking a look at top ones
head(significant_genes)

# Converting Ensembl gene IDs into readable gene names so I can recognise the
# genes instead of seeing ID codes.
#
# library(org.Hs.eg.db) being a lookup dictionary which maps gene IDs to
# gene names (human genes)
if (!require("org.Hs.eg.db"))
  BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# mapIds() looks up each ID and returns its readable name
#
# keys = my list of Ensembl IDs to look up
# column = what I want back - the gene name
# keytype = the type of ID I'm giving it
# multiVals = if more than one match, just take the first
significant_genes$gene_symbol = mapIds(org.Hs.eg.db,
                                       keys = rownames(significant_genes),
                                       column = "SYMBOL",
                                       keytype = "ENSEMBL",
                                       multiVals = "first")

# Taking another look, now with gene names included
head(significant_genes)

# Saving this table as a proper output file in my repo
# write.csv() saves a table as a .csv file - simple format anyone can
# open (e.g. in Excel)
write.csv(significant_genes, 
          "01_paper_replication/output/significant_genes_tumour_vs_normal.csv")

# Checking specifically whether TREM1 and CTGF (two key regulators paper
# identified) show up as significant in my own results
# %in% checks if each gene_symbol matches either "TREM1" or "CTGF"
significant_genes[significant_genes$gene_symbol %in% c("TREM1", "CTGF"), ]

# Volcano plot: shows every tested gene, plotting fold-change (x-axis) against
# statistical significance (y-axis).
# This creates a TRUE/FALSE column marking which genes meet my significance
# threshold
tumour_results_df$significant = tumour_results_df$pvalue < 0.001 &
  abs(tumour_results_df$log2FoldChange) >= 1

# ggplot() builds the plot
# geom_point() draw one dot per gene
# -log10(pvalue) flips small p-values into tall bars, so more significant genes
# appear higher up
ggplot(tumour_results_df, aes(x = log2FoldChange,
                              y = -log10(pvalue),
                              color = significant)) +
  geom_point(alpha = 0.5) +
  labs(title = "Primary Tumour vs Normal - Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10(pvalue)") + 
  theme_minimal()

# Saving the volcano plot as an image file in my output folder
# ggsave() saves the most recently created plot
# width/height set image size in inches
ggsave("01_paper_replication/output/volcano_plot_tumour_vs_normal.png",
       width = 7, height = 6)

# Heatmap: shows expression level of my top genes across all samples, so I can
# visually see the groups (normal/tumour) clustering apart
#
# Picking my top 20 most significant genes to keep the heatmap readable
# order() sorts by adjusted p-value (padj) with most significant first
# head(..., 20) keeps just the top 20 rows
top_genes = head(significant_genes[order(significant_genes$padj), ], 20)

# Getting normalised expression values for just these genes, across all samples
# counts() pulls out expression data from dds; normalized = TRUE adjusts for
# technical differences between samples so they are fairly comparable
normalised_counts = counts(dds, normalized = TRUE)

heatmap_data = normalised_counts[rownames(top_genes), ]
# Keeping the rows (genes) that match my top 20 genes list only

# Labelling rows with actual gene names instead of Ensembl IDs
rownames(heatmap_data) = top_genes$gene_symbol

# Drawing the heatmap
pheatmap (heatmap_data,
          scale = "row",
          show_colnames = FALSE,
          main = "Top 20 Genes: Primary Tumour vs Normal")
# scale = "row" makes each gene's colours relative to itself, so patterns are
# easier to see across genes with very different expression levels
#
# show_colnames = FALSE hides sample ID labels, sine there are 54 of them and
# they would be unreadable

# Saving heatmap as an image file in my output folder
png("01_paper_replication/output/heatmap_top20_tumour_vs_normal.png",
    width = 800,
    height = 600)
pheatmap(heatmap_data,
         scale = "row",
         show_colnames = FALSE,
         main = "Top 20 Genes: Primary Tumour vs Normal")
dev.off()
# png() opens a file connection to save into
# dev.off() closes it and finishes writing the file
# pheatmap doesn't work with ggsave since it's not a ggpplot object

# Stopping the console log for this script
txtStop()
