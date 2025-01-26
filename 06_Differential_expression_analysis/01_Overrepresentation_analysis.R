# ------------------------------------
#
# Part 7: Overrepresentation analysis
#
# ------------------------------------

# Install and load necessary packages:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

BiocManager::install("DESeq2")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("clusterProfiler")

library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# Set directory path:
setwd("~/UniBern/rnaseq_course/06_Differential_expression_analysis")

# ------------------------------------------------------------------------
#
# The following part is only necessary to run if the environment has been 
# reset after Part 6 (Differential expression analysis).
# 
# It consists of:
#   1. Creating the DESeqDataSet object (DESeqDataSetFromMatrix())
#   2. Pre-filtering rows with low gene counts (10 or less counts per row).
#   3. The differential expression analysis (DESeq2::DESeq())
#   4. Extraction of the results (DESeq2::results)
# 
# ------------------------------------------------------------------------

# Read in the reformatted counts data:
counts_data <- read.delim('../05_Exploratory_data_analysis/01_Reformat_Results/reformatted_reads_counts.tsv')

# Change the column names:
sample_names <- c("GeneID",
                  "SRR7821918",
                  "SRR7821919",
                  "SRR7821920",
                  "SRR7821921",
                  "SRR7821922",
                  "SRR7821937",
                  "SRR7821938",
                  "SRR7821939",
                  "SRR7821949",
                  "SRR7821950",
                  "SRR7821951",
                  "SRR7821952",
                  "SRR7821953",
                  "SRR7821968",
                  "SRR7821969",
                  "SRR7821970")

colnames(counts_data) <- sample_names

# Change the df such that the Geneid is now the row names:
rownames(counts_data) <- counts_data$GeneID
counts_data$GeneID <- NULL



# Read in the sample info:
coldata <- read.delim('../05_Exploratory_data_analysis/00_Other/Metadata.tsv')
coldata <- coldata[order(coldata$Sample), ]

# Change the df such that the sample name is now the row names:
rownames(coldata) <- coldata$Sample
coldata$Sample <- NULL



# Check that the sample names in 'counts_data' are the same as the sample names in 'coldata'.
# Check that the sample names are in the same order.
all(colnames(counts_data) %in% rownames(coldata))
all (colnames(counts_data) == rownames(coldata))


# Create the DESeqDataSet object using "DESeq2::DESeqDataSetFromMatrix()"
# Design = ~ condition
# To begin, we should do a 2 factor design. 

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_data,
                                      colData = coldata,
                                      design = ~ Group)
dds                                            # 78298



# Pre-filtering: Remove the rows with low gene counts:
# We set the threshold at 10 counts per row.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds                                            # 35743



# Run the DESeq2::DESeq() function and save to variable dds.
dds <- DESeq2::DESeq(dds)

# ------------------------------------------------------------

# Extract and checking results from the DESeq analyses (Lung and Blood)
res_lung <- DESeq2::results(dds, contrast = c("Group","Lung_WT_Case","Lung_WT_Control"))
res_blood <- DESeq2::results(dds, contrast = c("Group","Blood_WT_Case","Blood_WT_Control"))
head(res_lung); nrow(res_lung)                 # 35743
head(res_blood); nrow(res_blood)               # 35743

# Remove NA from results:
res_lung <- na.omit(res_lung)
res_blood <- na.omit(res_blood)
head(res_lung); nrow(res_lung)                 # 29369
head(res_blood); nrow(res_blood)               # 30062

# Convert the DESeqResults into a df:
res_lung.df <- as.data.frame(res_lung)
res_blood.df <- as.data.frame(res_blood)
head(res_lung.df); nrow(res_lung.df)           # 29369
head(res_blood.df); nrow(res_blood.df)         # 30062

# -----------------------------------------------------------------
#
# From this point on, the overrepresentation analysis analysis begins.
#
# -----------------------------------------------------------------

# Determine the Cutoff for log2FoldChange and padj:
pCutoff = 1e-05     # padj cutoff
FCcutoff = 2        # log2FoldChange cutoff

sigs_lung <- res_lung[(res_lung$log2FoldChange >= FCcutoff & res_lung$padj <= pCutoff) | 
                           (res_lung$log2FoldChange <= -FCcutoff & res_lung$padj <= pCutoff), ]
sigs_blood <- res_blood[(res_blood$log2FoldChange >= FCcutoff & res_blood$padj <= pCutoff |
                           res_blood$log2FoldChange <= -FCcutoff & res_blood$padj <= pCutoff), ]
head(sigs_lung); nrow(sigs_lung)               # 1406
head(sigs_blood); nrow(sigs_blood)             # 3788

# Convert the DESeqResults into a df:
sigs_lung.df <- as.data.frame(sigs_lung)
sigs_blood.df <- as.data.frame(sigs_blood)
head(sigs_lung.df); nrow(sigs_lung.df)         # 1406
head(sigs_blood.df); nrow(sigs_blood.df)       # 3788



# Overrepresentation analysis:
gene_lung <- rownames(sigs_lung.df)            # Ensembl IDs of all DE lung genes.
gene_blood <- rownames(sigs_blood.df)          # Ensembl IDs of all DE blood genes.
universe_lung <- rownames(res_lung.df)         # Ensembl IDs of all genes measured.
universe_blood <- rownames(res_blood.df)       # Ensembl IDs of all genes measured.



ego_lung <- enrichGO(
  gene          = gene_lung,
  universe      = universe_lung,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",         # Choose BP, MF, CC, or ALL
  keyType       = "ENSEMBL",
  pAdjustMethod = "BH",         # Benjamini-Hochberg adjustment
  qvalueCutoff  = pCutoff,      # Adjusted p-value threshold
  readable      = TRUE          # Convert Ensembl IDs to gene symbols
)

ego_blood <- enrichGO(
  gene          = gene_blood,
  universe      = universe_blood,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",         # Choose BP, MF, CC, or ALL
  keyType       = "ENSEMBL",
  pAdjustMethod = "BH",         # Benjamini-Hochberg adjustment
  qvalueCutoff  = pCutoff,      # Adjusted p-value threshold
  readable      = TRUE          # Convert Ensembl IDs to gene symbols
)

# Visualise the results using a dotplot:
dotplot(ego_lung, showCategory=10) + ggtitle("Top 10 Lung GO Terms by Biological Process")
dotplot(ego_blood, showCategory=10) + ggtitle("Top 10 Blood GO Terms by Biological Process")