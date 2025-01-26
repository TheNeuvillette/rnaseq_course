# ----------------------------------------
#
# Part 6: Differential expression analysis
#
# ----------------------------------------

# Install and load necessary packages:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
BiocManager::install('org.Mm.eg.db')

library(DESeq2)
library(EnhancedVolcano)
library(org.Mm.eg.db)

# Set directory path:
setwd("~/UniBern/rnaseq_course/06_Differential_expression_analysis")



# ------------------------------------------------------------------------
#
# The following part is only necessary to run if the environment has been 
# reset after part 5 (Exploratory data analysis).
# 
# It consists of:
#   1. Creating the DESeqDataSet object (DESeqDataSetFromMatrix())
#   2. Pre-filtering rows with low gene counts (10 or less counts per row).
#   3. The differential expression analysis (DESeq2::DESeq())
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

# -----------------------------------------------------------------
#
# From this point on, the differential expression analysis begins.
#
# -----------------------------------------------------------------

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



# Convert the Ensembl ID to the gene symbol:
res_lung.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res_lung.df), keytype = 'ENSEMBL', column = 'SYMBOL')
res_blood.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res_blood.df), keytype = 'ENSEMBL', column = 'SYMBOL')
head(res_lung.df); nrow(res_lung.df)           # 29369
head(res_blood.df); nrow(res_blood.df)         # 30062

# Remove NA from results:
res_lung.df <- na.omit(res_lung.df)
res_blood.df <- na.omit(res_blood.df)
head(res_lung.df); nrow(res_lung.df)           # 19551
head(res_blood.df); nrow(res_blood.df)         # 19712



# Determine the Cutoff for log2FoldChange and padj:
pCutoff = 1e-05     # padj cutoff
FCcutoff = 2        # log2FoldChange cutoff



# Make the Volcano Plots:
selectLab_lung <- c('Wars1', 'Gbp2', 'Gbp5', 'Tap1', 'Igtp', 'Upp1', 'Cd274', 
               'Irgm1', 'Serpina3g','Zbp1', 'Acod1', 'Ms4a3',
               'Cd79a', 'Fmod', 'Krt13', 'Krt4', 'Col15a1', 'Ly6d')

EnhancedVolcano(res_lung.df, x = 'log2FoldChange', y = 'padj',
                selectLab = selectLab_lung,
                lab = res_lung.df$symbol,
                pCutoff = pCutoff,
                FCcutoff = FCcutoff,
                title = "Volcano Plot of DE genes - Lung Case vs. Lung Control",
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

selectLab_blood <- c('Gbp2','Cd79a', 'Gbp5', 'Plin2', 'Upp1', 'Niban3', 'Batf2',
                     'Fpr2', 'Ermn', 'N4bp1', 'Cd79a', 'Dcaf12', 'Spib',
                     'Fcer2a', 'Cd79b', 'Sftpc', 'Ptprn2', 'Sox2', 'Ctnna2')

EnhancedVolcano(res_blood.df, x = 'log2FoldChange', y = 'padj',
                selectLab = selectLab_blood,
                lab = res_blood.df$symbol,
                pCutoff = pCutoff,
                FCcutoff = FCcutoff,
                title = "Volcano Plot of DE genes - Blood Case vs. Blood Control",
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.75)



# Count how many genes are differentially expressed (DE) with padj < pCutoff:
# Total genes DE with padj < pCutoff:
sum(res_lung$padj <= pCutoff)                  # 3365
sum(res_blood$padj <= pCutoff)                 # 6726

# Genes DE with padj < pCutoff and have a gene symbol:
sum(res_lung.df$padj <= pCutoff)               # 3065
sum(res_blood.df$padj <= pCutoff)              # 5986



# Count DE genes are up-regulated vs down-regulated:
# Total upregulated DE genes with log2FoldChange > FCcutoff:
sum(res_lung$log2FoldChange >= FCcutoff)      # 1682
sum(res_blood$log2FoldChange >= FCcutoff)     # 2621

# Total downregulated DE genes with log2FoldChange < -FCcutoff:
sum(res_lung$log2FoldChange <= -FCcutoff)     # 2538
sum(res_blood$log2FoldChange <= -FCcutoff)    # 6635

# Upregulated DE genes with log2FoldChange > FCcutoff and have a gene symbol:
sum(res_lung.df$log2FoldChange >= FCcutoff)   # 851
sum(res_blood.df$log2FoldChange >= FCcutoff)  # 1200

# Downregulated DE genes with log2FoldChange < -FCcutoff and have a gene symbol:
sum(res_lung.df$log2FoldChange <= -FCcutoff)  # 1655
sum(res_blood.df$log2FoldChange <= -FCcutoff) # 4264



# Count DE genes with padj < pCutoff that are up-regulated vs down-regulated:
# Total upregulated DE genes:
sum(res_lung$log2FoldChange >= FCcutoff & res_lung$padj <= pCutoff)             # 678
sum(res_blood$log2FoldChange >= FCcutoff & res_blood$padj <= pCutoff)           # 705

# Total downregulated DE genes: 
sum(res_lung$log2FoldChange <= -FCcutoff & res_lung$padj <= pCutoff)            # 728
sum(res_blood$log2FoldChange <= -FCcutoff & res_blood$padj <= pCutoff)          # 3083

# Upregulated DE genes and have a gene symbol:
sum(res_lung.df$log2FoldChange >= FCcutoff & res_lung.df$padj <= pCutoff)       # 553
sum(res_blood.df$log2FoldChange >= FCcutoff & res_blood.df$padj <= pCutoff)     # 604

# Downregulated DE genes and have a gene symbol:
sum(res_lung.df$log2FoldChange <= -FCcutoff & res_lung.df$padj <= pCutoff)      # 646
sum(res_blood.df$log2FoldChange <= -FCcutoff & res_blood.df$padj <= pCutoff)    # 2584



# Based on the original publication, select 2-3 genes that are of particular
# interest and investigate their expression level.

# Zbp1: https://www.uniprot.org/uniprotkb/Q9QY24/entry
res_lung.df[res_lung.df$symbol == 'Zbp1', ]
res_blood.df[res_blood.df$symbol == 'Zbp1', ]
Zbp1 <- plotCounts(dds, gene='ENSMUSG00000027514', intgroup="Group", returnData = TRUE)
boxplot(count ~ Group , data=Zbp1, main = "Expression of Zbp1")


# Gbp5: https://www.uniprot.org/uniprotkb/Q8CFB4/entry
res_lung.df[res_lung.df$symbol == 'Gbp5', ]
res_blood.df[res_blood.df$symbol == 'Gbp5', ]
Gbp5 <- plotCounts(dds, gene='ENSMUSG00000105504', intgroup="Group", returnData = TRUE)
boxplot(count ~ Group , data=Gbp5, main = "Expression of Zbp1")