# ---------------------------------
#
# Part 5: Exploratory data analysis
#
# ---------------------------------

# Install and load necessary packages:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

# Set directory path:
setwd("~/UniBern/rnaseq_course/05_Exploratory_data_analysis")



# Reformat the table from FeatureCounts to correspond to the format expected by DESeq2.
# Remove the first line and the columns containing Chr, Start, End, Strand and Length.
# In bash, you could achieve this by combining tail and cut.
# This step was done using the bash script '01_Result_reformatting.sh'.



# Read in the reformatted counts data:
counts_data <- read.delim('./01_Reformat_Results/reformatted_reads_counts.tsv')

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
coldata <- read.delim('./00_Other/Metadata.tsv')
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
dds



# Pre-filtering: Remove the rows with low gene counts across all samples.
# This decreses the size of the DESeq2 object and increases computation speed.
# We set the threshold at 10 counts per row.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds



# Run the DESeq2::DESeq() function and save to variable dds.
# Differential expression analysis.
dds <- DESeq2::DESeq(dds)



# Remove the dependence of the variance on the mean using DESeq2::vst() and set blind=TRUE.
vsd <- DESeq2::vst(dds, blind=T)



# Assess how the samples cluster using DESeq2::plotPCA().
DESeq2::plotPCA(vsd, intgroup='Group')
