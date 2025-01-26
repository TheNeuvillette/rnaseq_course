# RNA-sequencing course - Part 5: Exploratory data analysis


## Aim of part 5:
The aim of this part is to get a first impression and do quality control of the FeatureCount results. Specifically, we must check if the samples from the same experimental group show similar gene expression patterns.


## Procedure and organization of this part:
1. Reformat the FeatureCount results to make it compatible with DESeq2. For that, the header (1st line) and several columns are removed using a bash script.
2. In an R script, the reformatted FeatureCount results are used to do a PCA Plot using DESeq2. To do so, a DESeqDataSet object is first created from the reformatted FeatureCount results. Then, a differential expression analysis is done followed by removing the dependence of the variance. This then allows to plot the PCA.


## RNA-sequencing folder structure:
- Each step of the project has its own direcory.
- Folder structure created with 00_Setup_directories_of_03_Read_Mapping.sh
