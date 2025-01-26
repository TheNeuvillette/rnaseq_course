# RNA-Sequencing Course: Comparative Analysis of *Toxoplasma gondii* Infection on Gene Expression in Mouse Blood and Lung Tissues

Understanding how the immune response is reflected across various organs and the interplay of the organs upon disease is key for effective diagnoses and valuable insight into immune regulation as a system. In this project, the aim is to identify the effect of *Toxoplasma gondii* infection on the gene expression of mouse blood and lung tissues, and to identify gene ontology (GO) terms associated with differentially expressed genes (DE). Furthermore, the project aims to compare gene expression and GO terms of DE genes between blood and lung tissue. The differential expression analysis workflow is based on the sequencing data of the [Singhania, et al.](https://www.nature.com/articles/s41467-019-10601-6) study. The workflow consists of quality checking the sequencing data using FastQC, mapping the sequenced reads to the reference genome via hisat2 and samtools and then counting the number of reads per gene using featureCounts. The produced table of counts per gene is further analyzed and visualized by exploratory data analysis, DE analysis and overrepresentation analysis using DESeq2 as well as clusterProfiler.

## 1. Getting started

The DE analysis workflow in this project uses a subset of the RNA-seq data from the [Singhania, et al.](https://www.nature.com/articles/s41467-019-10601-6) study, with the fastq files available on Gene Expression Omnibus (GEO), [GEO accession GSE119855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119855). This project’s workflow uses the RNA-seq data of lung and blood samples from wild type mice (*Mus Musculus*) and mice infected with the alveolate *Toxoplasma gondii*. The dataset includes 3 mice control samples and 5 infected mice samples for both tissues (details [here](https://github.com/TheNeuvillette/rnaseq_course/blob/main/01_Getting_started/Sample_Information)).

## 2. Quality check
Quality assessment of the fastq files was performed using FastQC [Version 0.12.1] and summarized using MultiQC [Version 1.19], with the scripts being available under [02_Quality_check](https://github.com/TheNeuvillette/rnaseq_course/tree/main/02_Quality_check).

The quality check workflow consisted of several steps, all in their separate script. Firstly, the Illumina sequencing data (fastq files) was copied into the private folder. Secondly, a list containing the sample numbers as well as the path to the r1 and r2 file of the samples was created. Thirdly, the quality of all samples was controlled using FastQC. All files were run in parallel using SLURM job arrays. Finally, a MultiQC was done to summarize the FastQC of all samples.

## 3. Read mapping
To match the raw RNA-sequencing reads to their respective genes, the reads were mapped onto the [Ensembl reference genome of *Mus Musculus*](https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/) [Release 113]. Read mapping was done using HISAT2 (Hierarchical Indexing for Spliced Alignment of Transcripts) [Version 2.2.1] and post-processed using SAMtools [Version 1.20]. The scripts of the read mapping part are available under [03_Read_mapping](https://github.com/TheNeuvillette/rnaseq_course/tree/main/03_Read_mapping).

The read mapping workflow consisted of several steps combining HISAT2 and SAMtools. HISAT2 index files were first generated to perform alignment efficiently, followed by the alignment to the reference genome using the strandedness stetting RF (paired-end reads). Using SAMtools, the SAM files were converted to BAM files, facilitating work with the read alignment files. The SAM files were further sorted and indexed by genomic coordinates, enabling efficient analysis and access to the genome in the downstream analysis.

## 4. Read counting
Using the mapped reads alongside the [*Mus Musculus* genome annotation file](https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/) [Release 113] allowed the quantification of the gene expression levels of all genes. Read counting was done using featureCounts [Version 2.0.6] with the strandedness setting 2. The scripts being available under [04_Read_counting](https://github.com/TheNeuvillette/rnaseq_course/tree/main/04_Read_counting).

## 5. Exploratory data analysis
Exploratory data analysis (EDA) encompassed quality checking and summarization of the characteristics of the gene count data prior to more complex differential gene expression analyses. In particular, it was checked if samples of the same condition showed similar gene expression patterns using principle component analysis (PCA). EDA of the count data from RNA-seq was done in R [Version 4.4.1] using the package DESeq2 [Version 1.46.0]. The scripts of the exploratory data analysis part are available under [05_Exploratory_data_analysis](https://github.com/TheNeuvillette/rnaseq_course/tree/main/05_Exploratory_data_analysis).

The EDA workflow consisted of several steps beginning with the featureCounts output and resulting in a PCA visualization of the gene expression profile of all samples. The featureCounts output was first reformatted, making it compatible with DESeq2, followed by the creation of a DESeqDataSet object. This allowed for the DE analysis and following removal of the variance dependance on the mean, enabling the construction of the PCA plot. 

## 6. Differential expression analysis

Differential expression analysis includes the identification, investigation and visualization of DE genes between different experimental conditions, specifically between infected and control samples. The results of the previously run DE analysis using DESeq2 were extracted for both lung and blood pairwise contrasts (both times infected vs. control) and visualized using a volcano plot using the R package EnhancedVolcano [Version 1.24.0]. The R script used for the differential expression analysis is available under [06_Differential_expression_analysis](https://github.com/TheNeuvillette/rnaseq_course/tree/main/06_Differential_expression_analysis).

## 7. Overrepresentation analysis

To summarize the DE analysis results into gene ontology (GO) terms using the biological process (BP) subontology, the R package clusterProfiler [Version 4.14.4] was used. The R script used for this step is available under [07_Overrepresentation_analysis](https://github.com/TheNeuvillette/rnaseq_course/tree/main/07_Overrepresentation_analysis).

