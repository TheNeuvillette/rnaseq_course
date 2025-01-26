#!/usr/bin/env bash

#SBATCH --job-name=Reformatting_results
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end

# Copy the FeatureCount results:
cp ../04_Read_counting/01_Read_Counting/reads_counts.txt ./01_Reformat_Results/reads_counts.txt
cp ../04_Read_counting/01_Read_Counting/reads_counts.txt.summary ./01_Reformat_Results/reads_counts.txt.summary

# Reformat the table from FeatureCounts to correspond to the format expected by DESeq2.
# Specifically, the first line and the columns containing Chr, Start, End, Strand and Length must be removed.
# In bash, you could achieve this by combining tail and cut.

# Input and output file names
input_file=./01_Reformat_Results/reads_counts.txt
output_file=./01_Reformat_Results/reformatted_reads_counts.tsv

# Remove the first line and the specified columns
tail -n +2 "$input_file" | cut --complement -f2-6 > "$output_file"

