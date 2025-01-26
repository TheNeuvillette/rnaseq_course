#!/usr/bin/env bash

#SBATCH --job-name=Read_counting
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end

# Define the path to the diffrent directories:
ANNOTATION_FILE=./../03_Read_mapping/01_Reference_genome_sequence_and_annotation/Mus_musculus.GRCm39.113.gtf
INPUTDIR=./../03_Read_mapping/05_Samtool_Sorted_and_Indexed
OUTPUTDIR=./01_Read_Counting

#Index the coordinate sorted bam files
apptainer exec \
 --bind /data/ /containers/apptainer/subread_2.0.6.sif \
 featureCounts -p -s 2 -T 8 -t 'exon' -g 'gene_id' -a $ANNOTATION_FILE -o $OUTPUTDIR/reads_counts.txt $INPUTDIR/*.bam