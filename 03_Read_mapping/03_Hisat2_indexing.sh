#!/usr/bin/env bash

#SBATCH --job-name=Hisat2_indexing
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end

#Define the path to the diffrent directories:
INPUTDIR=./01_Reference_genome_sequence_and_annotation
OUTPUTDIR=./03_Hisat2_Index_files

#Produce all required index files with Hisat2:
apptainer exec \
 --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
 hisat2-build $INPUTDIR/Mus_musculus.GRCm39.dna.primary_assembly.fa $OUTPUTDIR/Hisat2_Index