#!/usr/bin/env bash

#SBATCH --job-name=setup_directories
#SBATCH --time=00:01:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=200M
#SBATCH --partition=pibu_el8

mkdir -p 00_Other
mkdir -p 01_Reference_genome_sequence_and_annotation
mkdir -p 02_Samplelist
mkdir -p 03_Hisat2_Index_files
mkdir -p 04_Hisat2_Mapping
mkdir -p 05_Samtool_Sorted_and_Indexed
