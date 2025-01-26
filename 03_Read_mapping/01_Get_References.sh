#!/usr/bin/env bash

#SBATCH --job-name=Get_References
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end

wget -P ./01_Reference_genome_sequence_and_annotation https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget -P ./01_Reference_genome_sequence_and_annotation https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz

gunzip ./01_Reference_genome_sequence_and_annotation/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip ./01_Reference_genome_sequence_and_annotation/Mus_musculus.GRCm39.113.gtf.gz