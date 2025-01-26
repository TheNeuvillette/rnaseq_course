#!/usr/bin/env bash

#SBATCH --job-name=Samtools_sorting
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=30G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --array=1-16

# Define the path to the diffrent directories:
INDEX_BASENAME_DIR=./03_Hisat2_Index_files/Hisat2_Index
SAMPLELIST=./02_Samplelist/samplelist.tsv
INPUTDIR=./04_Hisat2_Mapping
OUTPUTDIR=./05_Samtool_Sorted_and_Indexed

# Get the sample name alongside read1 and read2 file path out of samplelist.tsv
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

#Sort the bam files by genomic coordinates
apptainer exec \
 --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
 samtools sort -m 30G -@ 8 -o $OUTPUTDIR/${SAMPLE}_sorted.bam -T temp $INPUTDIR/${SAMPLE}.bam