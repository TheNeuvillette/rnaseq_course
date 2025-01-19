#!/usr/bin/env bash

#SBATCH --job-name=FastQC_Test
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --array=0-15

# define the path to the diffrent directories, including the container directory
WORKDIR=./01_Reads
OUTDIR=./99_FastQC
CONTAINER_PATH=/containers/apptainer/fastqc-0.12.1.sif

# defining read1 and read2
FILES=(./01_Reads)
READ1=${FILES[$SLURM_ARRAY_TASK_ID*2]}
READ2=${FILES[$SLURM_ARRAY_TASK_ID*2+1]}


# fastqc analysis from the container using the --bind flag to give the container the access to the sample files
apptainer exec --bind /data/users/jimhof/rnaseq_course/02_Quality_check/01_Reads $CONTAINER_PATH \
 fastqc -t 2 -o $OUTDIR $READ1 $READ2
