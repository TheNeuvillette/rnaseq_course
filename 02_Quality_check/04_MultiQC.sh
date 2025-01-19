#!/usr/bin/env bash

#SBATCH --job-name=FastQC
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end


# define the path to the diffrent directories, including the container directory
WORKDIR=./03_FastQC
OUTDIR=./04_MultiQC
CONTAINERDIR=/containers/apptainer/multiqc-1.19.sif

# run multiqc in the container using appainter
# the bind option gives access to the working directory 
apptainer exec --bind /data/users/jimhof/rnaseq_course/02_Quality_check/ $CONTAINERDIR \
 multiqc -d $WORKDIR -o $OUTDIR