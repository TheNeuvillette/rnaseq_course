#!/usr/bin/env bash

#SBATCH --job-name=FastQC
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --array=1-16


# define the path to the diffrent directories, including the container directory
WORKDIR=./01_Reads
SAMPLELIST=./02_Samplelist/samplelist.tsv
OUTDIR=./03_FastQC
CONTAINERDIR=/containers/apptainer/fastqc-0.12.1.sif

# get the sample name alongside read1 and read2 file path out of samplelist.tsv
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# define the output file
OUTFILE=$OUTDIR/${SAMPLE}.txt

# run fastqc in the container using appainter
# the bind option gives access to the working directory 
apptainer exec --bind /data/users/jimhof/rnaseq_course/02_Quality_check/01_Reads $CONTAINERDIR \
 fastqc -t 2 -o $OUTDIR $READ1 $READ2 >> $OUTFILE 2>&1