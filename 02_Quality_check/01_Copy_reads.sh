#!/usr/bin/env bash

#SBATCH --job-name=Copy_reads
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=end

cp /data/courses/rnaseq_course/toxoplasma_de/README /data/users/jimhof/rnaseq_course/02_Quality_check/00_Other
cp /data/courses/rnaseq_course/toxoplasma_de/reads/*.fastq.gz /data/users/jimhof/rnaseq_course/02_Quality_check/01_Reads