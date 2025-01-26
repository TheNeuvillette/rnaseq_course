#!/bin/bash

#SBATCH --job-name=Get_Samplelist
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end

FASTQ_FOLDER=./../02_Quality_check/01_Reads

for FILE in $FASTQ_FOLDER/*_*1.fastq.gz
do 
    PREFIX="${FILE%_*.fastq.gz}"
    SAMPLE=`basename $PREFIX`
    echo -e "${SAMPLE}\t$FILE\t${FILE%?.fastq.gz}2.fastq.gz" >> ./02_Samplelist/samplelist.tsv
done