#!/usr/bin/env bash

#SBATCH --job-name=Hisat2_mapping
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=janosch.imhof@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --array=1-16


# Define the path to the diffrent directories:
INDEX_BASENAME_DIR=./03_Hisat2_Index_files/Hisat2_Index
SAMPLELIST=./02_Samplelist/samplelist.tsv
OUTPUTDIR=./04_Hisat2_Mapping

# Get the sample name alongside read1 and read2 file path out of samplelist.tsv
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# Map the reads, must be done separately for each strand. Strandedness setting=RF
apptainer exec \
 --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2 \
 -x $INDEX_BASENAME_DIR -1 $READ1 -2 $READ2 -S $OUTPUTDIR/$SAMPLE.sam -p 8 \
 --rna-strandness RF --summary-file $OUTPUTDIR/${SAMPLE}_summary_file.txt

# Sam files to bam
apptainer exec \
 --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
 samtools view -hbS $OUTPUTDIR/$SAMPLE.sam > $OUTPUTDIR/$SAMPLE.bam

# Delete Sam files
rm $OUTPUTDIR/$SAMPLE.sam