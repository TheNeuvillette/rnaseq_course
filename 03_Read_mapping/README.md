# RNA-sequencing course - Part 3: Read mapping


## Aim of part 3:
The aim of this part is to map the reads (Illumina sequencing data) of the samples to the reference genome to find out which genes are expressed.


## Procedure and organization of this part:
1. The mouse reference genome and the associated annotation are downloaded from the [Ensembl ftp site](https://www.ensembl.org/info/data/ftp/index.html). To double check the file integrity, checksums is computed of both files using "sum 'filename'".
2. A list containing the sample numbers as well as the path to the r1 and r2 file of the samples was created.
3. All required index files were created using Hisat2.
4. All samples were mapped to the reference genome using Hisat2. In the same step, the sam files were converted into bam files and the sam files were deleted. The mapping of the multiple samples was run in parallel using SLURM job arrays.
5. The bam files were sorted and indexed using samtools. Again, the multiple samples were run in parallel using SLURM job arrays.


## RNA-sequencing folder structure:
- Each step of the project has its own direcory.
- Folder structure created with 00_Setup_directories_of_03_Read_Mapping.sh