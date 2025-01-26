#!/usr/bin/env bash

#SBATCH --job-name=setup_directories
#SBATCH --time=00:01:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=200M
#SBATCH --partition=pibu_el8

mkdir -p 00_Other
mkdir -p 01_Read_Counting
