#!/bin/bash
# Job name:
#SBATCH --job-name=amm_rna
#
# Account:
#SBATCH --account=co_rosalind
#
# Partition:
#SBATCH --partition=savio2_htc
#
# Quality of Service:
#SBATCH --qos=rosalind_htc2_normal
#
# Wall clock limit:
#SBATCH --time=2-00:00:00
#
# Nodes:
#SBATCH --nodes=1
# Specify number of tasks for use case (example):
#SBATCH --ntasks-per-node=10
## Command(s) to run:
module load fastqc/0.11.9
module load samtools/1.8
module load multiqc/1.9
snakemake -j --rerun-incomplete