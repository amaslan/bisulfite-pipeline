#!/bin/bash
# Job name:
#SBATCH --job-name=amm_rna
#
# Account:
#SBATCH --account=co_rosalind
#
# Partition:
#SBATCH --partition=savio
#
# Quality of Service:
#SBATCH --qos=rosalind_savio_normal
#
# Wall clock limit:
#SBATCH --time=2-00:00:00
#
# Nodes:
#SBATCH --nodes=1
## Command(s) to run:
module load fastqc/0.11.9
module load samtools/1.8
module load multiqc/1.9
snakemake -j --rerun-incomplete