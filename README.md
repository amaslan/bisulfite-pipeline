# bisulfite-pipeline
bismak bisulfite analysis

### Requirements
- update config.yaml with directory paths
- fastqc, bismark, bowtie2, samtools
- fasta of reference genome

### Overview
1. FastQC
2. MultiQC
3. Bismark genome prep --> genome is bisulfite converted and indexed
4. Bismark mapping in paired end mode
5. Reads are deduplicated with deduplicate_bismark
6. Methylation calls are extracted with bismark_methylation_extractor
