'''
Perform QC and bisulfite analysis with bismark

Requirements:
- update config.yaml
- fastqc, bismark, bowtie
- fasta of reference genome

Usage: snakemake -j
Include -j to run with number of available CPU cores in the machine

dry run:
snakemake -np

work flow diagram:
snakemake --forceall --dag | dot -Tpng > dag.png

'''
# refer to: https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
# module load bowtie2/2.3.4.1
# module load bismark

# Bismark requires: Perl, Bowtie/Bowtie2
# use --multicore for parallelization
# bowtie if <50, bowtie2 if >50 read length



from os.path import join, basename, dirname
import re

configfile: 'config.yaml'

GENOME_DIR = config['GENOME_DIR']
DATA_DIR = config['DATA_DIR']
OUT_DIR = config['OUT_DIR']



# FUNCTION: sample_dictionary
# PARAMETER: fastq -> path to a folder with fastq files (FASTQ_DIR)
# RETURN: 
    # dictionary -> key: sample name ; value: dictionary where keys are 'R1' and/or 'R2' and values are the sample's full path
    # list_of_R1_and_R2 --> list of all samples without '_001.fastq.gz' ending

dictionary = {}
list_of_R1_and_R2 = []

def sample_dictionary(fastq):
	for root, dirs, files in os.walk(fastq):
		global dictionary
		global list_of_R1_and_R2
		for sample_files in files:
			# change based on whether reference or other data depending on naming format
			if bool(re.search('.fq.gz', sample_files)):
				sample_complete = sample_files.split('.fq.gz')
				sample_info = sample_files.split('.')
				sample_name = sample_info[0]
				if sample_name in dictionary: 
					full_path = root + '/' + sample_files
					dictionary[sample_name][sample_info[1]] = [full_path]
					list_of_R1_and_R2.append(full_path)
				else: 
					dictionary[sample_name] = {}
					full_path = root + '/' + sample_files
					dictionary[sample_name][sample_info[1]] = [full_path]
					list_of_R1_and_R2.append(full_path) 
		return dictionary, list_of_R1_and_R2

FILES, SAMPLES_FULL_PATH = sample_dictionary(DATA_DIR)
SAMPLES = sorted(FILES.keys())

sample_names = []

def R1_R2_sample_names(list_r1_r2):
    global sample_names
    for i in range(0, len(list_r1_r2)):
        name_full = list_r1_r2[i].split('/')[-1]
        name = name_full.rstrip('.fq.gz')
        sample_names.append(name)
    return sorted(sample_names)

SAMPLES_FULL = R1_R2_sample_names(SAMPLES_FULL_PATH)

rule all:
	input:
		# fastqc
		[OUT_DIR + "/" + x for x in expand('fastQC_output/{sample_full}_fastqc.html', sample_full = SAMPLES_FULL)],
		# multiqc
		OUT_DIR + '/quality_control_metrics/multiqc/multiqc_report.html'

# 1. fastqc & multiqc
rule fast_qc:
	input:
		lambda wildcards: [s for s in SAMPLES_FULL_PATH if wildcards.sample_full in s]
	output:
		'{OUT_DIR}/fastQC_output/{sample_full}_fastqc.html'
	log:
		'{OUT_DIR}/logs/{sample_full}_fastqc.log',
	shell:
		"""
		fastqc -o {OUT_DIR}/fastQC_output {input} &> {log}
		"""

rule multi_qc_metrics:
    input:
        fastqc_files = expand(join(OUT_DIR, 'fastQC_output', '{sample_full}_fastqc.html'), sample_full = SAMPLES_FULL), 
    output:
        '{OUT_DIR}/quality_control_metrics/multiqc/multiqc_report.html'
    shell:
    	"""
		multiqc {OUT_DIR} -o '{OUT_DIR}/quality_control_metrics/multiqc'
		"""

# 2. bismark genome preparation --> genome is bisulfite converted and indexed
# bismark_genome_preparation --bowtie2 --verbose genome_folder 

# 3. Bismark mapping in paired end mode; reference GRCh37
# tolerating n non-bisulfite matches per read
# bismark --bowtie2 -n 1 -l read_length genome_folder -1 .fastq -2 .fastq

# 4. deduplicate reads

# 5. CpG methylation calls extracted from deduplicated mapping output using Bismark methylation extractor (v 0.14.4) in paired-end mode
# bismark_methylation_extractor -p sample.bismark.sam