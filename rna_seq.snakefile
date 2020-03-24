import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "rules/helpers.py"
SPECIES_VERSION = get_species_version(config['species'])
GENOME_DIR = os.path.join(config['STAR_indices'],config['species'],SPECIES_VERSION,"star_indices_overhang" + str(config['readLen']))
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])
#zip them into a directory to make getting the location easier
SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

rule all:
	input:
		GENOME_DIR + "/SA",
		expand(config['feature_counts_output_folder'] + "{name}_featureCounts_results.txt", name = SAMPLE_NAMES),
		expand(config['star_output_folder'] + "{name}.Aligned.sorted.out.bam",name = SAMPLE_NAMES),
		expand(config['star_output_folder'] + "{name}.Aligned.sorted.out.bam.bai", name = SAMPLE_NAMES),
		# expand(config["fastqc_output_folder"] + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)
		#expand(config['project_top_level'] + "multiqc_report.html")
	shadow: "minimal"




# include: "rules/fastqc.smk"
include: "rules/fastp.smk"
include: "rules/generate_star_index.smk"
include: "rules/star.smk"
include: "rules/feature_counts.smk"
