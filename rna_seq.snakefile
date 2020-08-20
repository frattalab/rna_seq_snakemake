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

#Construct all output directory strings needed for workflow definition
fastqc_outdir = get_output_dir(config["project_top_level"], config["fastqc_output_folder"])
interleaved_outdir = get_output_dir(config['project_top_level'], config['interleave_master_output_folder'])
star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])
feature_counts_outdir = get_output_dir(config["project_top_level"], config["feature_counts_output_folder"])

#save typing
workflow = config['workflow']

if workflow not in ['fastq_qc','align','interleave_fastq_qc']:
	raise ValueError("{0} is an invalid value for 'workflow' - must be one of {1}".format(workflow, ', '.join(['fastq_qc','align','interleave_fastq_qc'])))

elif workflow == "fastq_qc":
		rule all:
			input:
				expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)
			shadow: "minimal"

		include: "rules/fastp.smk"
		include: "rules/fastqc.smk"

elif workflow == "interleave_fastq_qc":
		rule all:
		    input:
		        expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS),
				expand(interleaved_outdir + "{name}_interleaved.fastq.gz", name = SAMPLE_NAMES)
			shadow: "minimal"

		include: "rules/fastp.smk"
		include: "rules/fastqc.smk"
		include: "rules/interleave_fastqs.smk"


elif workflow == "align":
		rule all:
			input:
				GENOME_DIR + "/SA",
				expand(feature_counts_outdir + "{name}_featureCounts_results.txt", name = SAMPLE_NAMES),
				expand(star_outdir + "{name}.Aligned.sorted.out.bam", name = SAMPLE_NAMES),
				expand(star_outdir + "{name}.Aligned.sorted.out.bam.bai", name = SAMPLE_NAMES),
				# expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)
				#expand(config['project_top_level'] + "multiqc_report.html")
			shadow: "minimal"

		#Not necessary to specify include rules each time but doing to make wf easier to follow
		include: "rules/fastp.smk"
		include: "rules/fastqc.smk"
		include: "rules/generate_star_index.smk"
		include: "rules/star.smk"
		include: "rules/feature_counts.smk"



#rule all:
#	input:
#		GENOME_DIR + "/SA",
#		expand(config['feature_counts_output_folder'] + "{name}_featureCounts_results.txt", name = SAMPLE_NAMES),
#		expand(config['star_output_folder'] + "{name}.Aligned.sorted.out.bam",name = SAMPLE_NAMES),
#		expand(config['star_output_folder'] + "{name}.Aligned.sorted.out.bam.bai", name = SAMPLE_NAMES),
#		# expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)
#		#expand(config['project_top_level'] + "multiqc_report.html")
#	shadow: "minimal"




#include: "rules/fastqc.smk"
#include: "rules/fastp.smk"
#include: "rules/generate_star_index.smk"
#include: "rules/star.smk"
#include: "rules/feature_counts.smk"
