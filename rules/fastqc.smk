
import os
include: "helpers.py"

configfile: "config/config.yaml"

#get output folder
fastqc_outdir = get_output_dir(config["project_top_level"], config["fastqc_output_folder"])

#make sure the output folder for fastqc exists before running anything
os.system("mkdir -p {0}".format(fastqc_outdir))


#fastq name is the sample_unit_fastqfile
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])

# FASTQ_NAME = sample ID, FILE_LOCATION = path to fastq file
# Fastqc gets it's filename by stripping ".fastq.gz" from the input FASTQ
# FASTQ_PREFIX is a list of <prefix> for each sample (from /path/to/<prefix>.fastq.gz)
FASTQ_PREFIX = [re.sub(".fastq.gz","", location.rpartition('/')[2]) for location in FILE_LOCATION]

# print("FASTQ_NAME values - {0}".format(", ".join(FASTQ_NAME)))
#zip them into a directory to make getting the location easier
ORDER_DICT = dict(zip(FASTQ_NAME, FILE_LOCATION))
# print("\nThis is ORDER_DICT\n")
# print(ORDER_DICT)

#

#first rule is a general rule that specifies the final output of everything, here we have the expected
#output of the individual fastqc's and the multiqc html file. snakemake will check to see if these output files
#exist if they dont it will keep going

# all_files = expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)

rule all_fstq:
	input:
		expand(fastqc_outdir + "{sample}/{unit}/{fastq_prefix}_fastqc.html",zip, sample = FASTQ_NAME, unit=UNITS, fastq_prefix = FASTQ_PREFIX),
		#fastqc_outdir + "fastqc_multiqc_report.html"

rule fastqc:
	input:
		fastq_file = lambda wildcards: return_fastq_location(wildcards.sample),
		#fq_name = lambda wildcards, input: re.sub(".fastq.gz","", ORDER_DICT[input.fastq_file].rpartition('/')[2]),
	output:
		out_fastqc = fastqc_outdir + "{sample}/{unit}/{fastq_prefix}_fastqc.html"
	params:
		#fq_name = lambda wildcards, input: re.sub(".fastq.gz","", ORDER_DICT[input.fastq_file].rpartition('/')[2]),
		outdir = fastqc_outdir + "{sample}/{unit}/"

	conda:
		"../env/align.yaml"

	shell:
		"""
		mkdir -p {params.outdir}
		fastqc {input.fastq_file} -o {params.outdir}
		"""
