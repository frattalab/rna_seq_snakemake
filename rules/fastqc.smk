
include: "helpers.py"

configfile: "config/config.yaml"

#get output folder
fastqc_outdir = get_output_dir(config["project_top_level"], config["fastqc_output_folder"])

#make sure the output folder for fastqc exists before running anything
os.system("mkdir -p {0}".format(fastqc_outdir))


#fastq name is the sample_unit_fastqfile
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])
#zip them into a directory to make getting the location easier
ORDER_DICT = dict(zip(FASTQ_NAME, FILE_LOCATION))
#first rule is a general rule that specifies the final output of everything, here we have the expected
#output of the individual fastqc's and the multiqc html file. snakemake will check to see if these output files
#exist if they dont it will keep going

all_files = expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)

rule all_fstq:
	input:
		expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS),
		#fastqc_outdir + "fastqc_multiqc_report.html"

rule fastqc:
	input:
		fastq_file = lambda wildcards: return_fastq_location(wildcards.fastq_name)
	output:
		out_fastqc = fastqc_outdir + "{unit}/{fastq_name}_fastqc.html"
	params:
		fq_name = lambda wildcards, input: re.sub(".fastq.gz","_fastqc", input.fastq_file.rpartition('/')[2]),
		outdir = fastqc_outdir
	log:
        "logs/{fastq_name}_fastqc.log"
	shell:
		"""
		mkdir -p {params.outdir}{wildcards.unit}
		{config[fastqc_path]} {input.fastq_file} -o {params.outdir}{wildcards.unit}
		"""
