
configfile: "config/config.yaml"
include: "helpers.py"

FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])

all_files = expand(config["fastqc_output_folder"] + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)
rule all_output:
	input:
		config["fastqc_output_folder"] + "multiqc_report.html"


rule multiqc:
    input:
        all_files
    output:
        config["fastqc_output_folder"] + "multiqc_report.html"
    log:
        "logs/multiqc.log"
    shell:
        "{config[multiqc_path]} -d {config[fastqc_output_folder]} -o {config[fastqc_output_folder]} {config[multiqc_configparams]}"