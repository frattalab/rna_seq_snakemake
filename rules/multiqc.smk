

FASTQ_NAME = get_fastq_names(config["sampleCSVpath"])
def check_all_done(FASTQ_NAME):

	
rule all:
	input:
		config["fastqc_output_folder"] + "multiqc.html"

rule multiqc:
    input:
        check_all_done
    output:
        config["fastqc_output_folder"] + "multiqc.html"
    log:
        "logs/multiqc.log"
    shell:
        "{config[multiqc_path]} {config[fastqc_output_folder]}"