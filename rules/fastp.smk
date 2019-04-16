#performs trimming used fastp
configfile: "config/config.yaml"


#fastq name is the sample_unit_fastqfile 
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])
#zip them into a directory to make getting the location easier
ORDER_DICT = dict(zip(FASTQ_NAME, FILE_LOCATION))

rule all_trimmed:
	input: 
		expand(config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq",zip, fastq_name=FASTQ_NAME, unit=UNITS)

rule fastp_trimming:
	input:
		fastp_input_parameters = return_fastp_inputs(config["sampleCSVpath"], config["end_type"])
	output:
		out_fastqc = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq"
	params:
		fastp_parameters: return_parsed_fastp_params(config['fastp_parameters'])
	run:
		shell("{config[fastp_path]} {input.fastp_input_parameters} {params.fastp_parameters}")

