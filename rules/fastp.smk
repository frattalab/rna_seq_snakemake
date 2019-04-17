
configfile: "config/config.yaml"

#fastq name is the sample_unit_fastqfile 
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])
#zip them into a directory to make getting the location easier
SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")

units = samples.unit
fast2 = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in samples.fast2]
ORDER_UNIT_2 = dict(zip(units,fast2))



rule all_trimmed:
	input: 
		expand(config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq",zip, fastq_name=FASTQ_NAME, unit=UNITS)

rule fastp_trimming:
	input:
	#if the input is single end, simply get the value in the fast1 column
		if config["end_type"] == "se":
			fastq_file = lambda wildcards: SAMPLES.loc[SAMPLES['unit'] == wildcards.unit]["fast1"].values[0]
		if config["end_type"] == "pe":
			fastq_file1 = lambda wildcards: SAMPLES.loc[SAMPLES['unit'] == wildcards.unit]["fast1"].values[0],
			fastq_file2 = lambda wildcards: SAMPLES.loc[SAMPLES['unit'] == wildcards.unit]["fast2"].values[0]
	output:
		out_fastqc = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq"
	params:
		fastp_parameters = return_parsed_fastp_params(config['fastp_parameters']),
		if config["end_type"] == "pe":
			 out_fastqc2 = lambda wildcards: config["fastp_trimmed_output_folder"] + "{unit}/" + ORDER_UNIT_2[wildcards.unit], + "_trimmed.fastq"
	run:
		if config["end_type"] == "se":
			shell("{config[fastp_path]} -i {input.fastq_file} -o {output.out_fastqc} {params.fastp_parameters}")
		if config["end_type"] == "pe":
			shell("{config[fastp_path]} -i {input.fastq_file1} -in2 {input.fastq_file1} -o {output.out_fastqc2} -o2 {params.out_fastqc2} {params.fastp_parameters}")
