import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "helpers.py"
#RULE ORDER DIRECTIVE soo cool. This says that you should try to run the paired end rule beofre the single end rule to produce the output, but if 
#the proper input isn't there for the pairend end e.g. two trimmed fastqs, then run the single end.
ruleorder: run_star_pe > run_star_se
#make sure the output folder for STAR exists before running anything
os.system("mkdir -p {0}".format(config["star_output_folder"]))

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
UNITS = SAMPLES['unit'].tolist()


rule all_star:
	input:
		expand(config['star_output_folder'] + "{name}/{unit}.bam",zip, unit = UNITS,name = SAMPLE_NAMES)

rule run_star_pe:
	input:
		one = lambda wildcards: get_trimmed(wildcards.unit, wildcards.name)[0],
		two = lambda wildcards: get_trimmed(wildcards.unit, wildcards.name)[1]

	output:
		config['star_output_folder'] + "{name}/{unit}.bam" 
	params:
		fastp_parameters = return_parsed_extra_params(config['fastp_parameters'])
	threads: 12

	run:
		cmd = ["{config[star_path]} \
		--genomeDir {params.genomeDir} \
		--readFilesIn {params.trimmed_fastqs} \
		--outFileNamePrefix /path/to/output/dir/prefix \
		--readFilesCommand zcat \
		--runThreadN {threads} \
		{params.extra_star_parameters}"]
		print(cmd)
		shell(cmd)

rule run_star_se:
	input:
		one = lambda wildcards: get_trimmed(wildcards.unit, wildcards.name)
	output:
		config['star_output_folder'] + "{name}/{unit}.bam" 
	params:
		fastp_parameters = return_parsed_extra_params(config['fastp_parameters'])
	threads: 12

	run:
		cmd = ["{config[star_path]} \
		--genomeDir {params.genomeDir} \
		--readFilesIn {params.trimmed_fastqs} \
		--outFileNamePrefix /path/to/output/dir/prefix \
		--readFilesCommand zcat \
		--runThreadN {threads} \
		{params.extra_star_parameters}"]
		print(cmd)
		shell(cmd)