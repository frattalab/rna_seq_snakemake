import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "helpers.py"

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()


rule merge_all_trimmed:
	input:
		expand(config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz", name = SAMPLE_NAMES)
		
if config['end_type'] == "pe":
	rule merge_trimmed:
		input:
			one = lambda wildcards: get_trimmed(wildcards.name)[0],
			two = lambda wildcards: get_trimmed(wildcards.name)[1]
		output:
			out_one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz", 
			out_two = config['merged_fastq_folder'] + "{name}_2.merged.fastq.gz"
		params:
			#taking the input files and putting them into a comma separated list
			one = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[0]),
			two = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[1])
		shell:
			"""
			cat {params.two} > {output.out_two}
			cat {params.one} > {output.out_one}
			"""
else:
	rule merge_trimmed:
		input:
			one = lambda wildcards: get_trimmed(wildcards.name)[0]
		output:
			out_one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz"
		params:
			#taking the input files and putting them into a comma separated list
			one = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[0]),
		shell:
			"""
			cat {params.one} > {output.out_one}
			"""

