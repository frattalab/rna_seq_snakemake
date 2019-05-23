

import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "helpers.py"

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()


rule all_samtools:
	input:
		expand(config['star_output_folder'] + "{name}/{name}.Aligned.sorted.out.bam",name = SAMPLE_NAMES),
		expand(config['star_output_folder'] + "{name}/{name}.Aligned.sorted.out.bam.bai",name = SAMPLE_NAMES)

rule sort_bams:
	input:
		config['star_output_folder'] + "{name}/{name}_2pass.Aligned.out.bam"
	output:
		config['star_output_folder'] + "{name}/{name}_2pass.Aligned.sorted.out.bam"
	shell:
		"""
		{config[samtools_path]} sort {input} -o {output}
		"""
rule sort_index_bams:
	input:
		config['star_output_folder'] + "{name}/{name}_2pass.Aligned.sorted.out.bam"
	output:
		config['star_output_folder'] + "{name}/{name}_2pass.Aligned.sorted.out.bam.bai"
	shell:
		"""
		{config[samtools_path]} index {input}
		"""