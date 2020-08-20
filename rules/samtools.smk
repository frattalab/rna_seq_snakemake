

import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])

rule all_samtools:
	input:
		star_outdir + "{name}.Aligned.sorted.out.bam",
		star_outdir + "{name}.Aligned.sorted.out.bam.bai"

rule sort_bams:
	input:
		star_outdir + "{name}.Aligned.out.bam"
	output:
		star_outdir + "{name}.Aligned.sorted.out.bam"
	shell:
		"""
		{config[samtools_path]} sort {input} -o {output}
		"""
rule sort_index_bams:
	input:
		star_outdir + "{name}.Aligned.sorted.out.bam"
	output:
		star_outdir + "{name}.Aligned.sorted.out.bam.bai"
	shell:
		"""
		{config[samtools_path]} index {input}
		"""
