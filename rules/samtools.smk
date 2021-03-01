import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])

flagstat_outdir = get_output_dir(config["project_top_level"], config["samtools_output_folder"])

rule all_samtools:
	input:
		star_outdir + "{name}.Aligned.sorted.out.bam",
		star_outdir + "{name}.Aligned.sorted.out.bam.bai",
		flagstat_outdir + "{name}.flagstat.txt"


rule sort_bams:
	input:
		star_outdir + "{name}.Aligned.out.bam"

	output:
		star_outdir + "{name}.Aligned.sorted.out.bam"

	conda:
		"../env/align.yaml"

	shell:
		"""
		samtools sort {input} -o {output}
		"""

rule sort_index_bams:
	input:
		star_outdir + "{name}.Aligned.sorted.out.bam"

	output:
		star_outdir + "{name}.Aligned.sorted.out.bam.bai"

	conda:
		"../env/align.yaml"

	shell:
		"""
		samtools index {input}
		"""


rule flagstat:
	input:
		star_outdir + "{name}.Aligned.sorted.out.bam",
		star_outdir + "{name}.Aligned.sorted.out.bam.bai"

	output:
		flagstat_outdir + "{name}.flagstat.txt"

	conda:
		"../env/align.yaml"

	shell:
		"""
		samtools flagstat {input[0]} > {output}
		"""
