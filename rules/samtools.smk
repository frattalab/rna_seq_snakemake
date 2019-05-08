






rule all_samtools:
	input:
		expand(config['star_output_folder'] + "{name}/{name}.Aligned.sorted.out.bam",name = SAMPLE_NAMES),
		expand(config['star_output_folder'] + "{name}/{name}.Aligned.sorted.out.bam.bai",name = SAMPLE_NAMES)

rule sort_bams:
	input:
	output:
	shell:
		"""
		{}
		"""