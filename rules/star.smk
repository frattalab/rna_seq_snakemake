import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

# RULE ORDER DIRECTIVE
# if paired end, use the paired end rule to run, if single end use the single end rule to run
if config['end_type'] == "pe":
	ruleorder: run_star_pe > run_star_se
else:
	ruleorder: run_star_se > run_star_pe


#make sure the output folder for STAR exists before running anything
star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])
os.system("mkdir -p {0}".format(star_outdir))

merged_outdir = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])


SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
SPECIES_VERSION = get_species_version(config['species'])

#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species
GENOME_DIR = os.path.join(config['STAR_indices'],config['species'],SPECIES_VERSION,"star_indices_overhang" + str(config['readLen']))

rule all_samtools:
	input:
		expand(star_outdir + "{name}.Aligned.sorted.out.bam",name = SAMPLE_NAMES),
		expand(star_outdir + "{name}.Aligned.sorted.out.bam.bai", name = SAMPLE_NAMES),
		expand(star_outdir + "{name}.flagstat.txt", name = SAMPLE_NAMES)

rule run_star_pe:
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	input:
		generated_index = GENOME_DIR + "/SA",
		one = lambda wildcards: get_processed_fastq(wildcards.name, pair=1),
		two = lambda wildcards: get_processed_fastq(wildcards.name, pair=2)
	output:
		star_outdir + "{name}.SJ.out.tab",
		star_outdir + "{name}.Log.final.out",
		temp(star_outdir + "{name}.Aligned.out.bam")
	params:
		extra_star_parameters = return_parsed_extra_params(config['extra_star_parameters']),
		genomeDir = GENOME_DIR,
		outTmpDir = os.path.join(star_outdir + "{name}_tmpdir"),
		outputPrefix = os.path.join(star_outdir + "{name}.")
	threads:
		4
	conda:
		"../env/align.yaml"

	shell:
		"""
		rm -rf {params.outTmpDir}
		STAR --genomeDir {params.genomeDir} \
		--readFilesIn {input.one} {input.two} \
		--outFileNamePrefix {params.outputPrefix} \
		--readFilesCommand zcat --runThreadN {threads} \
		{params.extra_star_parameters} \
		--outTmpDir {params.outTmpDir}
		"""

rule run_star_se:
	input:
		generated_index = GENOME_DIR + "/SA",
		one = lambda wildcards: get_processed_fastq(wildcards.name, pair=1)
	output:
		star_outdir + "{name}.SJ.out.tab",
		star_outdir + "{name}.Log.final.out",
		temp(star_outdir + "{name}.Aligned.out.bam")
	params:
		extra_star_parameters = return_parsed_extra_params(config['extra_star_parameters']),
		genomeDir = GENOME_DIR,
		outTmpDir = os.path.join(star_outdir + "{name}_tmpdir"),
		outputPrefix = os.path.join(star_outdir + "{name}.")
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	conda:
		"../env/align.yaml"
	threads:
		4
	shell:
		"""
		rm -rf {params.outTmpDir}
		STAR --genomeDir {params.genomeDir} \
		--readFilesIn {input.one} \
		--outFileNamePrefix {params.outputPrefix} \
		--readFilesCommand zcat --runThreadN {threads} \
		{params.extra_star_parameters} \
		--outTmpDir {params.outTmpDir}
		"""

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
		star_outdir + "{name}.flagstat.txt"

	conda:
		"../env/align.yaml"

	shell:
		"""
		samtools flagstat {input[0]} > {output}
		"""
