import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"


#include: "rules/fastp.smk"
#RULE ORDER DIRECTIVE
#if paired end, use the paired end rule to run, if single end use the single end rule to run
if config['end_type'] == "pe":
	ruleorder: run_star_pe > run_star_se
else:
	ruleorder: run_star_se > run_star_pe
#make sure the output folder for STAR exists before running anything
os.system("mkdir -p {0}".format(config["star_output_folder"]))

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
SPECIES_VERSION = get_species_version(config['species'])

#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species
GENOME_DIR = os.path.join(config['STAR_indices'],config['species'],SPECIES_VERSION,"star_indices_overhang" + str(config['readLen']))
print(GENOME_DIR)

rule all_samtools:
	input:
		config['star_output_folder'] + "{name}.Aligned.sorted.out.bam",
		config['star_output_folder'] + "{name}.Aligned.sorted.out.bam.bai"

rule run_star_pe:
	input:
		one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz",
		two = config['merged_fastq_folder'] + "{name}_2.merged.fastq.gz"
	output:
		config['star_output_folder'] + "{name}.SJ.out.tab",
		config['star_output_folder'] + "{name}.Log.final.out",
		temp(config['star_output_folder'] + "{name}.Aligned.out.bam")
	params:
		extra_star_parameters = return_parsed_extra_params(config['extra_star_parameters']),
		genomeDir = GENOME_DIR,
		outTmpDir = os.path.join(config['star_output_folder'] + "{name}_tmpdir"),
		outputPrefix = os.path.join(config['star_output_folder'] + "{name}."),
	threads:
		4
	shell:
		"""
		rm -rf {params.outTmpDir}
		{config[star_path]} --genomeDir {params.genomeDir} \
		--readFilesIn {input.one} {input.two} \
		--outFileNamePrefix {params.outputPrefix} \
		--readFilesCommand zcat --runThreadN {threads} \
		{params.extra_star_parameters} \
		--outTmpDir {params.outTmpDir}
		"""

rule run_star_se:
	input:
		one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz"
	output:
		config['star_output_folder'] + "{name}.SJ.out.tab",
		config['star_output_folder'] + "{name}.Log.final.out",
		temp(config['star_output_folder'] + "{name}.Aligned.out.bam")
	params:
		extra_star_parameters = return_parsed_extra_params(config['extra_star_parameters']),
		genomeDir = GENOME_DIR,
		outTmpDir = os.path.join(config['star_output_folder'] + "{name}_tmpdir"),
		outputPrefix = os.path.join(config['star_output_folder'] + "{name}."),
	threads:
		4
	shell:
		"""
		rm -rf {params.outTmpDir}
		{config[star_path]} --genomeDir {params.genomeDir} \
		--readFilesIn {input.one} \
		--outFileNamePrefix {params.outputPrefix} \
		--readFilesCommand zcat --runThreadN {threads} \
		{params.extra_star_parameters} \
		--outTmpDir {params.outTmpDir}
		"""
		
rule sort_bams:
	input:
		config['star_output_folder'] + "{name}.Aligned.out.bam"
	output:
		config['star_output_folder'] + "{name}.Aligned.sorted.out.bam"
	shell:
		"""
		{config[samtools_path]} sort {input} -o {output}
		"""
rule sort_index_bams:
	input:
		config['star_output_folder'] + "{name}.Aligned.sorted.out.bam"
	output:
		config['star_output_folder'] + "{name}.Aligned.sorted.out.bam.bai"
	shell:
		"""
		{config[samtools_path]} index {input}
		"""
