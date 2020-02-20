import os

configfile: "/SAN/vyplab/alb_projects/pipelines/rna_seq_snakemake/config/config.yaml"
cluster_config: "/SAN/vyplab/alb_projects/pipelines/rna_seq_snakemake/config/cluster.yaml"
include: "/SAN/vyplab/alb_projects/pipelines/rna_seq_snakemake/rules/helpers.py"

kallisto_output_folder = config["project_top_level"] + "kallisto/"

#make sure the output folder for Kallisto exists before running anything
os.system("mkdir -p {0}".format(kallisto_output_folder))

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
SPECIES_VERSION = get_species_version(config['species'])
GENOME_DIR = os.path.join(config['kallisto_indices'],config['species'],SPECIES_VERSION)
print(GENOME_DIR)


rule all_build:
	input:
		expand(os.path.join(kallisto_output_folder,"{name}", "abundance.h5"), name = SAMPLE_NAMES),
		GENOME_DIR + ".gencode.idx"

rule generate_genome:
	input:
		fasta = get_genome_fasta(config['species']),
		gtf = get_gtf(config['species'])
	output:
		GENOME_DIR + ".gencode.idx"
	threads:
		4
	shell:
		"""
		set +u;
		source activate kallisto
		kallisto index -i {output} {input.fasta}
		"""

if config['end_type'] == "pe":
	rule run_kallisto:
		input:
			generated_index = GENOME_DIR + ".gencode.idx",
			one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz",
			two = config['merged_fastq_folder'] + "{name}_2.merged.fastq.gz"
		output:
			os.path.join(kallisto_output_folder,"{name}", "abundance.h5"),
			os.path.join(kallisto_output_folder,"{name}", "abundance.tsv"),
			os.path.join(kallisto_output_folder,"{name}", "run_info.json")
		params:
			outputPrefix = os.path.join(kallisto_output_folder,"{name}")
		threads:
			4
		shell:
			"""
			set +u;
			source activate kallisto
			mkdir -p {params.outputPrefix}
			kallisto quant -i {input.generated_index} -o {params.outputPrefix} {input.one} {input.two}
			"""
if config['end_type'] == "se":
	rule run_kallisto:
		input:
			generated_index = GENOME_DIR + ".gencode.idx",
			one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz"
		output:
			os.path.join(kallisto_output_folder,"{name}", "abundance.h5"),
			os.path.join(kallisto_output_folder,"{name}", "abundance.tsv"),
			os.path.join(kallisto_output_folder,"{name}", "run_info.json")
		params:
			outputPrefix = os.path.join(kallisto_output_folder,"{name}"),
			strandness = get_kallisto_strand(config['feature_counts_strand_info'])
		threads:
			4
		shell:
			"""
			set +u;
			source activate kallisto
			mkdir -p {params.outputPrefix}
			kallisto quant -i {input.generated_index} -o {params.outputPrefix} --single {input.one}
			"""
