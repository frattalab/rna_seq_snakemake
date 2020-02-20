import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "/SAN/vyplab/alb_projects/pipelines/rna_seq_snakemake/rules/helpers.py"

kallisto_output_folder = config["project_top_level"] + "kallisto"

#make sure the output folder for Kallisto exists before running anything
os.system("mkdir -p {0}".format(kallisto_output_folder))

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
SPECIES_VERSION = get_species_version(config['species'])
GENOME_DIR = os.path.join(config['kallisto_indices'],config['species'],SPECIES_VERSION)
print(GENOME_DIR)
#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species

print(GENOME_DIR)

rule all_build:
	input:
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
