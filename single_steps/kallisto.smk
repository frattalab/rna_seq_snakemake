import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"



#make sure the output folder for Kallisto exists before running anything
os.system("mkdir -p {0}".format(config["kallisto_output_folder"]))

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
SPECIES_VERSION = get_species_version(config['species'])

#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species
GENOME_DIR = config['kallisto_idx']
print(GENOME_DIR)

rule all_kallisto:
	input:
		expand(os.path.join(config['kallisto_output_folder'],"{name}", "abundance.h5"), name = SAMPLE_NAMES)
if config['end_type'] == "pe":
	rule run_kallisto:
		input:
			generated_index = GENOME_DIR,
			one = config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz",
			two = config['merged_fastq_folder'] + "{name}_2.merged.fastq.gz"
		output:
			os.path.join(config['kallisto_output_folder'],"{name}", "abundance.h5"),
			os.path.join(config['kallisto_output_folder'],"{name}", "abundance.tsv"),
			os.path.join(config['kallisto_output_folder'],"{name}", "run_info.json")
		params:
			outputPrefix = os.path.join(config['kallisto_output_folder'],"{name}")
		threads:
			4
		shell:
			"""
			mkdir -p {params.outputPrefix}
			{config[kallisto_path]} quant -i {input.generated_index} -o {params.outputPrefix} {input.one} {input.two}
			"""
