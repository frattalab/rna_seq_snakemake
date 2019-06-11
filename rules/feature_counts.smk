import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "helpers.py"

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
#make sure the output folder for featureCounts exists before running anything
os.system("mkdir -p {0}".format(config["feature_counts_output_folder"]))
#this function uses the text file located in the config folder "star_genomes_species.csv" and 
#the config file species parameter to 
#give the correct genome for the species
REFERENCE_ANNOTATION = get_gtf(config['species'])
		
rule all_featurecounts:
	input:
		expand(config['feature_counts_output_folder'] + "{name}_featureCounts_results.txt", name = SAMPLE_NAMES)

rule feature_counts:
	input:
		aligned_bam = config['star_output_folder'] + "{name}.Aligned.sorted.out.bam",
		aligned_bai = config['star_output_folder'] + "{name}.Aligned.sorted.out.bam.bai"
	output:
		out_name = config['feature_counts_output_folder'] + "{name}_featureCounts_results.txt"
	params:
		ref_anno = REFERENCE_ANNOTATION
	shell:
		"""
		{config[feature_counts_path]} -a {params.ref_anno} -o {output.out_name} {input.aligned_bam}
		"""