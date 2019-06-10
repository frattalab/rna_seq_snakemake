import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "helpers.py"

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
#make sure the output folder for featureCounts exists before running anything
os.system("mkdir -p {0}".format(config["feature_counts_output_folder"]))

rule all_samtools:

rule feature_counts:
	input:
		config['star_output_folder'] + "{name}.Aligned.sorted.out.bam",
		config['star_output_folder'] + "{name}.Aligned.sorted.out.bam.bai"
	output:
		config['feature_counts_output_folder'] + "{name}_featureCounts_results.txt"