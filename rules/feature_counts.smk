import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

feature_counts_outdir = get_output_dir(config["project_top_level"], config["feature_counts_output_folder"])
star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])

#make sure the output folder for featureCounts exists before running anything
os.system("mkdir -p {0}".format(feature_counts_outdir))
#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species
REFERENCE_ANNOTATION = get_gtf(config['species'])

rule all_featurecounts:
	input:
		expand(feature_counts_outdir + "{name}_featureCounts_results.txt", name = SAMPLE_NAMES)

rule feature_counts:
	input:
		aligned_bam = star_outdir + "{name}.Aligned.sorted.out.bam",
		aligned_bai = star_outdir + "{name}.Aligned.sorted.out.bam.bai"
	output:
		out_name = feature_counts_outdir + "{name}_featureCounts_results.txt"
	params:
		ref_anno = REFERENCE_ANNOTATION,
		stranded = config['feature_counts_strand_info']
	run:
		if config["end_type"] == "pe":
			shell("{config[feature_counts_path]} -p -t exon -g gene_id -a {params.ref_anno}  --extraAttributes gene_name -o {output.out_name} {params.stranded} {input.aligned_bam}")
		if config["end_type"] == "se":
			shell("{config[feature_counts_path]} -a {params.ref_anno} --extraAttributes gene_name -o {output.out_name} {params.stranded} {input.aligned_bam}")
