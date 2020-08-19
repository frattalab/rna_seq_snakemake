import os
singlestep = TRUE

if singlestep == TRUE:
    project_folder =  "/SAN/vyplab/alb_projects/data/buratti_new_shsy5y/"
    end_type = "pe"
    suffix = ".Aligned.sorted.out"
    star_output_folder = project_folder + "STAR_aligned"
    SAMPLE_NAMES, = glob_wildcards(star_output_folder + "{sample}" + suffix + ".bam")
    print(SAMPLE_NAMES)
else:
    configfile: "config/config.yaml"
    cluster_config: "config/cluster.yaml"
    include: "helpers.py"
    suffix = ".Aligned.sorted.out"
    SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
    SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

    SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
    #make sure the output folder for featureCounts exists before running anything
    os.system("mkdir -p {0}".format(config["feature_counts_output_folder"]))
    star_output_folder = config['star_output_folder']
    end_type = config["end_type"]
#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species
REFERENCE_ANNOTATION = get_gtf(config['species'])

rule tpmcounts:
	input:
		expand(star_output_folder + "{name}" + suffix + "_genes.out", name = SAMPLE_NAMES)

rule tpmcalculator_path:
	input:
		aligned_bam = star_output_folder + "{name}.Aligned.sorted.out.bam",
		aligned_bai = star_output_folder + "{name}.Aligned.sorted.out.bam.bai"
	output:
		star_output_folder + "{name}" + suffix + "_genes.out"
	params:
		ref_anno = REFERENCE_ANNOTATION,
		stranded = config['feature_counts_strand_info']
	run:
		if config["end_type"] == "pe":
			shell("{config[tpmcalculator_path]} -g {params.ref_anno} -b {input.aligned_bam} -p -e -a")
		if config["end_type"] == "se":
			shell("{config[tpmcalculator_path]} -g {params.ref_anno} -b {input.aligned_bam} -e -a")
# rule move_tpm_output:
# 	input:
# 		aligned_bam = star_output_folder + "{name}.Aligned.sorted.out.bam",
# 		aligned_bai = star_output_folder + "{name}.Aligned.sorted.out.bam.bai"
# 	output:
# 		out_name = config['feature_counts_output_folder'] + "{name}_featureCounts_results.txt"
# 	params:
# 		ref_anno = REFERENCE_ANNOTATION,
# 		stranded = config['feature_counts_strand_info']
# 	run:
# 		if config["end_type"] == "pe":
# 			shell("{config[tpmcalculator_path]} -g {params.ref_anno} -b {input.aligned_bam} -p -e -a")
# 		if config["end_type"] == "se":
# 			shell("{config[tpmcalculator_path]} -g {params.ref_anno} -b {input.aligned_bam} -e -a")
