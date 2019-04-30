
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "rules/helpers.py"

#make sure the output folder for STAR exists before running anything
os.system("mkdir -p {0}".format(config["star_output_folder"]))

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
UNITS = SAMPLES['unit'].tolist()

rule all_star:
	input:
		expand(config['star_output'] + "{sample_name}/{unit}}_Aligned.out.bam", zip, sample_name=SAMPLE_NAMES, unit=UNITS)

rule run_star:
	input:
		config["fastp_trimmed_output_folder"] + "{unit}/{{fastq_name}}_trimmed.fastq.gz"
	output:
		config['star_output'] + "{sample_name}/{unit}_Aligned.out.bam",
        config['star_output'] + "{sample_name}/{unit}_ReadsPerGene.out.tab"
	params:
		extra_star_parameters = return_parsed_extra_params(config['fastp_parameters']),
		trimmed_fastqs = return_trimmed_fastqs(wilcards.sample_name, wildcards.unit, config['end_type']),
		genomeDir = return_genome_dir_by_species(config['species']) 

	threads: 12

	run:
		cmd = ["{config[star_path]} \
		--genomeDir {params.genomeDir} \
		--readFilesIn {params.trimmed_fastqs} \
		--outFileNamePrefix /path/to/output/dir/prefix \
		--readFilesCommand zcat \
		--runThreadN {threads} \
		{params.extra_star_parameters}"]
		print(cmd)
		shell(cmd)