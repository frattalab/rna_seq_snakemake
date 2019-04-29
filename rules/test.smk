
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "helpers.py"

#make sure the output folder for STAR exists before running anything
os.system("mkdir -p {0}".format(config["star_output_folder"]))

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
UNITS = SAMPLES['unit'].tolist()

rule all_star:
	input:
		expand(config['star_output_folder'] + "{sample_name}/{unit}}_temp.bam", zip, sample_name=SAMPLE_NAMES, unit=UNITS)
rule run_star:
	input:
		config["fastp_trimmed_output_folder"] + "{unit}/{{fastq_name}}_trimmed.fastq.gz"
	output:
		config['star_output'] + "{sample_name}/{unit}_temp.bam"
	params:
		fastp_parameters = return_parsed_extra_params(config['fastp_parameters'])
	run:
		shell("print {input}")
		shell("touch test.txt")
		shell("print {params.extra_star_parameters}")