
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "helpers.py"

#make sure the output folder for STAR exists before running anything
os.system("mkdir -p {0}".format(config["star_output_folder"]))

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
UNITS = SAMPLES['unit'].tolist()

rule all_sam:
	input:
		expand(config['star_output_folder'] + "{name}/{name}.Sorted.Merged.Aligned.out.bam",name = SAMPLE_NAMES)

rule samtools:
	input:
		expand(config['star_output_folder'] + "{name}/{unit}/{unit}_Aligned.out.bam",zip, unit = UNITS,name = SAMPLE_NAMES)
	output:
		config['star_output_folder'] + "{name}/{name}.Sorted.Merged.Aligned.out.bam"
	params:
		#collected_bams_by_sample = return_all_sample_aligned_bams(wildcards.name)

