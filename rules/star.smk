
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "rules/helpers.py"

#make sure the output folder for STAR exists before running anything
os.system("mkdir -p {0}".format(config["star_output_folder"]))

def get_fastq_files_by_sample():
	return SAMPLES["samples"][wildcards.sample]

rule all_star:
	input:
		expand(config['star_output'] + "{sample_name}", fastq_name=FASTQ_NAME)

rule run_star:
	input:
		get_fastq_files_by_sample
	output:
		"star/{sample}-{unit}/Aligned.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab",
	parameters:
		extra_star_parameters = return_parsed_extra_params(config['fastp_parameters']),
		input_files = return_trimmed_fastqs(wildcards.unit, wildcards.fastconfig['end_type'])
	# threads: 24
	# ${starexec} 
	# --readFilesIn $f1_total $f2_total 
	# --readFilesCommand zcat 
	# --genomeLoad ${memorymode} 
	# --genomeDir /SAN/vyplab/HuRNASeq/reference_datasets/STAR/mouse/GrCh38.96/
	# --runThreadN  4 $chimeraCommand 
	# --outFileNamePrefix ${SCRATCH_DIR}/${sample} 
	# --outSAMtype $STARoutput $twopass 
	# --outSAMunmapped Within 
	# --outSAMheaderHD ID:${sample} PL:Illumina
	 STAR --genomeDir Genome1/ --genomeLoad LoadAndKeep --readFilesIn SampleTest_R1_trimmed.1M.fastq.gz SampleTest_R2_trimmed.1M.fastq.gz --readFilesCommand zcat

	run:
		cmd = ["{config[star_path]} \
		--readFilesIn {params.trimmed_fastqs} \
		--readFilesCommand zcat \
		--runThreadN {cluster_config[threads} \
		--outSAMtype BAM SortedByCoordinate"]
		print(cmd)
		shell(cmd)



/cluster/project8/vyp/vincent/Software/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR --readFilesIn /SAN/vyplab/TDP43_RNA/4SU_mouse/fastp_trimmed/190313_NS500195_0484_AH55VFBGXB/F1A_S21_R1_001_trimmed.fastq.gz --readFilesCommand zcat --runThreadN 4 --genomeDir /SAN/vyplab/HuRNASeq/reference_datasets/STAR/Mouse/GRCm38.96/star_indices_overhang100 --outSAMattributes MD NH --alignEndsType EndToEnd --outFileNamePrefix /SAN/vyplab/TDP43_RNA/4SU_mouse/star_test/190313_NS500195_0484_AH55VFBGXB 