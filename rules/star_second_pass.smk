import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml" 
include: "helpers.py"

# os.system("mkdir -p {0}".format(config["star_output_folder"]))

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

ALL_OUTS = [config['star_output_folder'] + s + "/" + s + ".SJ.out.tab" for s in list(set(SAMPLES.sample_name))]
# UNITS = SAMPLES['unit'].tolist()
# #this function uses the text file located in the config folder "star_genomes_species.csv" and 
# #the config file species parameter to 
# #give the correct genome for the species
GENOME_DIR = get_genome_directory(config['species'])

rule all_star_second:
	input:
		expand(config['star_output_folder'] + "{name}/{name}_2pass.Aligned.out.bam",name = SAMPLE_NAMES)

rule filter_and_copy_splice_junctions:
	input:
	#input is all the splice tables
		config['star_output_folder'] + "{name}/{name}.SJ.out.tab"
	output:
		config['star_output_folder'] + "{name}/filtered_combined.SJ.out.tab"	
	#make a quick thing called chrm that has the standard chroms, loop through all the sj.out files
	#make sure the first column is in the standard chroms, make sure that it's not a 'non-canoical splice junction(5th column != 0)
	#and the junction is supported by at least one unique mapper, if that's all true, then write the line.
	run:
		chrm = [str(x) for x in range(1,23)] + ["X","Y"]
		output_name = config['star_output_folder'] + wildcards.name + "/filtered_combined.SJ.out.tab"
		with open(output_name, 'w') as outfile:
		    for fname in ALL_OUTS:
		        with open(fname) as infile:
		            for line in infile:
		                split_line = line.split("\t")
		                if (split_line[0] in chrm and \
		                    split_line[4] != "0" and \
		                    int(split_line[6]) > 10):
		                    outfile.write(line)

rule run_star_second_pass_pe:
	input:
		first_pass = expand(config['star_output_folder'] + "{name}/{name}.Aligned.out.bam",name = SAMPLE_NAMES),
		one = lambda wildcards: get_trimmed(wildcards.name)[0],
		two = lambda wildcards: get_trimmed(wildcards.name)[1]
	output:
		config['star_output_folder'] + "{name}/{name}_2pass.Aligned.out.bam"
	threads: 4
	params:
		extra_star_parameters = return_parsed_extra_params(config['extra_star_parameters']),
		genomeDir = GENOME_DIR,
		outTmpDir = os.path.join(config['star_output_folder'] + "{name}/_tmpdir"),
		outputPrefix = os.path.join(config['star_output_folder'] + "{name}/{name}_2pass."),
		#taking the input files and putting them into a comma separated list
		one = lambda wildcards: ','.join(get_trimmed(wildcards.name)[0]),
		two = lambda wildcards: ','.join(get_trimmed(wildcards.name)[1]),
		splice_junctions_tabs = config['star_output_folder'] + "{name}/filtered_combined.SJ.out.tab",
		sjdb = str(int(file_len(config['star_output_folder'] + "{name}/filtered_combined.SJ.out.tab") + 1))
	shell: 
		"""
		rm -f -r {params.outTmpDir}
		{config[star_path]} --genomeDir {params.genomeDir} \
		--readFilesIn {params.one} \
		--outFileNamePrefix {params.outputPrefix} \
		--readFilesCommand zcat --runThreadN {threads} \
		{params.extra_star_parameters} \
		--outTmpDir {params.outTmpDir} \
		--sjdbFileChrStartEnd {params.splice_junctions_tabs} \
		--limitSjdbInsertNsj {params.sjdb}
		""" 

# rule run_star_second_pass_se:
# 	input:
# 		first_pass = expand(config['star_output_folder'] + "{name}/{name}.Aligned.out.bam",name = SAMPLE_NAMES),
# 		one = lambda wildcards: get_trimmed(wildcards.name)[0]
# 	output:
# 		config['star_output_folder'] + "{name}/{name}_2pass.Aligned.out.bam"
# 	params:
# 		extra_star_parameters = return_parsed_extra_params(config['extra_star_parameters']),
# 		genomeDir = GENOME_DIR,
# 		outTmpDir = os.path.join(config['star_output_folder'] + "{name}/_tmpdir"),
# 		outputPrefix = os.path.join(config['star_output_folder'] + "{name}/{name}_2pass."),
# 		#taking the input files and putting them into a comma separated list
# 		one = lambda wildcards: ','.join(get_trimmed(wildcards.name)[0]),
# 		#a whole list of all the sjdb files produced in the first run
# 		splice_junctions_tabs = "--sjdbFileChrStartEnd " + return_first_passing_splice_junctions_list()
# 	shell: 
# 		"""
# 		{config[star_path]} --genomeDir {params.genomeDir} \
# 		--readFilesIn {params.one} \
# 		--outFileNamePrefix {params.outputPrefix} \
# 		--readFilesCommand zcat --runThreadN {threads} \
# 		{params.extra_star_parameters} \
# 		--outTmpDir {params.outTmpDir} \
# 		{params.splice_junctions_tabs}
# 		""" 

