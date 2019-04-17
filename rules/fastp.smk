
configfile: "config/config.yaml"
include: "helpers.py"

#zip them into a directory to make getting the location easier
fq, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])

#make sure the output folder for fastqc exists before running anything
os.system("mkdir -p {0}".format(config["fastp_trimmed_output_folder"]))
#read in a samples table
SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES.replace(np.nan, '', regex=True)
#so I want this rule to be run ONCE for every fast1, so the wild cards I'm giving are the 'name' of the fastq of the first read
FASTQ_NAME = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in SAMPLES['fast1'].tolist()]
FASTQ_NAME_2 = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) + "_trimmed.fastq.gz" for strpd in SAMPLES['fast2'].tolist()]

#I'm adding on the fastq name to the samples table so that it makes it easier to select on in later steps
SAMPLES['fast1_name'] = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in SAMPLES['fast1'].tolist()]
#here i'm defining the final output to be all the of the unit/fast1 trimmed files, this might actually become problematic for 
#paired end since I'm not also specifying 
rule all_trimmed:
	input: 
		expand(config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq.gz",zip, unit = UNITS, fastq_name=FASTQ_NAME),
		expand(config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name2}",zip, unit = UNITS, fastq_name=FASTQ_NAME_2)

#here i'm defining the final output to be all the of the unit/fast1 trimmed files

rule fastp_trimming:
	input:
	#get the value in the fast1 column
		fastq_file = lambda wildcards: SAMPLES.loc[(SAMPLES['fast1_name'] == wildcards.fastq_name) & (SAMPLES['unit'] == wildcards.unit)]["fast1"].values[0]
	output:
		out_fastqc = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq.gz",
		out_fastqc2 = 
	params:
		fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
		fastq_file2 = lambda wildcards: SAMPLES.loc[(SAMPLES['fast1_name'] == wildcards.fastq_name) & (SAMPLES['unit'] == wildcards.unit)]["fast2"].values[0]
	run:
		if config["end_type"] == "se":
			shell("{config[fastp_path]} -i {input.fastq_file} -o {output.out_fastqc} {params.fastp_parameters}")
		if config["end_type"] == "pe":
			shell("{config[fastp_path]} -i {input.fastq_file} -in2 {params.fastq_file2} -o {config[fastp_trimmed_output_folder]}{output.out_fastqc} -o2 {params.out_fastqc2} {params.fastp_parameters}")
