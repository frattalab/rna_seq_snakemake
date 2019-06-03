
configfile: "config/config.yaml"
include: "helpers.py"



#make sure the output folder for fastqc exists before running anything
os.system("mkdir -p {0}".format(config["fastp_trimmed_output_folder"]))
#read in a samples table
SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
#so I want this rule to be run ONCE for every fast1, so the wild cards I'm giving are the 'name' of the fastq of the first read
FASTQ_NAME = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in SAMPLES['fast1'].tolist()]

UNITS = SAMPLES['unit'].tolist()

if config['end_type'] == "pe":
	ruleorder: fastp_trimming_pe > fastp_trimming_se
else:
	ruleorder: fastp_trimming_se > fastp_trimming_pe
#here i'm defining the final output to be all the of the unit/fast1 trimmed files, this might actually become problematic for 
#paired end since I'm not also specifying 
rule all_trimmed:
	input: 
		expand(config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq.gz",zip, unit = UNITS, fastq_name=FASTQ_NAME),
		expand(config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.json",zip, unit = UNITS, fastq_name=FASTQ_NAME),
		expand(config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.html",zip, unit = UNITS, fastq_name=FASTQ_NAME)

#here i'm defining the final output to be all the of the unit/fast1 trimmed files
rule fastp_trimming_pe:
	input:
		fastq_file = lambda wildcards: return_fastq(wildcards.fastq_name,wildcards.unit,first_pair = True),
		fastq_file2 = lambda wildcards: return_fastq(wildcards.fastq_name,wildcards.unit,first_pair = False)
	output:
		out_fastqc = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq.gz",
		fastpjson = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.json",
		fastphtml = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.html",
	params:
		fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
		fastpjson = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.json",
		fastphtml = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.html",
		out_fastqc2 = lambda wildcards: return_fastq2_name(wildcards.fastq_name,wildcards.unit),

	run:
		shell("{config[fastp_path]} --in1 {input.fastq_file} --in2 {input.fastq_file2} --out1 {output.out_fastqc} --out2  {params.out_fastqc2} --json {output.fastpjson} --html {output.fastphtml} {params.fastp_parameters}")

rule fastp_trimming_se:
	input:
	#get the value in the fast1 column
		fastq_file = lambda wildcards: return_fastq(wildcards.fastq_name,wildcards.unit,first_pair = True)
	output:
		out_fastqc = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq.gz",
		fastpjson = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.json",
		fastphtml = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.html"
	params:
		fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
 		fastpjson = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.json",
		fastphtml = config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_fastp.html"
	run:
		shell("{config[fastp_path]} -i {input.fastq_file} -o {output.out_fastqc} --json {output.fastpjson} --html {output.fastphtml} {params.fastp_parameters}")
		