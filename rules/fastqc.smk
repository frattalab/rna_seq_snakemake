import re
import pandas as pd
import os
import numpy as np

configfile: "config/config.yaml"

def get_fastq_names(DATA):
    samples = pd.read_table(DATA, sep = ",")
    #take all the fasqs and combine them to a list remove any nas
    fastq_list = samples.fast1.dropna().tolist() + samples.fast2.dropna().tolist()
    #if there are any missing values pandas gets annoyed so replace nans with empty string
    samples = samples.replace(np.nan, '', regex=True)
    #get the associated sample and unit name for each fastq using and or operator to get either fast1 or fast2, works for both single-end and paired end this way
    unit_name = [samples.loc[(samples['fast1'] == fq)| (samples['fast2'] == fq)].unit.iloc[0] for fq in fastq_list]
    sample_name = [samples.loc[(samples['fast1'] == fq)| (samples['fast2'] == fq)].sample_name.iloc[0] for fq in fastq_list]
    #strip it down to just the name of the file bit
    stripped = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in fastq_list]
    #now combine these three lists to get a meaningfull and unique string for each fastq file
    fn = []
    for ind, value in enumerate(stripped):
        fn.append("_".join([sample_name[ind], unit_name[ind], value]))
    
    return(fn, fastq_list)

def return_fastq_location(wildcards):
	#return the file location from the fastq name
	return(ORDER_DICT[wildcards])

#make sure the output folder for fastqc exists before running anything
os.system("mkdir -p {0}".format(config["fastqc_output_folder"]))


#fastq name is the sample_unit_fastqfile 
FASTQ_NAME, FILE_LOCATION = get_fastq_names(config["sampleCSVpath"])
#zip them into a directory to make getting the location easier
ORDER_DICT = dict(zip(FASTQ_NAME, FILE_LOCATION))
#first rule is a general rule that specifies the final output of everything, here we have the expected
#output of the individual fastqc's and the multiqc html file. snakemake will check to see if these output files
#exist if they dont it will keep going
rule all:
	input: 
		expand(config["fastqc_output_folder"] + "{fastq_name}_fastqc.html",fastq_name=FASTQ_NAME)

rule fastqc:
	input:
		fastq_file = lambda wildcards: return_fastq_location(wildcards.fastq_name)
	output:
		out_fastqc = config["fastqc_output_folder"] + "{fastq_name}_fastqc.html",
		out_fastzip = config["fastqc_output_folder"] + "{fastq_name}_fastqc.zip"
	params:
		fq_name = lambda wildcards, input: re.sub(".gz","c", input.fastq_file.rpartition('/')[2])
	log:
        "logs/{fastq_name}_fastqc.log"
	shell:
		"""
		mkdir -p {config[fastqc_output_folder]}{wildcards.fastq_name}
		#{config[fastqc_path]} {input.fastq_file} -o {config[fastqc_output_folder]}{wildcards.fastq_name} 
		mv {config[fastqc_output_folder]}{wildcards.fastq_name}/{params.fq_name}_.html {output.out_fastqc} 
		mv {config[fastqc_output_folder]}{wildcards.fastq_name}/{params.fq_name}_.zip {output.out_fastzip}
		"""
