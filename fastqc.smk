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

def return_fastq_location(wildcards.fastq_name):
	#return the file location from the fastq name
	return(ORDER_DICT[wildcards.fastq_name])

#make sure the output folder for fastqc exists before running anything
os.system("mkdir -p {0}".format(config["fastqc_output_folder"]))


#fastq name is the sample_unit_fastqfile 
FASTQ_NAME, FILE_LOCATION = get_fastq_names(config["sampleCSVpath"])
ORDER_DICT = dict(zip(FASTQ_NAME, FILE_LOCATION))

rule all:
	input: 
		expand(config["fastqc_output_folder"] + "{fastq_name}_fastqc.html",fastq_name=FASTQ_NAME)

rule fastqc:
	input:
		fastq_file = lambda wildcards: return_fastq_location(wildcards.fastq_name)
	output:
		out_file = config["fastqc_output_folder"] + "{fastq_name}_fastqc.html"
	# params:
	# 	fastqc_call = config["fastqc_path"]
	shell:
		"{config[fastqc_path]} {input.fastq_file} -o {out_file}"