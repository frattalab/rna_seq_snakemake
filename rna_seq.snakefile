import pandas as pd
import os
import subprocess

import subprocess

configfile: "config/config.yaml"
include: "rules/helpers.py"


#zip them into a directory to make getting the location easier
fq, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])

SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES.replace(np.nan, '', regex=True)
FASTQ_NAME = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in SAMPLES['fast1'].tolist()]


rule all:
	input: 
		expand(config["fastp_trimmed_output_folder"] + "{unit}/{fastq_name}_trimmed.fastq.gz",zip, unit = UNITS, fastq_name=FASTQ_NAME)


include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/fastp.smk"