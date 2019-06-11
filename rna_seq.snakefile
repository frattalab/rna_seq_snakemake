import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "rules/helpers.py"


#zip them into a directory to make getting the location easier
SAMPLES = pd.read_table(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

rule all:
	input:
		expand(config['feature_counts_output_folder'] + "{name}_featureCounts_results.txt", name = SAMPLE_NAMES)


# include: "rules/fastqc.smk"
# include: "rules/multiqc.smk"
#include: "rules/fastp.smk"
include: "rules/merge_fastqc.smk"
include: "rules/star.smk"
include: "rules/samtools.smk"
include: "rules/feature_counts.smk"