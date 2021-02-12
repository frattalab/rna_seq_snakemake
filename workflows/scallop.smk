import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "../rules/helpers.py"

include: "../rules/scallop_run_sample.smk"
#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name']))

print(SAMPLE_NAMES)

SPECIES = config["species"]
GTF = get_gtf(SPECIES)

#make sure the output folder for STAR exists before running anything
star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])

rule all:
    input:
        expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS),
        expand(interleaved_outdir + "{name}_interleaved.fastq.gz", name = SAMPLE_NAMES)
