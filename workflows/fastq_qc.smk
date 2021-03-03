import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "../rules/helpers.py"
SPECIES_VERSION = get_species_version(config['species'])
GENOME_DIR = os.path.join(config['STAR_indices'],config['species'],SPECIES_VERSION,"star_indices_overhang" + str(config['readLen']))
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])
#zip them into a directory to make getting the location easier
SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

#Construct all output directory strings needed for workflow definition
fastqc_outdir = get_output_dir(config["project_top_level"], config["fastqc_output_folder"])


# This is to help multiqc know which files to track
workflow_str = "fastq_qc"


include: "../rules/fastp.smk"
include: "../rules/fastqc.smk"
include: "../rules/multiqc.smk"

rule all:
    input:
        expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)
    shadow: "minimal"
