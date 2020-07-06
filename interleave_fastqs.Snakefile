import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "rules/helpers.py"

FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

rule all:
    input:
        expand("{interleaved_outdir}{name}_interleaved.fastq.gz", interleaved_outdir = os.path.join(config['interleave_master_output_folder'],''),name = SAMPLE_NAMES)


include: "rules/fastp.smk"
include: "rules/interleave_fastqs.smk"
