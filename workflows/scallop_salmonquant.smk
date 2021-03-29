import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "../rules/helpers.py"

include: "../rules/scallop_run_sample.smk"
include: "../rules/scallop_extract_quantify.smk"
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
        expand(scallop_outdir + '{sample}' + ".gtf", sample = SAMPLE_NAMES),
        expand(scallop_outdir + "gffall.{sample}.gtf.tmap",sample = SAMPLE_NAMES),
        expand(scallop_outdir + "{sample}.unique.gtf",sample = SAMPLE_NAMES),
        expand(scallop_outdir + "{sample}/" + "quant.sf", sample = SAMPLE_NAMES),
        os.path.join(scallop_outdir,"scallop_merged.gtf"),
        os.path.join(scallop_outdir,"scallop_unique.fa"),
        os.path.join(scallop_outdir, "extended_transcriptome/seq.bin"),
        os.path.join(scallop_outdir, "extended_transcriptome/pos.bin")
