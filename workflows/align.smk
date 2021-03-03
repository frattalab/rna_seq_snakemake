import pandas as pd
import os
import subprocess

# This is to help multiqc know which files to track
workflow_str = "align"

configfile: "config/config.yaml"
include: "../rules/helpers.py"
include: "../rules/fastp.smk"
include: "../rules/fastqc.smk"
include: "../rules/generate_star_index.smk"
include: "../rules/star.smk"
include: "../rules/feature_counts.smk"
include: "../rules/tpmcalculator.smk"
include: "../rules/rseqc.smk"
include: "../rules/multiqc.smk"


SPECIES_VERSION = get_species_version(config['species'])
GENOME_DIR = os.path.join(config['STAR_indices'],config['species'],SPECIES_VERSION,"star_indices_overhang" + str(config['readLen']))
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])
#zip them into a directory to make getting the location easier
SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()



#Construct all output directory strings needed for workflow definition
fastqc_outdir = get_output_dir(config["project_top_level"], config["fastqc_output_folder"])
interleaved_outdir = get_output_dir(config['project_top_level'], config['interleave_master_output_folder'])
star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])
feature_counts_outdir = get_output_dir(config["project_top_level"], config["feature_counts_output_folder"])
tpm_outdir = get_output_dir(config['project_top_level'], config['tpmcalculator_output_folder'])
multiqc_output_folder = os.path.join(get_output_dir(config["project_top_level"], config["multiqc_output_folder"]), workflow_str, "")




rule all:
    input:
        GENOME_DIR + "/SA",
        expand(feature_counts_outdir + "{name}_featureCounts_results.txt", name = SAMPLE_NAMES),
        expand(tpm_outdir + "{name}" + suffix + "_genes.out", name = SAMPLE_NAMES),
        expand(star_outdir + "{name}.Aligned.sorted.out.bam", name = SAMPLE_NAMES),
        expand(star_outdir + "{name}.Aligned.sorted.out.bam.bai", name = SAMPLE_NAMES),
        os.path.join(multiqc_output_folder, "multiqc_report.html")
        # expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)
        #expand(config['project_top_level'] + "multiqc_report.html")
    shadow: "minimal"
