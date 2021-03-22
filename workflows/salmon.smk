import pandas as pd
import os
import subprocess

configfile: "config/config.yaml"
include: "../rules/helpers.py"



#GENOME_DIR = os.path.join(config['STAR_indices'],config['species'],SPECIES_VERSION,"star_indices_overhang" + str(config['readLen']))
INDEX_DIR = config["salmon_indices"]
SPECIES = config["species"]
SPECIES_VERSION = get_species_version(SPECIES)
ANNOTATION_VERSION = get_annotation_version(SPECIES)
DECOY_TYPE = config["salmon_index_type"]
KMER_SIZE = config["salmon_index_kmer_size"]

TXOME_DIR = salmon_target_index(INDEX_DIR, SPECIES, SPECIES_VERSION, DECOY_TYPE, ANNOTATION_VERSION, KMER_SIZE)

#zip them into a directory to make getting the location easier
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])
SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

# This is to help multiqc know which files to track
workflow_str = "salmon"


#Construct all output directory strings needed for workflow definition
fastqc_outdir = get_output_dir(config["project_top_level"], config["fastqc_output_folder"])
salmon_outdir = get_output_dir(config["project_top_level"], config["salmon_output_folder"])

include: "../rules/fastp.smk"
include: "../rules/fastqc.smk"
include: "../rules/generate_salmon_index.smk"
include: "../rules/salmon_quant.smk"
include: "../rules/multiqc.smk"



rule all:
    input:
        expand(salmon_outdir + "{sample}/" + "quant.sf", sample = SAMPLE_NAMES),
        os.path.join(multiqc_output_folder, "multiqc_report.html")
        #os.path.join(TXOME_DIR, "seq.bin"),
        #os.path.join(TXOME_DIR, "pos.bin")
        #expand(fastqc_outdir + "{unit}/{fastq_name}_fastqc.html",zip, fastq_name=FASTQ_NAME, unit=UNITS)
    shadow: "minimal"
