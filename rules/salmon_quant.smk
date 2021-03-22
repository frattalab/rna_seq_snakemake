import os
import pandas as pd
import numpy as np

configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"


# RULE ORDER DIRECTIVE
# if paired end, use the paired end rule to run, if single end use the single end rule to run
if config['end_type'] == "pe":
    ruleorder: salmon_quant_pe > salmon_quant_se
else:
    ruleorder: salmon_quant_se > salmon_quant_pe


## Define target salmon transcriptome index

INDEX_DIR = config["salmon_indices"]
SPECIES = config["species"]
SPECIES_VERSION = get_species_version(SPECIES)
ANNOTATION_VERSION = get_annotation_version(SPECIES)
DECOY_TYPE = config["salmon_index_type"]
KMER_SIZE = config["salmon_index_kmer_size"]

# dir for Index generated for each species's genome assembly, each type of decoy sequence (homologous only or whole genome), provided txome annotation and kmer size
if config['custom_salmon_index'] == False:
    TXOME_DIR = salmon_target_index(INDEX_DIR, SPECIES, SPECIES_VERSION, DECOY_TYPE, ANNOTATION_VERSION, KMER_SIZE)
else:
    TXOME_DIR = config['custom_salmon_index']


#####
## Define input/output folders, sample names
#####

FASTQ_DIR = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

# Master output dir where salmon quant files are stored
OUTPUT_DIR = get_output_dir(config["project_top_level"], config["salmon_output_folder"])

if not os.path.exists(OUTPUT_DIR):
    os.system("mkdir -p {}".format(OUTPUT_DIR))

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

# # a top level folder where the bams reside
# project_dir = "/SAN/vyplab/alb_projects/data/4su_full_ward_tdp_kd_ipsc/"
# out_spot = "salmon/"
# fastq_spot = "merged_fastqs/"
#
# salmon_strand_info = "-l ISF"
# end_type = "pe"
#
#
# # =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------
#
# output_dir = os.path.join(project_dir,out_spot)
# fastq_dir = os.path.join(project_dir,fastq_spot)
#
# SAMPLES, = glob_wildcards(fastq_dir + "{sample}_1.merged.fastq.gz")
# print(SAMPLES)

######
## Main pipeline definition
######


rule salmon_quant_all:
    input:
        expand(OUTPUT_DIR + "{sample}/" + "quant.sf", sample = SAMPLE_NAMES)


rule salmon_quant_pe:
    input:
        fast1 = FASTQ_DIR  + "{sample}_1.merged.fastq.gz",
        fast2 = FASTQ_DIR  + "{sample}_2.merged.fastq.gz",
        index = os.path.join(TXOME_DIR, "seq.bin")

    output:
        os.path.join(OUTPUT_DIR, "{sample}", "quant.sf")

    params:
        salmon = config["salmon_path"],
        index_dir = TXOME_DIR,
        output_dir = os.path.join(OUTPUT_DIR, "{sample}"),
        libtype = get_salmon_strand(config["feature_counts_strand_info"]),
        gtf = get_gtf(SPECIES),
        extra_params = return_parsed_extra_params(config["extra_salmon_parameters"])

    # threads: 4

    shell:
        """
        {params.salmon} quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --geneMap {params.gtf} \
        {params.extra_params} \
        -o {params.output_dir} \
        """

rule salmon_quant_se:
    input:
        fast1 = merged_outdir + "{name}_1.merged.fastq.gz",
        index = os.path.join(TXOME_DIR, "seq.bin")

    output:
        os.path.join(OUTPUT_DIR, "{sample}", "quant.sf")

    params:
        salmon = config["salmon_path"],
        index_dir = TXOME_DIR,
        output_dir = os.path.join(OUTPUT_DIR, "{sample}"),
        libtype = get_salmon_strand(config["feature_counts_strand_info"]),
        gtf = get_gtf(SPECIES),
        extra_params = return_parsed_extra_params(config["extra_salmon_parameters"])

    shell:
        """
        {params.salmon} quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        -r {input.fast1} \
        --geneMap {params.gtf} \
        {params.extra_params} \
        -o {params.output_dir} \
        """
