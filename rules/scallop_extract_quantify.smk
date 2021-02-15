import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "helpers.py"

SPECIES = config["species"]
GENOME_FA = get_genome_fasta(SPECIES)

#make sure the output folder for STAR exists before running anything
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])
print(scallop_outdir)
txome_fa = get_transcriptome_fasta(SPECIES)

# 1 decoys (& merged txome + decoys FA ) file generated for each genome assembly + transcriptome annotation version
# Used as input to salmon index
INDEX_DIR = config["salmon_indices"]
DECOYS_DIR = os.path.join(INDEX_DIR, SPECIES, SPECIES_VERSION, "decoys", DECOY_TYPE, ANNOTATION_VERSION, "")
print(DECOYS_DIR)

rule extraction_quantification:
    input:
        os.path.join(scallop_outdir,"scallop_unique.fa"),
        os.path.join(scallop_outdir, "extended_transcriptome/seq.bin"),
        os.path.join(scallop_outdir, "extended_transcriptome/pos.bin")

rule get_cnda:
    input:
        os.path.join(scallop_outdir,"scallop_merged.gtf")
    output:
        os.path.join(scallop_outdir,"scallop_unique.fa")
    params:
        gffread = config['gffread']
    shell:
        "{params.gffread} {input} -g {GENOME_FA} -w {output}"

rule build_extended_cdna:
    input:
        os.path.join(scallop_outdir,"scallop_merged.gtf")
    output:
        os.path.join(scallop_outdir,"scallop_union.fa")
    params:
        gffread = config['gffread']
    shell:
        "cat unique.fa {txome_fa} > {output}"

rule salmon_index_extended:
    input:
        extended_fa = os.path.join(scallop_outdir,"scallop_union.fa"),
        decoys = os.path.join(DECOYS_DIR, "decoys.txt")

    output:
        os.path.join(scallop_outdir, "extended_transcriptome/seq.bin"),
        os.path.join(scallop_outdir, "extended_transcriptome/pos.bin")

    params:
        salmon = config["salmon_path"],
        k = KMER_SIZE,
        outdir = os.path.join(scallop_outdir, "extended_transcriptome"),
        gencode = "--gencode" if config["transcriptome_source"] == "gencode" else ""

    threads:
        4

    shell:
        """
        {params.salmon} index \
        -t {input.gentrome_fa} \
        -i {params.outdir} \
        --decoys {input.decoys} \
        -k {params.k} \
        {params.gencode} \
        -p {threads}
        """
