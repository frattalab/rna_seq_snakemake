import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "helpers.py"

SPECIES = config["species"]
GENOME_FA = get_genome_fasta(SPECIES)
SPECIES_VERSION = get_species_version(SPECIES)
INDEX_DIR = config["salmon_indices"]
ANNOTATION_VERSION = get_annotation_version(SPECIES)
KMER_SIZE = config["salmon_index_kmer_size"]
DECOYS_DIR = os.path.join(INDEX_DIR, SPECIES, SPECIES_VERSION, "decoys", "full", ANNOTATION_VERSION, "")
print(DECOYS_DIR)
FASTQ_DIR = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

#make sure the output folder for STAR exists before running anything
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])
print(scallop_outdir)
txome_fa = get_transcriptome_fasta(SPECIES)

# 1 decoys (& merged txome + decoys FA ) file generated for each genome assembly + transcriptome annotation version
# Used as input to salmon index


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
        os.path.join(scallop_outdir,"scallop_unique.fa")
    output:
        os.path.join(scallop_outdir,"scallop_union.fa")
    params:
        gffread = config['gffread']
    shell:
        "cat {input} {txome_fa} > {output}"

rule create_tab_delimited_gene_txt_reference:
    input:
        os.path.join(scallop_outdir,"scallop_merged.gtf")
    output:
        os.path.join(scallop_outdir,"scallop_union.fa")
    params:
        gffread = config['gffread']
    shell:
        """
        set +u;
        source activate salmon
        python3 {params.quick_script} --gtf {input.scallop}
        """

rule create_tab_delimited_gene_txt_scallop:
    input:
        os.path.join(scallop_outdir,"scallop_merged.gtf")
    output:
        os.path.join(scallop_outdir,"scallop_union.fa")
    params:
        gffread = config['gffread']
    shell:
        """
        set +u;
        source activate salmon
        python3 {params.quick_script} --gtf {input.scallop}
        """

rule generate_full_decoy:
    input:
        genome_fa = get_genome_fasta(SPECIES),
        txome_fa = os.path.join(scallop_outdir,"scallop_unique.fa"),
    output:
        gentrome_fa = os.path.join(scallop_outdir, "gentrome.fa"),
        decoys = os.path.join(scallop_outdir, "decoys.txt")
    params:
        outdir = DECOYS_DIR
    shell:
        """
        grep "^>" {input.genome_fa} | cut -d " " -f 1 > {output.decoys}
        sed -i.bak -e 's/>//g' {output.decoys}

        cat {input.txome_fa} {input.genome_fa} > {output.gentrome_fa}
        """
rule salmon_index_extended:
    input:
        extended_fa = os.path.join(scallop_outdir, "gentrome.fa"),
        decoys = os.path.join(scallop_outdir, "decoys.txt")
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
        -t {input.extended_fa} \
        -i {params.outdir} \
        --decoys {input.decoys} \
        -k {params.k} \
        {params.gencode} \
        -p {threads}
        """

rule salmon_quant:
    input:
        fast1 = FASTQ_DIR  + "{sample}_1.merged.fastq.gz",
        fast2 = FASTQ_DIR  + "{sample}_2.merged.fastq.gz",
        index = os.path.join(scallop_outdir, "extended_transcriptome/seq.bin")
    output:
        os.path.join(OUTPUT_DIR, "{sample}", "quant.sf")
    params:
        salmon = config["salmon_path"],
        index_dir = os.path.join(scallop_outdir, "extended_transcriptome/")",
        output_dir = os.path.join(OUTPUT_DIR, "{sample}"),
        libtype = get_salmon_strand(config["feature_counts_strand_info"]),
        gtf = get_gtf(SPECIES),
        extra_params = return_parsed_extra_params(config["extra_salmon_parameters"])
    threads: 4
    shell:
        """
        {params.salmon} quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --geneMap {params.gtf} \
        --threads {threads} \
        {params.extra_params} \
        -o {params.output_dir} \
        """
