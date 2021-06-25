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

GTF = get_gtf(SPECIES)
TAB_GTF = GTF.replace(".gtf",".tx_gene.tsv")
print(TAB_GTF)
#make sure the output folder for STAR exists before running anything
scallop_outdir = get_output_dir(config["project_top_level"], config['scallop_output'])
print(scallop_outdir)
txome_fa = get_transcriptome_fasta(SPECIES)


SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
# 1 decoys (& merged txome + decoys FA ) file generated for each genome assembly + transcriptome annotation version
# Used as input to salmon index


rule extraction_quantification:
    input:
        os.path.join(scallop_outdir,"scallop_unique.fa"),
        os.path.join(scallop_outdir, "extended_transcriptome/seq.bin"),
        os.path.join(scallop_outdir, "extended_transcriptome/pos.bin"),
        os.path.join(scallop_outdir,"scallop_ref.tx_gene.tsv"),
        expand(scallop_outdir + "{sample}/" + "quant.sf", sample = SAMPLE_NAMES),
        TAB_GTF


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

rule create_tab_delimited_tx_gene_reference:
    input:
        GTF
    output:
        TAB_GTF
    params:
        quick_script = "scripts/write_tx_gene.py"
    shell:
        """
        set +u;
        source activate salmon
        python3 {params.quick_script} --gtf {input} --output {output}
        """

rule create_tab_delimited_tx_gene_scallop:
    input:
        scallop_gtf = os.path.join(scallop_outdir,"scallop_merged.gtf")
    output:
        os.path.join(scallop_outdir,"scallop.tx_gene.tsv")
    params:
        quick_script = "scripts/write_tx_gene.py"
    shell:
        """
        set +u;
        source activate salmon
        python3 {params.quick_script} --gtf {input.scallop_gtf} --output {output}
        """

rule cat_tabs:
    input:
        scallop_tx_gene = os.path.join(scallop_outdir,"scallop.tx_gene.tsv"),
        ref_tx_gene = TAB_GTF
    output:
        os.path.join(scallop_outdir,"scallop_ref.tx_gene.tsv")
    shell:
        """
        cat {input.scallop_tx_gene} {input.ref_tx_gene} > {output}
        """

rule generate_full_decoy:
    input:
        genome_fa = get_genome_fasta(SPECIES),
        scallop_txome_fa = os.path.join(scallop_outdir,"scallop_unique.fa"),
        ref_txome_fa = get_transcriptome_fasta(SPECIES),

    output:
        gentrome_fa = os.path.join(scallop_outdir, "gentrome.fa"),
        decoys = os.path.join(scallop_outdir, "decoys.txt")
    params:
        outdir = DECOYS_DIR
    shell:
        """
        grep "^>" {input.genome_fa} | cut -d " " -f 1 > {output.decoys}
        sed -i.bak -e 's/>//g' {output.decoys}

        cat {input.scallop_txome_fa} {input.ref_txome_fa} {input.genome_fa} > {output.gentrome_fa}
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
        fast1 = lambda wildcards: get_processed_fastq(wildcards.sample, pair=1),
        fast2 = lambda wildcards: get_processed_fastq(wildcards.sample, pair=1),
        index = os.path.join(scallop_outdir, "extended_transcriptome/seq.bin"),
        scallop_ref = os.path.join(scallop_outdir,"scallop.tx_gene.tsv")
    output:
        os.path.join(scallop_outdir, "{sample}", "quant.sf")
    params:
        salmon = config["salmon_path"],
        index_dir = os.path.join(scallop_outdir, "extended_transcriptome/"),
        output_dir = os.path.join(scallop_outdir, "{sample}"),
        libtype = get_salmon_strand(config["feature_counts_strand_info"]),
        extra_params = return_parsed_extra_params(config["extra_salmon_parameters"])
    threads: 4
    shell:
        """
        {params.salmon} quant \
        --gcBias \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --geneMap {input.scallop_ref} \
        --threads {threads} \
        {params.extra_params} \
        -o {params.output_dir} \
        """
