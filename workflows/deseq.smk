import os
import pandas as pd
import numpy as np
import yaml

configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "../rules/helpers.py"

with open("config/config.yaml", "r") as stream:
    config = yaml.safe_load(stream)

BASES, CONTRASTS = return_bases_and_contrasts('config/DESeq2comparisons.yaml')
print(BASES)
print(CONTRASTS)

DESEQ_DIR = os.path.join(config['project_top_level'] , config['DESeq2_output'])
SALMON_DIR = os.path.join(config['project_top_level'] , config['salmon_output_folder'])

gtf = get_gtf(config['species'])
tx_gene = gtf.rstrip('\annotation.gtf') + '.tx2gene.csv'

rule all:
    input:
        expand(os.path.join(DESEQ_DIR, "{bse}_{contrast}" + ".DESEQ2_results.csv"), zip, bse = BASES, contrast = CONTRASTS),
        tx_gene = gtf.rstrip('\annotation.gtf') + '.tx2gene.csv'

rule write_tx2gene:
     input:
        gtf = get_gtf(config['species'])

     output:
        tx_gene = gtf.rstrip('\annotation.gtf') + '.tx2gene.csv'


     shell:
        """
        python scripts/write_tx_gene.py \
        --gtf {input.gtf} \
        --output {output.tx_gene} \
        """

rule run_deseq:
    input:
        base_group = lambda wildcards: salmon_files_from_contrast(wildcards.bse),
        contrast_group = lambda wildcards: salmon_files_from_contrast(wildcards.contrast)

    wildcard_constraints:
        bse="|".join(BASES),
        contrast="|".join(CONTRASTS)

    output:
        os.path.join(DESEQ_DIR, "{bse}_{contrast}" + ".DESEQ2_results.csv")

    params:
        samplesheet = config['sampleCSVpath'],
        salmon_dir = SALMON_DIR.rstrip('\/'),
        deseq_dir = DESEQ_DIR.rstrip('\/'),
        tx2gene = tx_gene

    conda:
        "../env/deseq2.yaml"

    shell:
        """
        Rscript scripts/quick_salmon_deseq_command_line.R \
        --salmon_quant {params.salmon_dir} \
        --samplesheet {params.samplesheet} \
        --tx2gene {params.tx2gene} \
        --outputdir {params.deseq_dir} \
        """
