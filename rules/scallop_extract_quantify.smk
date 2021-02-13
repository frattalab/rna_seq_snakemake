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

rule extraction_quantification:
    input:
        os.path.join(scallop_outdir,"scallop_unique.fa")

rule get_cnda:
    input:
        os.path.join(scallop_outdir,"scallop_merged.gtf")
    output:
        os.path.join(scallop_outdir,"scallop_unique.fa")
    params:
        gffread = config['gffread']
    shell:
        "{params.gffread} {input} -g {GENOME_FA} -w {output}"
