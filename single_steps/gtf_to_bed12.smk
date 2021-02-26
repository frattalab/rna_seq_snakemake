import os
configfile: "config/single_steps_config.yaml"



'''
Convert a GTF entries to BED12 format
Uses UCSC tools. We have scripts on cluster
Edit single_steps_config with target gtf and out spot
'''

options = config["gtf_to_bed12"]

rule all:
    input:
        options["bed12"]


rule gtf_to_bed12:
    input:
        options["gtf"]

    output:
        options["bed12"]

    params:
        gtf2pred = os.path.join(options["ucsc_tools_dir"], "gtfToGenePred"),
        pred2bed = os.path.join(options["ucsc_tools_dir"], "genePredToBed"),
        gene_pred = options["temp_gene_pred"]

    shell:
        """
        {params.gtf2pred} {input} {params.gene_pred}
        {params.pred2bed} {params.gene_pred} {output}
        rm {params.gene_pred}
        """
