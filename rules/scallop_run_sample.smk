import pandas as pd
import os
import subprocess
import yaml
configfile: "config/config.yaml"
include: "helpers.py"
localrules: compose_gtf_list

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
print(scallop_outdir)
rule all_scallop:
    input:
        expand(scallop_outdir + '{sample}' + ".gtf", sample = SAMPLE_NAMES),
        expand(scallop_outdir + "gffall.{sample}.gtf.map",sample = SAMPLE_NAMES)


rule scallop_per_samp:
    input:
        bam_file = lambda wildcards: star_outdir + '{sample}' + config['bam_suffix']
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    output:
        os.path.join(scallop_outdir,'{sample}' + ".gtf")
    params:
        scallop_path = config['scallop_path'],
        scallop_out_folder = scallop_outdir,
        scallop_extra_config = return_parsed_extra_params(config['scallop_extra_parameters'])
    shell:
        """
        mkdir -p {params.scallop_out_folder}
        {params.scallop_path} -i {input.bam_file} -o {output} {params.scallop_extra_config}
        """

rule compare_reference:
    input:
        os.path.join(scallop_outdir,'{sample}' + ".gtf")
    output:
        os.path.join(scallop_outdir, "gffall.{sample}.gtf.map")
    params:
        ref_gtf = GTF,
        gffcompare = "/SAN/vyplab/alb_projects/tools/gffcompare-0.11.6.Linux_x86_64/gffcompare"
    shell:
        """
        {params.gffcompare} -o gffall -r {params.ref_gtf} {input}
        """

# rule compose_gtf_list:
#     input:
#         expand(scallop_outdir + '{sample}.gtf', sample=SAMPLE_NAMES)
#     output:
#         txt = os.path.join(scallop_outdir,"gtf_list.txt")
#     run:
#         with open(output.txt, 'w') as out:
#             print(*input, sep="\n", file=out)
# rule merge_scallop_gtfs:
#     input:
#         gtf_list = os.path.join(scallop_outdir,"gtf_list.txt")
#     output:
#         merged_gtf = os.path.join(scallop_outdir,"scallop_merged.gtf")
#     params:
#         gtfmerge = '/SAN/vyplab/alb_projects/tools/rnaseqtools-1.0.3/gtfmerge/gtfmerge'
#     shell:
#         """
#         {params.gtfmerge} union {input.gtf_list} {output.merged_gtf} -t 2 -n
#         """
