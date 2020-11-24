# Quick and dirty snakemake rule to downsample a pair of fastq files
# I used this to generate smaller FASTQs for the example  sample table and configs
# so testing is nice and quick

# Uses bbmap to get pairs of fastq files with certain number of read pairs

import os
import pandas as pd
include: "../rules/helpers.py"

sample_tbl_path = "/SAN/vyplab/vyplab_reference_genomes/pipeline_test_data/ward_tiny_fastqs_sample_tbl.csv"
n_reads = 100000 #100k should be plenty...
output_dir = "/SAN/vyplab/vyplab_reference_genomes/pipeline_test_data/ward_tiny_fastqs/" #make sure trailing slash
end_type = "pe" #"pe" or something not "pe" if single-ended data


## Below shouldn't change

bbmap_path = "/SAN/vyplab/alb_projects/tools/BBMap_38_86/bbmap/"
SAMPLES = pd.read_csv(sample_tbl_path, sep=",")
UNITS = SAMPLES['unit'].tolist()
SAMPLE_NAMES = SAMPLES['sample_name'].tolist()


#print(SAMPLES)
#print("\n----THIS IS SAMPLE_NAMES----\n")
#print(SAMPLE_NAMES)
#print(SAMPLES.columns.values)
#print(SAMPLES.loc[(SAMPLES['sample_name'] == 'Cont-A_S1') & (SAMPLES['unit'] == 'Cont-A_S1')]["fast1"].values[0])


os.system("mkdir -p {0}".format(output_dir))

wildcard_constraints:
            unit="|".join(UNITS),
            name="|".join(SAMPLE_NAMES)


rule all:
    input:
        expand(output_dir + "{unit}_{name}_R1_subsampled.fastq.gz",zip, unit = UNITS, name = SAMPLE_NAMES),
        expand(output_dir + "{unit}_{name}_R2_subsampled.fastq.gz" if end_type == "pe" else [], zip, unit = UNITS, name = SAMPLE_NAMES)


if end_type == "pe":
    rule subsample_fastqs:
        input:
            fq1 = lambda wildcards: SAMPLES.loc[(SAMPLES['sample_name'] == wildcards.name) & (SAMPLES['unit'] == wildcards.unit)]["fast1"].values[0],
            fq2 = lambda wildcards: SAMPLES.loc[(SAMPLES['sample_name'] == wildcards.name) & (SAMPLES['unit'] == wildcards.unit)]["fast2"].values[0]

        output:
            out_fq1 = output_dir + "{unit}_{name}_R1_subsampled.fastq.gz",
            out_fq2 = output_dir + "{unit}_{name}_R2_subsampled.fastq.gz"

        params:
            script = bbmap_path + "reformat.sh",
            n_pairs = n_reads

        threads:
            4

        shell:
            """
            bash {params.script} in={input.fq1} in2={input.fq2} \
            out={output.out_fq1} out2={output.out_fq2} \
            srt={params.n_pairs}
            """

else:
    rule subsample_fastqs:
        input:
            fq1 = lambda wildcards: SAMPLES.loc[(SAMPLES['sample_name'] == wildcards.name) & (SAMPLES['unit'] == wildcards.unit)]["fast1"].values[0]

        output:
            out_fq1 = output_dir + "{unit}_{name}_R1_subsampled.fastq.gz",

        params:
            script = bbmap_path + "reformat.sh",
            n = n_reads

        threads:
            4

        shell:
            """
            bash {params.script} in={input.fq1} \
            out={output.out_fq1} \
            srt={params.n}
            """
