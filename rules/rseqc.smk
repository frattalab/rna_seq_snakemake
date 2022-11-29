import os
import pandas as pd

configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"



'''
Want to run:
inner_distance
read_distribution
gene_body_coverage
infer_experiment
junction_saturation

These will be collated by MultiQC

RSeQC requires transcript/gene models in BED12 format
For human and mouse, these will be specified in the refence_files_species.csv
There should be a single step to generate these if need to recreate (single_steps/gtf_to_bed12.smk)
'''

### Input parameters

STAR_OUTDIR = get_output_dir(config["project_top_level"], config['star_output_folder'])
RSEQC_OUTDIR = get_output_dir(config["project_top_level"], config["rseqc_output_folder"])
SPECIES = config["species"]

#make sure the output folder for rseqc exists before running anything
os.system("mkdir -p {0}".format(RSEQC_OUTDIR))

# RSeQC wants annotation in BED12 format. This can be pulled from refence_files_species
BED12 = get_bed12(SPECIES)


SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()


rule all_rseqc:
    input:
        expand(RSEQC_OUTDIR + "{sample}.geneBodyCoverage.txt", sample = SAMPLE_NAMES), #gene_body_coverage
        expand(RSEQC_OUTDIR + "{sample}.infer_experiment.txt", sample = SAMPLE_NAMES), #infer_experiment
        expand(RSEQC_OUTDIR + "{sample}.inner_distance_freq.txt", sample = SAMPLE_NAMES), # inner_distance
        expand(RSEQC_OUTDIR + "{sample}.junctionSaturation_plot.r", sample = SAMPLE_NAMES), # junction_saturation
        expand(RSEQC_OUTDIR + "{sample}.read_distribution.txt", sample = SAMPLE_NAMES) # read_distribution


##


rule gene_body_coverage:
    input:
        bam = STAR_OUTDIR + "{sample}.Aligned.sorted.out.bam",
        idx = STAR_OUTDIR + "{sample}.Aligned.sorted.out.bam.bai"

    output:
        os.path.join(RSEQC_OUTDIR, "{sample}.geneBodyCoverage.txt")

    params:
        # samples = lambda wildcards, input: ",".join(input.bams),
        annotation = BED12,
        prefix = os.path.join(RSEQC_OUTDIR, "{sample}")

    conda:
        "../env/align.yaml"

    shell:
        """
        geneBody_coverage.py \
        -i {input.bam} \
        -r {params.annotation} \
        -l 100 \
        -o {params.prefix}
        """


rule infer_experiment:
    input:
        bam = os.path.join(STAR_OUTDIR, "{sample}.Aligned.sorted.out.bam"),
        idx = os.path.join(STAR_OUTDIR, "{sample}.Aligned.sorted.out.bam")

    output:
        os.path.join(RSEQC_OUTDIR, "{sample}.infer_experiment.txt")

    params:
        annotation = BED12,
        sample_size = 200000, # Number of reads to sample from BAM
        min_qual = 30 # min map qual to be considered 'uniquely mapped'

    conda:
        "../env/align.yaml"

    shell:
        """
        infer_experiment.py \
        -i {input.bam} \
        -r {params.annotation} \
        -s {params.sample_size} \
        -q {params.min_qual} > {output}
        """

rule inner_distance_freq:
    input:
        bam = os.path.join(STAR_OUTDIR, "{sample}.Aligned.sorted.out.bam"),
        idx = os.path.join(STAR_OUTDIR, "{sample}.Aligned.sorted.out.bam")

    output:
        os.path.join(RSEQC_OUTDIR, "{sample}.inner_distance_freq.txt")

    params:
        annotation = BED12,
        prefix = os.path.join(RSEQC_OUTDIR, "{sample}")

    conda:
        "../env/align.yaml"

    shell:
        """
        inner_distance.py \
        -i {input.bam} \
        -r {params.annotation} \
        -o {params.prefix}
        """

rule junction_saturation:
    input:
        bam = os.path.join(STAR_OUTDIR, "{sample}.Aligned.sorted.out.bam"),
        idx = os.path.join(STAR_OUTDIR, "{sample}.Aligned.sorted.out.bam")

    output:
        os.path.join(RSEQC_OUTDIR, "{sample}.junctionSaturation_plot.r")

    params:
        annotation = BED12,
        prefix = os.path.join(RSEQC_OUTDIR, "{sample}")

    conda:
        "../env/align.yaml"

    shell:
        """
        junction_saturation.py \
        -i {input.bam} \
        -r {params.annotation} \
        -o {params.prefix}
        """

rule read_distribution:
    input:
        bam = os.path.join(STAR_OUTDIR, "{sample}.Aligned.sorted.out.bam"),
        idx = os.path.join(STAR_OUTDIR, "{sample}.Aligned.sorted.out.bam")

    output:
        os.path.join(RSEQC_OUTDIR, "{sample}.read_distribution.txt")

    params:
        annotation = BED12

    conda:
        "../env/align.yaml"

    shell:
        """
        read_distribution.py \
        -i {input.bam} \
        -r {params.annotation} > {output}
        """
