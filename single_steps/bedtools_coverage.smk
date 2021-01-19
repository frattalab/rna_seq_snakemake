configfile: "config/single_steps_config.yaml"
import os


'''
Quick and dirty rule to output BAM mean coverage over intervals from BED file with bedtools
Using bedtools coverage, consult docs for how to modify command line to see what else it can do and chop and change at this
As I've currently set up:
1. Reports per-base strand-specific coverage on each nucleotide in each interval in BED (or GFF) file
2. Additionally computes mean & median (separately) coverage for each interval (bedtools groupby)
'''

#########---------------
## Input parameters
#########---------------
options_dict = config["bedtools_coverage"]
config_path = "config/single_steps_config.yaml" # this is for copying later

bam_dir = options_dict["bam_dir"]

# Pick out bams from bam_dir with this bam_suffix
# Also used to strip from file to get sample name
# e.g. Cont-B_S2.pass2Aligned.sortedByCoord.out.bam becomes Cont-B_S2
bam_suffix = options_dict["bam_suffix"]

bed_path = options_dict["bed_path"]

#Make sure trailing slash
output_dir = options_dict["output_dir"]

#########---------------
## BEDTOOLS PARAMETERS
#########---------------

# If don't want to use these options, assign them to "" (unless told otherwise)

# -s - only report hits overlapping on the same strand
# -S - only report hits overlapping on opposite strands
strandedness = options_dict["strandedness"]

# Report the depth at each position in each interval
depth = options_dict["depth"]
# what operations to summarise coverage in each intervals from bed_path? All below are valid strings
# sum, count, count_distinct, min, max, last
# median, mode, antimode, stdev, sstdev, collapse
# distinct, concat, freqasc, freqdesc, first, last
# If just want per-base coverages, assign to an empty list with []
operations = options_dict["operations"]

#1-based - which column should operation be performed on? (with depth = -d & a 6 col BED, this is the 8th column)
operation_column = options_dict["operations_column"]

# No need to change these unless necessary

sorted = options_dict["sorted"]
bedtools_path = options_dict["bedtools_path"]
samtools_path = options_dict["samtools_path"]
########-----------------

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(bam_dir) if f.endswith(bam_suffix)]

print(SAMPLES)

if not os._exists(output_dir):
    os.system("mkdir -p {0}".format(output_dir))


########-----------------
localrules: all, copy_config

rule all:
    input:
        expand(output_dir + "{sample}.coverage.per_base.tsv", sample = SAMPLES),
        expand(output_dir + "{sample}.coverage.{operation}.tsv" if len(operations) > 0 else [], sample = SAMPLES, operation = operations)
        os.path.join(output_dir, config_path)

# Get order of chromosome reference names from header of BAM file
# Enables sorting BED for each sample, so can use -sorted option
# Error otherwise thrown (similar to https://github.com/arq5x/bedtools/issues/109)
# -sorted option is less memory intensive on large files
# https://www.biostars.org/p/70795/ - finswimmer answer (minus a unmatched quote)

rule bam_chrom_order:
    input:
        os.path.join(bam_dir, "{sample}" + bam_suffix)

    output:
        temp(output_dir + "{sample}.genome.txt")

    params:
        samtools_path

    group: "group1"

    shell:
        """
        {params} view -H {input} | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > {output}
        """

rule bedtools_coverage:
    input:
        bam = os.path.join(bam_dir, "{sample}" + bam_suffix),
        bed = bed_path,
        chr_order = output_dir + "{sample}.genome.txt"

    output:
        output_dir + "{sample}.coverage.per_base.tsv"

    params:
        path = bedtools_path,
        strand = strandedness,
        per_base = depth,
        sorted = sorted,
        genome = " ".join(["-g", os.path.join(output_dir + "{sample}.genome.txt")]) if sorted == "-sorted" else ""

    group: "group1"

    shell:
        """
        {params.path} sort -i {input.bed} -g {input.chr_order} | \
        {params.path} coverage \
        -a stdin \
        -b {input.bam} \
        {params.strand} \
        {params.per_base} \
        {params.genome} \
        {params.sorted} > {output}
        """

rule bedtools_groupby:
    input:
        output_dir + "{sample}.coverage.per_base.tsv",

    output:
        output_dir + "{sample}.coverage.{operation}.tsv"

    params:
        path = bedtools_path,
        op_col = operation_column

    group: "group1"

    shell:
        """
        {params.path} groupby -i {input} \
        -g 1,2,3 \
        -c {params.op_col} \
        -o {wildcards.operation} > {output}
        """


rule copy_config:
    input:
        conf = config_path,
        expand(output_dir + "{sample}.coverage.per_base.tsv", sample = SAMPLES),
        expand(output_dir + "{sample}.coverage.{operation}.tsv" if len(operations) > 0 else [], sample = SAMPLES, operation = operations)

    output:
        os.path.join(output_dir, config_path)

    shell:
        """
        cp {input.conf} {output}
        """
