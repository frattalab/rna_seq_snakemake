#configfile: "config.yaml"
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
bam_dir = "/SAN/vyplab/alb_projects/data/ward_bams/bam_files/"

# Pick out bams from bam_dir with this bam_suffix
# Also used to strip from file to get sample name
# e.g. Cont-B_S2.pass2Aligned.sortedByCoord.out.bam becomes Cont-B_S2
bam_suffix = ".pass2Aligned.sortedByCoord.out.bam"

bed_path = "/SAN/vyplab/sbs_projects/data/nanopore/IPSC_TDP43_KD_direct_rna/intergenic_mean_coverage/simple_gene_downstream_bins_batch_1.bed"

#Make sure trailing slash
output_dir = "/SAN/vyplab/sbs_projects/data/nanopore/IPSC_TDP43_KD_direct_rna/bedtools_strand_specific_intergenic_coverage/"

#########---------------
## BEDTOOLS PARAMETERS
#########---------------

# If don't want to use these options, assign them to "" (unless told otherwise)

# -s - only report hits overlapping on the same strand
# -S - only report hits overlapping on opposite strands
strandedness = "-s"

# Report the depth at each position in each interval
depth = "-d"

# what operations to summarise coverage in each intervals from bed_path? All below are valid strings
# sum, count, count_distinct, min, max, last
# median, mode, antimode, stdev, sstdev, collapse
# distinct, concat, freqasc, freqdesc, first, last
# If just want per-base coverages, assign to an empty list with []
operations = ["mean","median"]
#1-based - which column should operation be performed on? (with depth = -d, this is the 6th column)
operation_column = 6

# No need to change these two
sorted = "-sorted"
bedtools_path = "/SAN/vyplab/alb_projects/tools/bedtools"

########-----------------

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(bam_dir) if f.endswith(bam_suffix)]

print(SAMPLES)

if not os._exists(output_dir):
    os.system("mkdir -p {0}".format(output_dir))


########-----------------

rule all:
    input:
        expand(output_dir + "{sample}.coverage.per_base.tsv", sample = SAMPLES),
        expand(output_dir + "{sample}.coverage.{operation}.tsv" if len(operations) > 0 else [], sample = SAMPLES, operation = operations)

# Sort on chromosome and then by start position - allows use of -sorted option which is less memory intensive on large files
rule sort_bed:
    input:
        bed_path

    output:
        temp(output_dir + "sorted.bed")

    shell:
        """
        sort -k1,1 -k2,2n {input} > {output}
        """

rule bedtools_coverage:
    input:
        bam = os.path.join(bam_dir, "{sample}" + bam_suffix),
        bed = output_dir + "sorted.bed"

    output:
        output_dir + "{sample}.coverage.per_base.tsv"

    params:
        path = bedtools_path,
        strand = strandedness,
        per_base = depth,
        sorted = sorted

    shell:
        """
        {params.path} coverage \
        -a {input.bed} \
        -b {input.bam} \
        {params.strandedness} \
        {params.per_base} \
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

    shell:
        """
        {params.path} groupby -i {input} \
        -g 1,2,3 \
        -c {params.op_col} \
        -o {wildcards.operation} > {output}
        """
