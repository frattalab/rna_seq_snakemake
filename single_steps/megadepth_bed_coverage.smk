#configfile: "config.yaml"
import os


'''
Quick and dirty rule to output BAM mean coverage over intervals from BED file with megadepth
Run megadepth at command line to see what else it can do and chop and change at this
'''

bam_dir = "/SAN/vyplab/alb_projects/data/ward_bams/bam_files/"

# Pick out bams from bam_dir with this bam_suffix
# Also used to strip from file to get sample name
# e.g. Cont-B_S2.pass2Aligned.sortedByCoord.out.bam becomes Cont-B_S2
bam_suffix = ".pass2Aligned.sortedByCoord.out.bam"

bed_path = "/SAN/vyplab/sbs_projects/data/nanopore/IPSC_TDP43_KD_direct_rna/intergenic_mean_coverage/simple_gene_downstream_bins_batch_1.bed"

#sum or mean - what operation to return on coverage in intervals from bed_path?
operation = "mean"

#Make sure trailing slash
output_dir = "/SAN/vyplab/data/nanopore/IPSC_TDP43_KD_direct_rna/intergenic_mean_coverage/"

megadepth_path = "/SAN/vyplab/alb_projects/tools/megadepth"


SAMPLES = [f.rstrip(bam_suffix) for f in os.listdir(bam_dir) if f.endswith(bam_suffix)]

print("SAMPLES")


rule all:
    input:
        expand(output_dir + "{sample}.annotation.tsv", sample = SAMPLES)


rule megadepth_coverage:
    input:
        bam = os.path.join(bam_dir, "{sample}" + bam_suffix),
        bed = bed_path

    output:
        output_dir + "{sample}.annotation.tsv"

    params:
        path = megadepth_path,
        prefix = output_dir + "{sample}",
        op = operation,
        threads = 4,
        no_std_out = "--no-annotation-stdout"

    shell:
        """
        {params.path} {input.bam} \
        --prefix {params.prefix} \
        {params.no_std_out} \
        --threads {params.threads} \
        --annotation {input.bed} \
        --op {params.op}
        """
