import os
# a top level folder where the bams reside
project_dir = "/SAN/vyplab/alb_projects/data/muscle/analysis/"
out_spot = "feature_counts_strand_other/"
bam_spot = "STAR_aligned/"

# mouse and human gtf, comment dependening on your species
# gtf =  "/SAN/vyplab/vyplab_reference_genomes/annotation/mouse/gencode/gencode.vM22.annotation.gtf"
# gtf =  "/SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v31.annotation.gtf"
gtf =  "/SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v31.no_chr.annotation.gtf"

feature_counts_strand_info = "-s 1"
end_type = "pe"
# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------
feature_counts_path = "/SAN/vyplab/alb_projects/tools/subread-1.6.4-Linux-x86_64/bin/featureCounts"

output_dir = os.path.join(project_dir,out_spot)
bam_dir = os.path.join(project_dir,bam_spot)

SAMPLES, = glob_wildcards(bam_dir + "{sample}.bam")
print(SAMPLES)

rule all:
  input:
    expand(output_dir + "{sample}_featureCounts_results.txt", sample = SAMPLES)

rule feature_counts:
    input:
        aligned_bam = bam_dir + "{sample}.bam"
    output:
        out_name = output_dir + "{sample}_featureCounts_results.txt"
    params:
        ref_anno = gtf,
        stranded = feature_counts_strand_info
    run:
        shell("mkdir -p {output_dir}")
        if end_type == "pe":
            shell("{feature_counts_path} -p -t exon -g gene_id -a {params.ref_anno} -o {output.out_name} {params.stranded} {input.aligned_bam}")
        if end_type == "se":
            shell("{feature_counts_path} -a {params.ref_anno} -o {output.out_name} {params.stranded} {input.aligned_bam}")