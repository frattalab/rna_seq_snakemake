import os
project_dir = "/SAN/vyplab/alb_projects/data/muscle/analysis/"

bam_dir = os.path.join(project_dir,"STAR_aligned/")

# TODO ACTUALLY HAVE THIS READ FROM THE CONFIG
te_gtf = "/SAN/vyplab/vyplab_reference_genomes/annotation/human/transposable_elements/GRCh38_rmsk_TE.gtf"
gtf =  "/SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v31.no_chr.annotation.gtf"

# te_gtf = "/SAN/vyplab/vyplab_reference_genomes/annotation/mouse/transposable_elements/mm10_rmsk_TE.gtf"
# gtf =  "/SAN/vyplab/vyplab_reference_genomes/annotation/mouse/gencode/gencode.vM22.annotation.gtf"

output_dir = os.path.join(project_dir,"te_count_strand/")

SAMPLES, = glob_wildcards(bam_dir + "{sample}.Aligned.sorted.out.bam")

print(SAMPLES)
strand = "--stranded reverse"
rule all:
  input:
    expand(output_dir + "{sample}.cntTable", sample = SAMPLES)


rule run_te:
    input:
        sample_bam = lambda wildcards: bam_dir + "{sample}.Aligned.sorted.out.bam"
    output:
        output_dir + "{sample}.cntTable"
    params:
        strandness = strand
    shell:
        """
        set +u;
        source activate tetranscripts
        TEcount -b {input.sample_bam} \
        --sortByPos --GTF  {gtf}\
        --TE {te_gtf}\
        --project {output_dir}{wildcards.sample} \
        {params.strandness}
        """
