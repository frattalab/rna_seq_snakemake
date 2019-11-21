import os
project_dir = "/SAN/vyplab/alb_projects/data/4su_tdp_f210i/"

bam_dir = os.path.join(project_dir,"STAR_aligned_redone/")
te_gtf = "/SAN/vyplab/vyplab_reference_genomes/annotation/mouse/transposable_elements/mm10_rmsk_TE.gtf"
gtf =  "/SAN/vyplab/vyplab_reference_genomes/annotation/mouse/gencode/gencode.vM22.annotation.gtf"


output_dir = os.path.join(project_dir,"te_count/")

SAMPLES, = glob_wildcards(bam_dir + "{sample}.Aligned.sorted.out.bam")

rule all:
  input:
    expand(output_dir + "{sample}.cntTable", sample = SAMPLES)


rule run_te:
    input:
        sample_bam = lambda wildcards: bam_dir + "{sample}.Aligned.sorted.out.bam"
    output:
        output_dir + "{sample}.cntTable"
    shell:
        """
        set +u;
        source activate tetranscripts
        TEcount -b {input.sample_bam} \
        --sortByPos --GTF  {gtf}\
        --TE {te_gtf}\
        --project {output_dir}{wildcards.sample}
        """
