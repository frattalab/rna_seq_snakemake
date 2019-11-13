
bam_dir = "/SAN/vyplab/alb_projects/data/ibm_geo/STAR_aligned/"
te_gtf = "/Users/annaleigh/Documents/data/italian_muscle/testing_te/GRCh38_rmsk_TE.gtf"
gtf =  "/Users/annaleigh/Documents/data/italian_muscle/testing_te/gencode.v31.no_chr.annotation.gtf"

SAMPLES, = glob_wildcards(bam_dir + "{sample}.Aligned.sorted.out.bam")

rule all:
  input:


rule run_te:
    input:
        sample_bam = expand("{bam_dir}{sample}.Aligned.sorted.out.bam", sample=SAMPLES)
    output:
    shell:
        """
        set +u;
        source activate tetranscripts
        TEcount -b {input.sample_bam} \
        --sortByPos --GTF \
        --TE  .gz \
        --project
