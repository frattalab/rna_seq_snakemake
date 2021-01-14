import os
# a top level folder where the bams reside
project_dir = "/SAN/vyplab/alb_projects/data/4su_full_ward_tdp_kd_ipsc/"
out_spot = "salmon/"
fastq_spot = "merged_fastqs/"

salmon_strand_info = "-l ISF"
end_type = "pe"
# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------

output_dir = os.path.join(project_dir,out_spot)
fastq_dir = os.path.join(project_dir,fastq_spot)

SAMPLES, = glob_wildcards(fastq_dir + "{sample}_1.merged.fastq.gz")
print(SAMPLES)
rule all:
  input:
    expand(output_dir + "{sample}_" + "quant.sf", sample = SAMPLES)
rule salmon_quant:
    input:
        fast1 = fastq_dir  + "{sample}_1.merged.fastq.gz",
        fast2 = fastq_dir  + "{sample}_2.merged.fastq.gz",
    output:
        output_dir + "{sample}_" + "quant.sf"
    params:
        tempout = output_dir + "{sample}",
        tempout2 = output_dir + "{sample}" + "/quant.sf"
    threads: 2
    shell:
        """
        /SAN/vyplab/alb_projects/tools/salmon-latest_linux_x86_64/bin/salmon quant\
         -i /SAN/vyplab/vyplab_reference_genomes/salmon/transcriptome_index \
         -l {salmon_strand_info} \
         -1 {input.fast1} \
         -2 {input.fast2} \
         -o {params.tempout}\
         --geneMap /SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v34.annotation.gtf
         mv {params.tempout2} {output}
        """
