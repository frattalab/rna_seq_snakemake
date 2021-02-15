import os
# a top level folder where the bams reside
project_dir = "/home/annbrown/data/ursa_mouse"
out_spot = "name_sortbams/"
bam_spot = "ori_bams/"
fastq_dir = "pulled_fastq/"
bam_suffix = ".bam"
end_type = "pe"
# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------

output_dir = os.path.join(project_dir,out_spot)
bam_dir = os.path.join(project_dir,bam_spot)
fastq_dir = os.path.join(project_dir, fastq_dir)

SAMPLES, = glob_wildcards(bam_dir + "{sample}" + bam_suffix)
print(SAMPLES)

rule all:
  input:
    expand(fastq_dir + "{sample}_1.merged.fastq.gz", sample = SAMPLES)

rule name_sort:
    input:
        aligned_bam = bam_dir + "{sample}" + bam_suffix
    output:
       temp(output_dir + "{sample}_namesorted.bam")
    shell:
        """
        mkdir -p {output_dir}
        samtools sort -n -@ 2 {input.aligned_bam} -o {output}
        """
if end_type == "pe":
  rule bam_to_fastq:
      input:
          name_sort_bam = temp(output_dir + "{sample}_namesorted.bam")
      output:
          one = temp(fastq_dir + "{sample}_1.merged.fastq"),
          two = temp(fastq_dir + "{sample}_2.merged.fastq")
      shell:
          """
          bedtools bamtofastq -i {input} \
                        -fq {output.one} \
                        -fq2 {output.two}
          """
  rule gunzip_fastq:
      input:
          one = temp(fastq_dir + "{sample}_1.merged.fastq"),
          two = temp(fastq_dir + "{sample}_2.merged.fastq")
      output:
          one_out = fastq_dir + "{sample}_1.merged.fastq.gz",
          two_out = fastq_dir + "{sample}_2.merged.fastq.gz"
      shell:
          """
          gzip {input.one}
          gzip {input.two}
          """
else:
  rule bam_to_fastq:
      input:
          name_sort_bam = temp(output_dir + "{sample}_namesorted.bam")
      output:
          one = temp(fastq_dir + "{sample}_1.merged.fastq")
      shell:
          """
          bedtools bamtofastq -i {input} \
                        -fq {output.one}
          """
  rule gunzip_fastq:
      input:
          one = fastq_dir + "{sample}_1.merged.fastq"
      output:
          one_out = fastq_dir + "{sample}_1.merged.fastq.gz"
      shell:
          """
          gzip {input.one}
          """
