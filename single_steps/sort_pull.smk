import os

configfile: "config/single_steps_config.yaml"
# a top level folder where the bams reside

# Get a dict of options specific to sort_pull
options = config["sort_pull"]

nsort_out_spot = options["name_sort_outspot"]
bam_dir = options["bam_spot"]
fastq_dir = options["out_spot"]
bam_suffix = options["bam_suffix"]
end_type = options["end_type"]


SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(bam_dir) if f.endswith(bam_suffix)]

print(SAMPLES)

rule all:
    input:
        expand(os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(fastq_dir, "{sample}_2.merged.fastq.gz") if end_type == "pe" else [], sample = SAMPLES)

rule name_sort:
    input:
        aligned_bam = os.path.join(bam_dir, "{sample}" + bam_suffix)

    output:
        temp(os.path.join(nsort_out_spot, "{sample}_namesorted.bam"))

    params:
        samtools = options["samtools_path"]

    shell:
        """
        mkdir -p {output_dir}
        {params.samtools} sort -n -@ 2 {input.aligned_bam} -o {output}
        """


if end_type == "pe":
    rule bam_to_fastq:
        input:
            name_sort_bam = os.path.join(nsort_out_spot, "{sample}_namesorted.bam")
        output:
            one = temp(os.path.join(fastq_dir, "{sample}_1.merged.fastq")),
            two = temp(os.path.join(fastq_dir, "{sample}_2.merged.fastq"))

        params:
            bedtools = options["bedtools_path"]

        shell:
            """
            {params.bedtools} bamtofastq \
            -i {input} \
            -fq {output.one} \
            -fq2 {output.two}
            """

    rule gunzip_fastq:
        input:
            one = os.path.join(fastq_dir, "{sample}_1.merged.fastq"),
            two = os.path.join(fastq_dir, "{sample}_2.merged.fastq")
        output:
            one_out = os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz"),
            two_out = os.path.join(fastq_dir, "{sample}_2.merged.fastq.gz")
        shell:
            """
            gzip {input.one}
            gzip {input.two}
            """
else:
    rule bam_to_fastq:
        input:
            name_sort_bam = os.path.join(nsort_out_spot, "{sample}_namesorted.bam")

        output:
            one = temp(os.path.join(fastq_dir, "{sample}_1.merged.fastq"))

        params:
            bedtools = options["bedtools_path"]

        shell:
            """
            {params.bedtools} bamtofastq -i {input} \
            -fq {output.one}
            """

    rule gunzip_fastq:
        input:
            one = os.path.join(fastq_dir, "{sample}_1.merged.fastq")

        output:
            one_out = os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz")

        shell:
            """
            gzip {input.one}
            """
