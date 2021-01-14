rule salmon_index:
    input:
        "/SAN/vyplab/vyplab_reference_genomes/sequence/human/gencode/gencode.v34.transcripts.fa"
    output:
        directory("/SAN/vyplab/vyplab_reference_genomes/salmon/transcriptome_index")
    threads: 2
    shell:
        """
        /SAN/vyplab/alb_projects/tools/salmon-latest_linux_x86_64/bin/salmon index -t {input} -i {output} --threads {threads} --gencode
        """

