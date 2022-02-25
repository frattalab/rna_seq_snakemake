'''
rules to run
'''

rnaseqc_outdir = get_output_dir(config["project_top_level"], config[""])
log_outdir = get_output_dir(config["project_top_level"],
                            config["logs_output_folder"])

rule rnaseqc:
    '''
    '''
    input:
        bam = rules.sort_bam.output,
        bai = rules.index_sorted_bam.output,
        collapsed_gtf = rules.gtf_to_collapsed_gtf.output.collapsed_gtf

    output:
        metrics = os.path.join(rnaseqc_outdir, "{sample}.metrics.tsv"),
        coverage = os.path.join(rnaseqc_outdir, "{sample}.coverage.tsv")

    params:
        sample_id = "{sample}",
        outdir = get_output_dir(config["project_top_level"], config["rnaseqc_output_folder"]),
        verbose = "-v",
        unpaired = "--unpaired" if config["end_type"] != "pe" else "",
        coverage = "--coverage" # output coverage stats for each transcript

    log:
        os.path.join(log_outdir,
                     "rnaseqc",
                     "{sample}.rnaseqc.log")

    conda:
        ""

    shell:
        """
        rnaseqc {input.collapsed_gtf} \
        {input.bam} \
        {params.outdir} \
        -s {params.sample_id} \
        {params.verbose} \
        {params.coverage} \
        {params.unpaired} \
        &> {log}
        """
