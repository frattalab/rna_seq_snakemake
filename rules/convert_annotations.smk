import re

'''
Steps to convert GTFs to alternate annotation formats required by certain tools
Each of these steps will output the alternate file to the same location (w/ a new file extension)
'''


rule gtf_to_refFlat:
    '''
    Rule to generate refFlat from GTF file (required by Picard...)
    Lifted from here: https://github.com/broadinstitute/picard/issues/805#issue-224851540
    '''
    input:
        gtf = get_gtf(config['species'])

    output:
        refflat = re.sub(".gtf$", ".refFlat.txt", get_gtf(config['species']))

    params:
        script = os.path.join(options["ucsc_tools_dir"], "gtfToGenePred"),
    shell:
        '''
        {params.script} \
        -genePredExt \
        -geneNameAsName2 \
        -ignoreGroupsWithoutExons \
        {input.gtf} \
        /dev/stdout | \
        awk 'BEGIN {{ OFS="\t"}} {{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}' \
        > {output.refflat}
        '''


rule gtf_to_bed12:
    '''
    Rule to generate BED12 file from GTF
    '''
    input:
        get_gtf(config['species'])

    output:
        bed12 = re.sub(".gtf$", ".bed12", get_gtf(config['species']))

    params:
        gtf2pred = os.path.join(options["ucsc_tools_dir"], "gtfToGenePred"),
        pred2bed = os.path.join(options["ucsc_tools_dir"], "genePredToBed"),

    shell:
        """
        {params.gtf2pred} {input} {output.bed12}.tmp
        {params.pred2bed} {output.bed12}.tmp {output}
        rm {output.bed12}.tmp
        """


rule gtf_to_collapsed_gtf:
    '''
    RNA-SeQC requires a 'collapsed' GTF - no overlapping tx on same strand, each gene has 1 transcript
    Generated using GTEx collapse annotation script
    '''
    input:
        get_gtf(config['species'])

    output:
        collapsed_gtf = re.sub(".gtf$", ".gtex_collapsed.gtf", get_gtf(config['species']))

    params:
        stranded = "--stranded", # RNA-SeQC says only wants non-overlapping on same strand
        #

    conda:
        # need numpy, bx-python & pandas
        ""

    shell:
        """
        python collapse_annotation.py \
        {input} \
        {output} \
        {params.stranded}
        """
