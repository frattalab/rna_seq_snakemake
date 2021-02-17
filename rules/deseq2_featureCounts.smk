import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

#reading in the samples and dropping the samples to be excluded in order to get a list of sample names
samples = pd.read_csv(config['sampleCSVpath'])
samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
SAMPLE_NAMES = list(set(samples2['sample_name'] + config['bam_suffix']))
GROUPS = list(set(samples2['group']))

BASES, CONTRASTS = return_bases_and_contrasts('config/DESeq2comparisons.yaml')
print(BASES)
print(CONTRASTS)

DESEQ2_DIR = get_output_dir(config['project_top_level'], config['deseq2'])
DESEQ2_DIR = DESEQ2_DIR + "featureCounts/"


option_list = list(
    make_option(c("-f", "--folder_of_featurecounts"), type="character", default=NULL,
                help="a folder containing featurecounts output", metavar="character"),
    make_option(c("-b", "--base_grep"), type="character", default=NULL,
        help="a folder containing featurecounts output", metavar="character"),
    make_option(c("-c", "--contrast_grep"), type="character", default=NULL,
                help="output file name", metavar="character")
    make_option(c("-o", "--out"), type="character", default="out.txt",
          help="output file name; will output 3 files", metavar="character")
    make_option(c("-o", "--out"), type="character", default="out.txt",
          help="output file name; will output 3 files", metavar="character")
);

rule deseqOutput:
    input:
        expand(os.path.join(DESEQ2_DIR,"{bse}_{contrast}" + ".psi.tsv"),zip, bse = BASES,contrast = CONTRASTS)

rule run_standard_deseq:
    input:
        base_group = lambda wildcards: sample_names_from_contrast(wildcards.bse),
        contrast_group = lambda wildcards: sample_names_from_contrast(wildcards.contrast)
    output:
    params:
    shell:
    """
    Rscript standard_deseq2_command_line.R \
    --folder_of_featurecounts {input.csv} \
    --base_grep THIS NEEDS TO BE GENERATED AS A PASTED TOGETHER of the sample names \
    -- contrast_grep
    --out {output}
    """
