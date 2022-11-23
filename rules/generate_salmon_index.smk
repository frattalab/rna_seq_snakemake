import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

INDEX_DIR = config["salmon_indices"]
SPECIES = config["species"]
SPECIES_VERSION = get_species_version(SPECIES)
ANNOTATION_VERSION = get_annotation_version(SPECIES)
DECOY_TYPE = config["salmon_index_type"]
KMER_SIZE = config["salmon_index_kmer_size"]

# dir for Index generated for each species's genome assembly, each type of decoy sequence (homologous only or whole genome), provided txome annotation and kmer size
TXOME_DIR = salmon_target_index(INDEX_DIR, SPECIES, SPECIES_VERSION, DECOY_TYPE, ANNOTATION_VERSION, KMER_SIZE)
print(TXOME_DIR)

if not os.path.exists(TXOME_DIR):
    os.system("mkdir -p {}".format(TXOME_DIR))


# 1 decoys (& merged txome + decoys FA ) file generated for each genome assembly + transcriptome annotation version
# Used as input to salmon index
DECOYS_DIR = os.path.join(INDEX_DIR, SPECIES, SPECIES_VERSION, "decoys", DECOY_TYPE, ANNOTATION_VERSION, "")
print(DECOYS_DIR)

if not os.path.exists(DECOYS_DIR):
    os.system("mkdir -p {}".format(DECOYS_DIR))


###### RULE ORDER DIRECTIVE
## To make deciding the input to index easier, I gave index the same input file path
## Since both partial and full output.produce the same file, I need to enforce rule order depending on DECOY_TYPE to prevent an AmbiguousRuleException

if DECOY_TYPE == "full":

    ruleorder: generate_full_decoys > generate_partial_decoys

elif DECOY_TYPE == "partial":

    ruleorder: generate_partial_decoys > generate_full_decoys

else:
    raise ValueError("{} is invalid value for salmon_index_type. Must be one of 'full' or 'partial'".format(DECOY_TYPE))


###### WORKFLOW DEFINITION

rule all_salmon_index:
    input:
        os.path.join(TXOME_DIR, "seq.bin"),
        os.path.join(TXOME_DIR, "pos.bin") # Two targets just to be extra


# https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode
# https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh
rule generate_partial_decoys:
    input:
        genome_fa = get_genome_fasta(SPECIES),
        txome_fa = get_transcriptome_fasta(SPECIES),
        gtf = get_gtf(SPECIES)

    output:
        os.path.join(DECOYS_DIR, "gentrome.fa"),
        os.path.join(DECOYS_DIR, "decoys.txt")

    params:
        outdir = DECOYS_DIR,
        script = config["salmon_decoy_shell_script"],
        bedtools = config["bedtools_path"],
        mashmap = config["mashmap_path"]

    threads:
        4

    shell:
        """
        bash {params.script} \
        -j {threads} \
        -a {input.gtf} \
        -g {input.genome_fa} \
        -t {input.txome_fa} \
        -b {params.bedtools} \
        -m {params.mashmap} \
        -o {params.outdir}
        """

# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
rule generate_full_decoys:
    input:
        genome_fa = get_genome_fasta(SPECIES),
        txome_fa = get_transcriptome_fasta(SPECIES),

    output:
        gentrome_fa = os.path.join(DECOYS_DIR, "gentrome.fa"),
        decoys = os.path.join(DECOYS_DIR, "decoys.txt")

    params:
        outdir = DECOYS_DIR

    shell:
        """
        grep "^>" {input.genome_fa} | cut -d " " -f 1 > {output.decoys}
        sed -i.bak -e 's/>//g' {output.decoys}

        cat {input.txome_fa} {input.genome_fa} > {output.gentrome_fa}
        """


rule salmon_index:
    input:
        gentrome_fa = os.path.join(DECOYS_DIR, "gentrome.fa"),
        decoys = os.path.join(DECOYS_DIR, "decoys.txt")

    output:
        os.path.join(TXOME_DIR, "seq.bin"),
        os.path.join(TXOME_DIR, "pos.bin")

    params:
        k = KMER_SIZE,
        outdir = TXOME_DIR,
        gencode = "--gencode" if config["transcriptome_source"] == "gencode" else ""

    threads:
        4

    conda:
        "../env/align.yaml"

    shell:
        """
        salmon index \
        -t {input.gentrome_fa} \
        -i {params.outdir} \
        --decoys {input.decoys} \
        -k {params.k} \
        {params.gencode} \
        -p {threads}
        """
