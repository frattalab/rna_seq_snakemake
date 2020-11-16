import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

feature_counts_outdir = get_output_dir(config["project_top_level"], config["feature_counts_output_folder"])
star_outdir = get_output_dir(config["project_top_level"], config['star_output_folder'])

#make sure the output folder for featureCounts exists before running anything
os.system("mkdir -p {0}".format(feature_counts_outdir))
#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species
REFERENCE_ANNOTATION = get_gtf(config['species'])

BASES, CONTRASTS = return_bases_and_contrasts('config/DESeq2comparisons.yaml')

rule deseqOutput:
    input:
        expand(os.path.join(config['majiq_top_level'],"delta_psi_voila_tsv","{bse}_{contrast}" + ".psi.tsv"),zip, bse = BASES,contrast = CONTRASTS)
