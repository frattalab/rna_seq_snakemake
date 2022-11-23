import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

SPECIES_VERSION = get_species_version(config['species'])

GENOME_DIR = os.path.join(config['STAR_indices'],config['species'],SPECIES_VERSION,"star_indices_overhang" + str(config['readLen']))
print(GENOME_DIR)

rule all_star:
	input:
		GENOME_DIR + "/SA"

rule generate_genome:
	input:
		fasta = get_genome_fasta(config['species']),
		gtf = get_gtf(config['species'])
	output:
		GENOME_DIR + "/SA",
		GENOME_DIR + "/Genome"
	params:
		sjdbOverhang = config['readLen'] - 1
	threads:
		4
	conda:
		"../env/align.yaml"
	shell:
		"""
		STAR \
	    --runThreadN {threads} \
	    --runMode genomeGenerate \
	    --genomeDir {GENOME_DIR} \
	    --genomeFastaFiles {input.fasta} \
	    --sjdbGTFfile {input.gtf} \
	    --sjdbOverhang {params.sjdbOverhang} \
		--genomeSAsparseD 10
		"""
