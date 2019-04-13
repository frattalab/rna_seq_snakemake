import pandas as pd
import os
configfile: "config/config.yaml"

#because the config file is so variable, I want to write out a copy of the config file along with system time everytime we run
def write_config_copy():
	os.system("mkdir -p snakemake_configs")
	os.system("")

rule all:
	input:
		config["fastqc_output_folder"] + "multiqc.html"

include: "rules/fastqc.smk"
include: "rules/multiqc.smk"