import pandas as pd
import os
import subprocess

import subprocess

from datetime import datetime

configfile: "config/config.yaml"

rule all:
	input:
		config["fastqc_output_folder"] + "multiqc_report.html"

include: "rules/fastqc.smk"
include: "rules/multiqc.smk"