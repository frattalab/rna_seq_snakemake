configfile: "config.yaml"

rule fastqc:
	input:
		get.all.fastqfiles
	output:

rule star:
	input:
