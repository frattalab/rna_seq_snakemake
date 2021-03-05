#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_vmem=8G,h_rt=24:00:00,tmem=8G


# join stdout and stderr output
#$ -j y
#$ -R y



snakemake -s  snakemake -s upload_to_ena.smk \
--nolock \
--rerun-incomplete \
--latency-wait 100
