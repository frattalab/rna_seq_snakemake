#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash
#$ -h_rt 72:00:00
#$ -cwd
#$ -V
#$ -l h_vmem=2G,h_rt=20:00:00,tmem=2G


# join stdout and stderr output
#$ -j y
#$ -R y



snakemake -s  snakemake -s upload_to_ena.smk \
--nolock \
--rerun-incomplete \
--latency-wait 100
