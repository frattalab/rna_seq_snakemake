#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G

# join stdout and stderr output
#$ -j y
#$ -R y


if [ "$1" != "" ]; then
    RUN_NAME=$1
else
    RUN_NAME=$""
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p $FOLDER
cp single_steps/te_count.smk $FOLDER/$RUN_NAME_te_count.smk

snakemake -s single_steps/te_count.smk \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -R y -l h_vmem=14G,h_rt=14G -o $FOLDER" \
-j 500 \
--nolock \
--rerun-incomplete \
--latency-wait 100
