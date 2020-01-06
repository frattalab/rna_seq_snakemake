#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y

if [ "$1" != "" ]; then
    RUN_NAME=$1
else
    RUN_NAME=$""
fi

FOLDER=$(date +"%Y%m%d%H%M")
WRITEFOLDER=submission/$FOLDER
mkdir -p $WRITEFOLDER
cp config/config.yaml $WRITEFOLDER/$RUN_NAME.config.yaml

snakemake -s single_steps/feature_counts.smk \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -R y -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -o $WRITEFOLDER" \
-j 500 \
--nolock \
--rerun-incomplete \
--latency-wait 100
