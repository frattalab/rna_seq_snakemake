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
WRITEFOLDER=submissions/$FOLDER
mkdir -p $WRITEFOLDER


snakemake -s rules/tpmcalculator.smk \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -R y -l h_vmem=8G,h_rt=8G -o $WRITEFOLDER" \
-j 50 \
--nolock \
--rerun-incomplete \
--latency-wait 100
