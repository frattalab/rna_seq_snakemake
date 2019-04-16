#!/bin/bash
#Submit to the cluster, give it a unique name

if [ "$0" != "" ]; then
    RUN_NAME=$0
else
    RUN_NAME=$""
fi

FOLDER=$(date +"%Y%m%d%H%M")
mkdir -p $FOLDER
cp config/config.yaml $FOLDER/$RUN_NAME.config.yaml

snakemake -s rna_seq.snakefile \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -o $FOLDER" -j 100