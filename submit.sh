#Submit to the cluster using the following command

FOLDER=$(date +"%Y%m%d%H%M")
mkdir -p $FOLDER
cp config/cluster.yaml $FOLDER/cluster.yaml

snakemake -s rna_seq.snakefile \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -o $FOLDER" -j 100