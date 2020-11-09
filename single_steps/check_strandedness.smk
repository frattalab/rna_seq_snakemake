

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


WORKFLOW="workflows/${1}.smk"

if [ "$2" != "" ]; then
    RUN_NAME="$1"_"$2"
else
    RUN_NAME=$1
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")
#read in a samples table
SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
#so I want this rule to be run ONCE for every fast1, so the wild cards I'm giving are the 'name' of the fastq of the first read
#name = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in SAMPLES['fast1'].tolist()]
name1 = SAMPLES['fast1'][0]
if config['end_type'] == 'pe':
    name1 = SAMPLES['fast2'][0]

check_strandedness  --gtf /SAN/vyplab/vyplab_reference_genomes/annotation/human/ensembl/Homo_sapiens.GRCh38.100.gtf --transcripts /SAN/vyplab/vyplab_reference_genomes/sequence/human/ensembl/Homo_sapiens.GRCh38.cdna.all.fa.gz --reads_1 /SAN/vyplab/alb_projects/data/liu_facs_neurons/raw_data/SRR8571951_1.fastq.gz --reads_2 /SAN/vyplab/alb_projects/data/liu_facs_neurons/raw_data/SRR8571951_2.fastq.gz -n 10000


rule strandedOutput:

rule check_strandedness:
    input:
        sample_bam = lambda wildcards: bam_dir + "{sample}.Aligned.sorted.out.bam"
    output:
        output_dir + "{sample}.cntTable"
    params:
        strandness = strand
    shell:
        """
        mkdir -p {output_dir}
        set +u;
        source activate tetranscripts
        TEcount -b {input.sample_bam} \
        --sortByPos --GTF  {gtf}\
        --TE {te_gtf}\
        --project {output_dir}{wildcards.sample} \
        {params.strandness}
        """
