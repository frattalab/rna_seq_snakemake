
See links to instruction sections below:

- [Designing a sample sheet](#how-to-make-a-sample-sheet)
- [Running the pipeline](#how-to-run-the-pipeline)
- [Testing and updating the pipeline](#testing-and-updating-the-pipeline)
- [Checklist before pull requests](#what-to-do-before-pull-requests-or-feature-additions)

The only thing you should need to edit is the config.yaml file (and DESeq2comparisons.yaml file if you are running DESeq2). All directory paths should end with a trailing /

### How to make a sample sheet
The necessary columns are
`sample_name`	`unit`	`fast1`	`fast2`	`group`	`exclude_sample_downstream_analysis`

After that, you can include as many columns with sample metadata as you want.


#### What to put in each columns
`sample_name` - the name that you want each sample to have, for clarities sake I recommend using snake_case. Please no spaces, and don't start with a number, and don't use an periods in the sample name. 


`unit` - this is for the situation when there are multiple fastq's per sample, such as as the following case

sample_name | unit | fast1 | fast2 | group | exclude_sample_downstream_analysis
-- | -- | -- | -- | -- | --
sample_one | WTCHG_598911_202118 | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_598911_202118_1.fastq.gz | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_598911_202118_2.fastq.gz | OPMD |
sample_one | WTCHG_606179_202118 | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_606179_202118_1.fastq.gz | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_606179_202118_2.fastq.gz | OPMD |
sample_two | WTCHG_598911_203130 | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_598911_203130_1.fastq.gz | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_598911_203130_2.fastq.gz | IBM |
sample_two | WTCHG_606179_203130 | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_606179_203130_1.fastq.gz | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_606179_203130_2.fastq.gz | IBM |
sample_two | WTCHG_606180_203130 | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_606180_203130_1.fastq.gz | /SAN/vyplab/IoN_RNAseq/Bilal_Muscle_biopsies/fastqs/WTCHG_606180_203130_2.fastq.gz | IBM |


This is to allow trimming to occur on each fastq independently and then the merging is done later, if you don't have multiple fastqs per sample
please just put a placeholder there. If downloaded from SRA, I recommend the accession, but just don't leave it empty, you could fill  it with "placeholder".

sample_name | unit | fast1 | fast2 | group | exclude_sample_downstream_analysis
-- | -- | -- | -- | -- | --
axonal_control_1 | SRR11430624 | /SAN/vyplab/alb_projects/data/briese_tdp43_mouse_motorneuron/raw_data/SRR11430624.fastq |  | axonal_control |
axonal_control_2 | SRR11430625 | /SAN/vyplab/alb_projects/data/briese_tdp43_mouse_motorneuron/raw_data/SRR11430625.fastq |  | axonal_control |
axonal_control_3 | SRR11430626 | /SAN/vyplab/alb_projects/data/briese_tdp43_mouse_motorneuron/raw_data/SRR11430626.fastq |  | axonal_control |
axonal_control_4 | SRR11430627 | /SAN/vyplab/alb_projects/data/briese_tdp43_mouse_motorneuron/raw_data/SRR11430627.fastq |  | axonal_control |

`fast1`	`fast2` - paths to the fast1 and fast2 files, if data is single end leave fast2 blank

`group` - a group name, useful for downstream analysis

`exclude_sample_downstream_analysis` - if a sample ends up being dirty or you don't want to analyze it, place a 1, otherwise leave blank

(Should be covered above - feel free to skip this paragraph...)

In the samples csv sheet, the unit column is used for when the fastq's for a single sample have been split. E.G if one sample has multiple fastqs this will tell you that. If you only have one fastq per sample, the unit column must not be empty. Please fill it out with either the sample name or just some place holder text like "a". It will not work if you just put a number.

(End of redundant paragraph)

### How to run the pipeline

The pipeline has specific defined workflows. These are currently:

#### fastq_qc
1. Trim reads with fastp
2. Generate QC reports with FASTQC

#### interleave_fastq_qc
1. Trim reads with fastp
2. Generate QC reports with FASTQC
3. Combine paired end read files into single, interleaved FASTQ file with bbmap

#### align
1. Trim reads with fastp
2. Generate QC reports with FASTQC
3. Align reads to genome with STAR

#### salmon
1. Trim reads with fastp
2. Generate QC reports with FASTQC
3. Quantify transcripts with Salmon

**NB: salmon workflow currently only supports PAIRED-END READS**

#### DESeq2
1. Generate tx2gene mapping table for mapping genes to transcripts for DE from Salmon counts
2. Run differential expression test with DESeq2

#### Minimal dependencies

**You should have Snakemake executable from your PATH**. Check that this is the case before submission by typing in "which snakemake" at the command line. It should say something like this:
"/share/apps/python-3.7.2-shared/bin/snakemake""

To add Snakmake to your PATH. You'll also need to be using the correct version of Python. I recommend doing the following.

Open your .bash_profile:
"nano ~/.bash_profile"

Add this line to your .bash_profile to source the Cluster Folk's file for setting your Python version

"source /share/apps/examples/source_files/python/python-3.7.2.source"

Close the .bash_profile, source it, "source ~/.bash_profile"

Confirm that snakemake is callable

"which snakemake" should now say "/share/apps/python-3.7.2-shared/bin/snakemake"

#### Submitting to cluster

Before submission I recommend that you test that everything looks correct with a dry run first. This can be done with:

```
snakemake -n -p -s workflows/{workflow}.smk
```

Submit to the cluster using the following command. The **first argument** should be your **workflow of choice** (*{workflow}*), followed by a specific **run name for the job** (*{run_name}*, optional but recommended).

```
source submit.sh {workflow} {run_name}
```


### Testing and Updating the pipeline

If you make changes to the pipeline, **be sure to test the pipeline works correctly before committing/submitting a pull request**. I've made test files and datasets for *paired end* and *single end* runs, with corresponding downsampled FASTQs stored in our workspace on the cluster.

The configs and sample tables are ready to go, so to run the test datasets follow instructions below. If you perform a dry run and get a message like `nothing to be done`, go to the output test directory in vyplab_reference_genomes (you can find path in `config/test_{type}_config.yaml`) and delete the output. I recommend using `align` workflow, as this should include all possible rules in the pipeline.

#### Single end

**Dry run**

```
snakemake -n -p -s workflows/align.smk --configfile config/test_se_config.yaml
```

**Submit to cluster**

```
source submit_test_se.sh align {optional_run_name}
```

#### Paired end

**Dry run**

```
snakemake -s workflows/align.smk --configfile test_data_configs/test_pe_config.yaml -c 4 --use-conda
```

**Run locally on an interactive node**

```
qrsh -l tmem=16G,h_vmem=16G

snakemake -s workflows/align.smk --configfile test_data_configs/test_pe_config.yaml -c 4 --use-conda

```


### What to do before pull requests or feature additions

Please make sure you have done the following before committing changes/submitting a pull request:
- Ran pipeline (with test datasets) without workflow errors
- Updated all **config/*_config.yaml** (i.e. including test examples) with new parameters
- Updated README or config with additional instructions/comments (if adding functionality)

If you're adding a **new workflow**, also make sure to include a **submit script**.

(this is mostly a reminder for me...(*Sam*))


 multiqc         -p         -o test_data_analyzed/paired_end/multiqc/align/ test_data_analyzed/paired_end/qc/fastqc/ test_data_analyzed/paired_end/fastp_trimmed/ test_data_analyzed/paired_end/STAR_aligned/ test_data_analyzed/paired_end/feature_counts/
