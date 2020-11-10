The only thing you should need to edit is the config.yaml file. All directory paths should end with a trailing /

### How to make a sample sheet
The necessary columns are
`sample_name`	`unit`	`fast1`	`fast2`	`group`	`exclude_sample_downstream_analysis`

After that, you can include as many columns with sample metadata as you want. 


#### What to put in each columns
`sample_name` - the name that you want each sample to have, for clarities sake I recommend using snake_case. Please no spaces, and don't start with a number. 
`unit` - this is for the situation when there are multiple fastq's per sample, just as the following case


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
4. Generate gene (raw) read count tables with FeatureCounts
5. Generate gene TPM values with TPMcalculator


Submit to the cluster using the following command. The **first argument** should be your **workflow of choice** (*{workflow}*), followed by a specific **run name for the job** (*{run_name}*, optional but recommended).
`source submit.sh {workflow} {run_name}`


You should have Snakemake executable from your path. Check that this is the case before submission by typing in "which snakemake" at the command line. It should say this:
"/share/apps/python-3.7.2-shared/bin/snakemake""

To add Snakmake to your PATH. You'll also need to be using the correct version of Python. I recommend doing the following.

Open your .bash_profile:
"nano ~/.bash_profile"

Add this line to your .bash_profile to source the Cluster Folk's file for setting your Python version

"source /share/apps/examples/source_files/python/python-3.7.2.source"

Close the .bash_profile, source it, "source ~/.bash_profile"

Confirm that snakemake is callable

"which snakemake" should now say "/share/apps/python-3.7.2-shared/bin/snakemake"

In the samples csv sheet, the unit column is used for when the fastq's for a single sample have been split. E.G if one sample has multiple fastqs this will tell you that. If you only have one fastq per sample, the unit column must not be empty. Please fill it out with either the sample name or just some place holder text like "a". It will not work if you just put a number.

Before submission I recommend that you test that everything looks correct with a dry run first. This can be done with
"snakemake -s rna_seq.snakefile -n -p"
