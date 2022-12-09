import re
import pandas as pd
import os
import numpy as np
import subprocess
from subprocess import PIPE
import yaml

def get_fastq_names(DATA):
    samples = pd.read_csv(DATA, sep = ",")
    #take all the fasqs and combine them to a list remove any nas
    fastq_list = samples.fast1.dropna().tolist() + samples.fast2.dropna().tolist()
    # print("get_fastq_names - fastq_list values - {}".format(", ".join(fastq_list)))
    #if there are any missing values pandas gets annoyed so replace nans with empty string
    samples = samples.replace(np.nan, '', regex=True)

    #get the associated sample and unit name for each fastq using and or operator to get either fast1 or fast2, works for both single-end and paired end this way
    unit_name = [samples.loc[(samples['fast1'] == fq)| (samples['fast2'] == fq)].unit.iloc[0] for fq in fastq_list]

    #strip it down to just the name of the file bit
    #stripped = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in fastq_list]
    sample_names = samples.sample_name.tolist()
    #now combine these three lists to get a meaningfull and unique string for each fastq file

    return(sample_names, fastq_list, unit_name)

def return_fastq_location(wildcards):
	#return the file location from the fastq name
	return(ORDER_DICT[wildcards])

def return_parsed_extra_params(extra_params):
    #starting blank
    cmd = ""
    #for key in extra parameters
    for key in extra_params:
        #append the key value pair if it's a parmeter that needed something
        if extra_params[key]:
            cmd += " --{0} {1}".format(key,extra_params[key])
        else: #otherwise if it's parameter that's just a flag, append just the flag
            cmd += " --{0}".format(key)
    return(cmd)

def return_fastq(sample_name,unit, first_pair = True):
    SAMPLES = pd.read_csv(config["sampleCSVpath"])
    SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
    SAMPLES['fast1_name'] = [strpd.rpartition('/')[2].split(".")[0] for strpd in SAMPLES['fast1'].tolist()]

    if first_pair:
        return(SAMPLES.loc[(SAMPLES['sample_name'] == sample_name) & (SAMPLES['unit'] == unit)]["fast1"].values[0])
    else:
        return(SAMPLES.loc[(SAMPLES['sample_name'] == sample_name) & (SAMPLES['unit'] == unit)]["fast2"].values[0])

def get_trimmed(name):
    #the trimmed file is the output, and the unit, we find it from the sample and and the unit which snakemake wildcards are going through
    SAMPLES = pd.read_csv(config["sampleCSVpath"])
    SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

    fastp_outdir = get_output_dir(config['project_top_level'], config["fastp_trimmed_output_folder"])



    unit_fastqs = SAMPLES.loc[(SAMPLES.sample_name == name),'unit'].tolist()
    unit_fastqs2 = [u + "_" + name for u in unit_fastqs]

    trimmed_1 = [os.path.join(fastp_outdir,u + "_1.trimmed.fastq.gz") for u in unit_fastqs2]


    #if we have paired end data there will also be a trimmed 2, same thing, using the fast2 column instead
    if config['end_type'] == "pe":
        trimmed_2 = [os.path.join(fastp_outdir,u + "_2.trimmed.fastq.gz") for u in unit_fastqs2]

        #trimmed files is a list of the two
        trimmed_files = [trimmed_1, trimmed_2]

    else:
        trimmed_files = [trimmed_1]

    return(trimmed_files)



def return_all_trimmed(SAMPLES, pair = 1):
    #make a list of all the files trimme
    if pair == 1:
        list_of_trimmed = [config['fastp_trimmed_output_folder'] + \
        SAMPLES.loc[(SAMPLES.fast1 == fq),'unit'].tolist()[0] + \
        "/" +
        re.sub(".fastq.gz","",fq.rpartition('/')[2]) +\
        "_trimmed.fastq.gz" for fq in SAMPLES.fast1]

        return(list_of_trimmed)
    else:
        if config["end_type"] == "pe":

            list_of_trimmed = [config['fastp_trimmed_output_folder'] + \
            SAMPLES.loc[(SAMPLES.fast2 == fq),'unit'].tolist()[0] + \
            "/" +
            re.sub(".fastq.gz","",fq.rpartition('/')[2]) +\
            "_trimmed.fastq.gz" for fq in SAMPLES.fast2]

            return(list_of_trimmed)
        else:
            return(" ")

def get_processed_fastq(sample_name, pair=1, pair_prefix="_", trimmed_suffix=".trimmed.fastq.gz", merged_suffix=".merged.fastq.gz"):
    '''
    Return path to target processed FASTQ file for given sample and mate pair
    If sample only had a single FASTQ file in sample sheet, there is no need to do unnecessary 'merge_trimmed'
    --> return path to trimmed FASTQ
    If sample has multiple units (reads split across files), need to point to merged file
    '''

    # _1 or _2 after sample_name in filename
    rpair_str = pair_prefix + str(pair)

    fastp_outdir = get_output_dir(config['project_top_level'], config["fastp_trimmed_output_folder"])
    merged_outdir = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

    SAMPLES = pd.read_csv(config["sampleCSVpath"])
    SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
    # SAMPLES = SAMPLES.set_index("sample_name")

    if SAMPLES.loc[SAMPLES["sample_name"] == sample_name, "unit"].nunique() > 1:
        # sample had multiple fastq files, so need to point to merged FASTQ
        target_fastq = os.path.join(merged_outdir, sample_name + rpair_str + merged_suffix)

    else:
        # Sample had a single fastq file/unit, so need to point to trimmed FASTQ

        # Trimmed fastq names have {unit}_{sample_name}_{pair}{trimmed_suffix}
        unit = SAMPLES.set_index("sample_name").loc[sample_name, "unit"]
        target_fastq = os.path.join(fastp_outdir, unit + "_" + sample_name + rpair_str + trimmed_suffix)


    return target_fastq


def get_species_version(species):
    temp = pd.read_csv("config/reference_files_species.csv",sep = ",")
    return(temp.species_version[temp.species == species].tolist()[0])


def get_annotation_version(species):
    temp = pd.read_csv("config/reference_files_species.csv",sep = ",")
    return (temp.annotation_version[temp.species == species].tolist()[0])

def get_genome_fasta(species):
    temp = pd.read_csv("config/reference_files_species.csv",sep = ",")
    return(temp.genome_fa[temp.species == species].tolist()[0])

def get_transcriptome_fasta(species):
    temp = pd.read_csv("config/reference_files_species.csv",sep = ",")
    return (temp.transcriptome_fa[temp.species == species].tolist()[0])

def get_gtf(species):
    temp = pd.read_csv("config/reference_files_species.csv",sep = ",")
    return(temp.gtf[temp.species == species].tolist()[0])

def get_bed12(species):
    temp = pd.read_csv("config/reference_files_species.csv",sep = ",")
    return (temp.bed12[temp.species == species].tolist()[0])

# takes the featurcounts strand and returns the interpretation for kallisto_output_folder
def get_kallisto_strand(fcStrand):
    if fcStrand == "-s 0":
        return("")
    elif fcStrand == "-s 1":
        return("--fr-stranded")
    elif fcStrand == "-s 2":
        return("--rf-stranded")


def get_scallop_strand(fcStrand):
    '''
    Return string for scallop libtype denoting strandedness/orientation of library.
    Parses featureCounts strand and returns corresponding Scallop libtype string
    '''
    if fcStrand == "-s 0":
        return "unstranded"

    elif fcStrand == "-s 1":
        return "second"

    elif fcStrand == "-s 2":
        return "first"

def get_salmon_strand(fcStrand):
    '''
    Return string for salmon libtype denoting strandedness/orientation of library.
    Parses featureCounts strand and returns corresponding salmon libtype string
    (should include option to return the let Salmon infer for you string)
    '''
    if fcStrand == "-s 0":
        return "IU"

    elif fcStrand == "-s 1":
        return "ISF"

    elif fcStrand == "-s 2":
        return "ISR"

def get_collectRnaSeq_strand(fcStrand):
    if fcStrand == "-s 0":
        return("")
    elif fcStrand == "-s 1":
        return("FIRST_READ_TRANSCRIPTION_STRAND")
    elif fcStrand == "-s 2":
        return("SECOND_READ_TRANSCRIPTION_STRAND")


def get_output_dir(project_top_level, rule_output_path):
    '''
    Return path to output directory (with trailing slash)
    if rule_output_path is relative, then it is appended to project_top_level
    if rule_output_path is absolute then rule_output_path is returned
    '''

    if os.path.isabs(rule_output_path):
        if rule_output_path.endswith('/'):
            return rule_output_path
        else:
            return rule_output_path + "/"

    else:
        return os.path.join(project_top_level, rule_output_path, '')

def return_bases_and_contrasts(comparison_yaml):
    """
    returns all the bases and contrasts from the comparisons.yaml
    """

    with open(comparison_yaml, 'r') as stream:
        try:
            compare_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    base_keys = []
    contrast_keys = []

    for key in compare_dict:
        temp = compare_dict[key]
        for ind, k2 in enumerate(temp):
            if ind == 1:
                base_keys.append(k2)
            if ind == 2:
                contrast_keys.append(k2)
    return(base_keys,contrast_keys)


def salmon_target_index(txome_dir, species, species_version, decoy_type, annot_version, kmer_size):
    '''
    Which transcriptome index should Salmon use/generate?

    params
        txome_dir str
        path to master directory where salmon indexes are stored.

        decoy_type str ["full", "partial"]
        Which decoys sequences should be used alongside target transcriptome?
        full - whole genome treated as decoy
        partial - sequences in genome with high similarity to transcriptome sequences (>= 80 % sequence homology)

        species str
        What species is sequencing library generated from. most likely "human" or "mouse"

        kmer_size int
        Size of k-mers used to generate index. 31 recommended for 75+bp datasets.

    returns str
        path to directory (W/O TRAILING SLASH) of target decoy-aware transcriptome. concatenation of txome_dir, species, decoy_type and kmer_size
    '''

    if decoy_type not in ["full", "partial"]:

        raise ValueError("{0} is invalid value for decoy_type. Must be one of 'full' or 'partial'".format(decoy_type))

    else:
        return os.path.join(txome_dir, species, species_version, decoy_type, ".".join([annot_version, "kmer_" + str(kmer_size)]))



def multiqc_target_files(workflow_str, sample_names, fastq_prefixes, units):
    '''
    Returns list of target files for multiqc depending on the workflow being ran
    This is mostly manually defined (like our workflows) for now until I find a better solution
    Use this as an 'input function'
    Indexes of sample names, prefixes & units should match
    '''

    # Last one is a dummy value
    valid_workflows = ["fastq_qc", "align", "salmon", "multiqc_output"]

    if workflow_str not in valid_workflows:
        raise ValueError("{0} is not a supported value for 'workflow' - use one of {1}".format(workflow, ",".join(valid_workflows)))

    else:
        out_targets = []

        # Define all possible output directories for different qc metrics
        fastqc_outdir = get_output_dir(config["project_top_level"], config["fastqc_output_folder"])
        fastp_outdir = get_output_dir(config["project_top_level"], config["fastp_trimmed_output_folder"])
        star_outdir = get_output_dir(config["project_top_level"], config["star_output_folder"])
        salmon_outdir = get_output_dir(config["project_top_level"], config["salmon_output_folder"])
        rseqc_outdir = get_output_dir(config["project_top_level"], config["rseqc_output_folder"])
        # Define target files for each step
        targets_fastqc = expand(fastqc_outdir + "{sample}/{unit}/{fastq_prefix}_fastqc.html",zip, sample=sample_names, unit=units, fastq_prefix = fastq_prefixes)
        # print("fastqc targets - {0}".format(", ".join(targets_fastqc)))
        # print(targets_fastqc)
        # print(type(targets_fastqc))

        targets_fastp = expand(fastp_outdir + "{unit}_{name}_fastp.json", zip, name = sample_names, unit = units)

        # Created in same dir as STAR logs (but after bams generated)
        targets_star = expand(star_outdir + "{name}.flagstat.txt", name = sample_names)
        targets_salmon = expand(salmon_outdir + "{sample}/" + "quant.sf", sample = sample_names)

        rseq_target_suffixes = [".geneBodyCoverage.txt", ".infer_experiment.txt", ".inner_distance_freq.txt", ".junctionSaturation_plot.r", ".read_distribution.txt"]
        targets_rseqc = expand(rseqc_outdir + "{name}" + "{suffix}", name = sample_names, suffix = rseq_target_suffixes)
        # targets_rseqc.append(os.path.join(rseqc_outdir, "all_bams_output.geneBodyCoverage.txt"))


        if workflow_str == "fastq_qc":
            # Only need output from fastqc & fastp
            out_targets.extend(targets_fastqc)
            out_targets.extend(targets_fastp)

        elif workflow_str == "align":
            # Basically everything but salmon
            out_targets.extend(targets_fastqc)
            # print(out_targets)
            out_targets.extend(targets_fastp)
            out_targets.extend(targets_star)
            out_targets.extend(targets_rseqc)

            # print("out_targets for align - {0}".format(", ".join(out_targets)))

        elif workflow_str == "salmon":
            # Just fastq QC & salmon
            out_targets.extend(targets_fastqc)
            out_targets.extend(targets_fastp)
            out_targets.extend(targets_salmon)

        elif workflow_str == "multiqc_output":
            # This is dummy for if run independently (i.e. has no dependent rules so no targets)
            pass

    # print("this is out_targets for multiqc_target_files - {0}".format(",".join(out_targets)))

    return out_targets


def multiqc_target_dirs():
    '''
    Returns list of paths to directories for multiqc to scan for log files
    Since it scan recursively through dirs, and only penalty to searching extra dirs is added run-time
    For simplicity, this returns paths to all potential directories of different workflows, provided they exist / have been created prior to the DAG being generated
    '''

    outdir_keys = ["fastqc_output_folder", "fastp_trimmed_output_folder", "star_output_folder", "salmon_output_folder", "rseqc_output_folder"]

    # List of all output directories specified in config
    all_dir_paths = [get_output_dir(config["project_top_level"], config[x]) for x in outdir_keys]

    # print("this is all_dir_paths for multiqc {0}".format(",".join(all_dir_paths)))
    # Return only potential directories that have already exist
    # The directory must be made in the snakefile prior to the rule being run for this to work (otherwise the path doesn't exist when the DAG is generated)
    target_dir_paths = [p for p in all_dir_paths if os.path.exists(p)]

    # print("this is target_dir_paths for multiqc - {0}".format(",".join(target_dir_paths)))

    return target_dir_paths


def sample_names_from_contrast(grp):
    """
    given a contrast name or list of groups return a list of the files in that group
    """
    #reading in the samples
    samples = pd.read_csv(config['sampleCSVpath'])
    #there should be a column which allows you to exclude samples
    samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
    #read in the comparisons and make a dictionary of comparisons, comparisons needs to be in the config file
    compare_dict = load_comparisons()
    #go through the values of the dictionary and break when we find the right groups in that contrast
    grps, comparison_column = return_sample_names_group(grp)
    #take the sample names corresponding to those groups
    if comparison_column == "":
        print(grp)
        return([""])
    grp_samples = "|".join(set(list(samples2[samples2[comparison_column].isin(grps)].sample_name)))
    print(grp_samples)
    return(grp_samples)


def salmon_files_from_contrast(grp):
    """
    given a contrast name or list of groups return a list of the files in that group
    """
    #reading in the samples
    samples = pd.read_csv(config['sampleCSVpath'])
    #there should be a column which allows you to exclude samples
    samples2 = samples.loc[samples.exclude_sample_downstream_analysis != 1]
    #read in the comparisons and make a dictionary of comparisons, comparisons needs to be in the config file
    compare_dict = load_comparisons()
    #go through the values of the dictionary and break when we find the right groups in that contrast
    grps, comparison_column = return_sample_names_group(grp)
    #take the sample names corresponding to those groups
    if comparison_column == "":
        return([""])
    grp_samples = list(set(list(samples2[samples2[comparison_column].isin(grps)].sample_name)))
    salmon_outdir = get_output_dir(config["project_top_level"], config["salmon_output_folder"])
    salmon_suffix = "quant.sf"

    #build a list with the full path from those sample names
    fc_files = [os.path.join(salmon_outdir, x , salmon_suffix) \
                   for x in grp_samples]
    fc_files = list(set(fc_files))
    print(fc_files)
    return(fc_files)


def load_comparisons():
    comparisons = "config/DESeq2comparisons.yaml"
    with open(comparisons, 'r') as stream:
        try:
            compare_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return(compare_dict)


def check_key(dict, key):
     """
     simple funciton to check if a key is in a dictionary, if it's not returns false and if it is returns the value
     """
     if key in dict:
         return(dict[key])
     else:
         return(False)


def return_sample_names_group(grp):
    """
    given a group, return the names and the column_name associated with that
    """
    compare_dict = load_comparisons()
    for key in compare_dict.keys():
        grp_names = check_key(compare_dict[key],grp)
        if grp_names:
            column_name = compare_dict[key]['column_name'][0]
            return(grp_names,column_name)
    return("","")
