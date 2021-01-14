import re
import pandas as pd
import os
import numpy as np
import subprocess
from subprocess import PIPE

def get_fastq_names(DATA):
    samples = pd.read_csv(DATA, sep = ",")
    #take all the fasqs and combine them to a list remove any nas
    fastq_list = samples.fast1.dropna().tolist() + samples.fast2.dropna().tolist()
    #if there are any missing values pandas gets annoyed so replace nans with empty string
    samples = samples.replace(np.nan, '', regex=True)
    #get the associated sample and unit name for each fastq using and or operator to get either fast1 or fast2, works for both single-end and paired end this way
    unit_name = [samples.loc[(samples['fast1'] == fq)| (samples['fast2'] == fq)].unit.iloc[0] for fq in fastq_list]
    #strip it down to just the name of the file bit
    stripped = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in fastq_list]
    #now combine these three lists to get a meaningfull and unique string for each fastq file

    return(stripped, fastq_list, unit_name)

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

    trimmed_1 = [os.path.join(fastp_outdir,u + "_R1_trimmed.fastq.gz") for u in unit_fastqs2]


    #if we have paired end data there will also be a trimmed 2, same thing, using the fast2 column instead
    if config['end_type'] == "pe":
        trimmed_2 = [os.path.join(fastp_outdir,u + "_R2_trimmed.fastq.gz") for u in unit_fastqs2]

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

# takes the featurcounts strand and returns the interpretation for kallisto_output_folder
def get_kallisto_strand(fcStrand):
    if fcStrand == "-s 0":
        return("")
    elif fcStrand == "-s 1":
        return("--fr-stranded")
    elif fcStrand == "-s 2":
        return("--rf-stranded")


def get_salmon_strand(fcStrand):
    '''
    Return string for salmon libtype denoting strandedness/orientation of library.
    Parses featureCounts strand and returns corresponding salmon libtype string
    (include option to return the let Salmon infer for you string)
    '''
    return "placeholder"


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
