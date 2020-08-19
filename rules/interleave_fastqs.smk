import os
import pandas as pd

configfile: "config/config.yaml"
include: "helpers.py"


sample_tbl = pd.read_csv(config["sampleCSVpath"], sep = ",")
sample_tbl = sample_tbl.replace(np.nan, '', regex=True)

SAMPLE_NAMES = sample_tbl['sample_name'].tolist()
UNITS = SAMPLES['unit'].tolist()

#make sure master output_dir & subdirectories are created prior to running
#(save typing later on000)
temp_dir = os.path.join(config['interleave_master_output_folder'],"temp",'')
singletons_dir = os.path.join(config['interleave_master_output_folder'],"singletons",'')

os.system("mkdir -p {0}".format(config['interleave_master_output_folder']))
os.system("mkdir -p {0}".format(temp_dir))
os.system("mkdir -p {0}".format(singletons_dir))


rule interleave_all_merged:
    input:
        expand(config['merged_fastq_folder'] + "{name}_1.merged.fastq.gz", name = SAMPLE_NAMES),
        expand(config["merged_fastq_folder"] + "{name}_2.merged.fastq.gz" if config["end_type"] == "pe" else [],name = SAMPLE_NAMES),
        expand("{interleaved_outdir}{name}_interleaved.fastq.gz", interleaved_outdir = os.path.join(config['interleave_master_output_folder'],''),name = SAMPLE_NAMES),
        expand("{singletons_subdir}{name}_singletons.fastq.gz", singletons_subdir = singletons_dir, name = SAMPLE_NAMES)

#First run bbmap repair.sh to make sure readA.1 is first in fq.1 & readA.2 is first in fq.2 etc

rule validate_fastq_pairing:
    input:
        fastq_file = os.path.join(config['merged_fastq_folder'],"{name}_1.merged.fastq.gz"),
        fastq_file2 = os.path.join(config['merged_fastq_folder'],"{name}_2.merged.fastq.gz"),
        #fastq_file = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True)
        #fastq_file2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False)

    output:
        paired_fastq_file = temp(os.path.join(temp_dir,"{name}_fixed_1.fastq.gz")),
        paired_fastq_file2 = temp(os.path.join(temp_dir,"{name}_fixed_2.fastq.gz")),
        singletons = os.path.join(singletons_dir,"{name}_singletons.fastq.gz")

    params:
        repair_script = os.path.join(config['bbmap_path'],"repair.sh")

    threads : 4

    shell:
        '''
        bash {params.repair_script} in1={input.fastq_file} in2={input.fastq_file2} \
        out1={output.paired_fastq_file} out2={output.paired_fastq_file2} \
        outs={output.singletons} \
        t={threads} repair
        '''

rule interleave_fastqs:
    input:
        fastq_file = os.path.join(temp_dir,"{name}_fixed_1.fastq.gz"),
        fastq_file2 = os.path.join(temp_dir,"{name}_fixed_2.fastq.gz")

    output:
        os.path.join(config['interleave_master_output_folder'],"{name}_interleaved.fastq.gz")

    params:
        reformat_script = os.path.join(config['bbmap_path'],"reformat.sh"),

    threads: 4

    shell:
        '''
        bash {params.reformat_script} in1={input.fastq_file} in2={input.fastq_file2} \
        out={output} t={threads}
        '''
