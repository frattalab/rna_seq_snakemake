
configfile: "config/config.yaml"
include: "helpers.py"

#get output folders
fastp_outdir = get_output_dir(config['project_top_level'], config["fastp_trimmed_output_folder"])
merged_outdir = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

#make sure the output folder for fastptrimming exists before running anything
os.system("mkdir -p {0}".format(fastp_outdir))
#read in a samples table
SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)
#so I want this rule to be run ONCE for every fast1, so the wild cards I'm giving are the 'name' of the fastq of the first read
#name = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in SAMPLES['fast1'].tolist()]


UNITS = SAMPLES['unit'].tolist()
SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

rule merge_all_trimmed:
    input:
        expand(fastp_outdir + "{unit}_{name}_1.trimmed.fastq.gz",zip, unit = UNITS,name = SAMPLE_NAMES),
        expand(fastp_outdir + "{unit}_{name}_2.trimmed.fastq.gz" if config["end_type"] == "pe" else [],zip, unit = UNITS,name = SAMPLE_NAMES),
        expand(merged_outdir + "{name}_1.merged.fastq.gz", name = SAMPLE_NAMES),
        expand(merged_outdir + "{name}_2.merged.fastq.gz" if config["end_type"] == "pe" else [],name = SAMPLE_NAMES)
    wildcard_constraints:
        name="|".join(SAMPLE_NAMES),
        unit="|".join(UNITS)

if config['end_type'] == "pe":
    rule fastp_trimming:
        input:
        #get the value in the fast1 column
            fastq_file = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True),
            fastq_file2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False)
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES),
            unit="|".join(UNITS)
        conda:
            "../env/align.yaml"
        output:
            out_fastqc = fastp_outdir + "{unit}_{name}_1.trimmed.fastq.gz",
            out_fastqc2 = fastp_outdir + "{unit}_{name}_2.trimmed.fastq.gz",
            fastpjson = fastp_outdir + "{unit}_{name}_fastp.json",
            fastphtml = fastp_outdir + "{unit}_{name}_fastp.html",
        params:
            fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
            fastpjson = fastp_outdir + "{unit}_{name}_fastp.json",
            fastphtml = fastp_outdir + "{unit}_{name}_fastp.html"

        log:
            os.path.join(config['project_top_level'], "logs", "{unit}_{name}.fastp_stdout.log")

        shell:
            """
            #free -h
            fastp \
            --in1 {input.fastq_file} \
            --in2 {input.fastq_file2} \
            --out1 {output.out_fastqc} \
            --out2 {output.out_fastqc2} \
            --json {output.fastpjson} \
            --html {output.fastphtml} \
            {params.fastp_parameters} \
            2> {log}
            """
else:
        rule fastp_trimming:
            input:
            #get the value in the fast1 column
                fastq_file = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True)
            wildcard_constraints:
                name="|".join(SAMPLE_NAMES),
                unit="|".join(UNITS)
            conda:
                "../env/align.yaml"
            output:
                out_fastqc = fastp_outdir + "{unit}_{name}_1.trimmed.fastq.gz",
                fastpjson = fastp_outdir + "{unit}_{name}_fastp.json",
                fastphtml = fastp_outdir + "{unit}_{name}_fastp.html",
            params:
                fastp_parameters = return_parsed_extra_params(config['fastp_parameters']),
                fastq_file2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False),
            #out_fastqc2 = lambda wildcards: return_fastq2_name(wildcards.name,wildcards.unit),
                fastpjson = fastp_outdir + "{unit}_{name}_fastp.json",
                fastphtml = fastp_outdir + "{unit}_{name}_fastp.html"

            log:
                os.path.join(config['project_top_level'], "logs", "{unit}_{name}.fastp_stdout.log")

            shell:
                """
                fastp -i {input.fastq_file} \
                -o {output.out_fastqc} \
                --json {output.fastpjson} \
                --html {output.fastphtml} \
                {params.fastp_parameters} \
                2> {log}
                """

if config['end_type'] == "pe":
    rule merge_trimmed:
        input:
            one = lambda wildcards: get_trimmed(wildcards.name)[0],
            two = lambda wildcards: get_trimmed(wildcards.name)[1]
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES)
        output:
            out_one = merged_outdir + "{name}_1.merged.fastq.gz",
            out_two = merged_outdir + "{name}_2.merged.fastq.gz"
        params:
            #taking the input files and putting them into a comma separated list
            one = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[0]),
            two = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[1])
        shell:
            """
            cat {params.one} > {output.out_one}
            cat {params.two} > {output.out_two}
            """
else:
    rule merge_trimmed:
        input:
            one = lambda wildcards: get_trimmed(wildcards.name)[0]
        wildcard_constraints:
            name="|".join(SAMPLE_NAMES)
        output:
            out_one = merged_outdir + "{name}_1.merged.fastq.gz"
        params:
            #taking the input files and putting them into a comma separated list
            one = lambda wildcards: ' '.join(get_trimmed(wildcards.name)[0]),
        shell:
            """
            cat {params.one} > {output.out_one}
            """
