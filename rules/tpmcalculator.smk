import os
import subprocess

include: "helpers.py"
configfile: "config/config.yaml"

singlestep = "true"

if singlestep == "true":
    project_folder =  "/SAN/vyplab/alb_projects/data/tdp_ko_collection/"
    end_type = "pe"
    suffix = ".Aligned.sorted.out"
    star_output_folder = project_folder + "all_bams_sym"
    SAMPLE_NAMES, = glob_wildcards(star_output_folder + "{sample}" + suffix + ".bam")

    tpm_output_folder = project_folder + "TPMcalculator"

    print(SAMPLE_NAMES)
else:
    project_folder = config["project_top_level"]
    suffix = ".Aligned.sorted.out"
    SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
    SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

    SAMPLE_NAMES = SAMPLES['sample_name'].tolist()
    #make sure the output folder for featureCounts exists before running anything
    star_output_folder = config['star_output_folder']
    end_type = config["end_type"]
    tpm_output_folder = project_folder + "TPMcalculator"
#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species
REFERENCE_ANNOTATION = get_gtf(config['species'])
os.system("mkdir -p {0}".format("tpm_output_folder"))


rule tpmcounts:
    input:
        expand(tpm_output_folder + "{name}" + suffix + "_genes.out", name = SAMPLE_NAMES)


rule tpmcalculator_path:
    input:
        aligned_bam = star_output_folder + "{name}.Aligned.sorted.out.bam",
        aligned_bai = star_output_folder + "{name}.Aligned.sorted.out.bam.bai"
    output:
        tpm_output_folder + "{name}" + suffix + "_genes.out"
    params:
        ref_anno = REFERENCE_ANNOTATION
    run:
        shell("echo moving into {tpm_output_folder}")
        shell("cd {tpm_output_folder}")
        if config["end_type"] == "pe":
            shell("{config[tpmcalculator_path]} -g {params.ref_anno} -b {input.aligned_bam} -p -e -a")
            subprocess.run(["ls", "-l"])
        if config["end_type"] == "se":
            shell("{config[tpmcalculator_path]} -g {params.ref_anno} -b {input.aligned_bam} -e -a")
