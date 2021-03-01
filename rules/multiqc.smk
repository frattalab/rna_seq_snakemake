import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

# Need units from here for fastq files
FASTQ_NAME, FILE_LOCATION, UNITS = get_fastq_names(config["sampleCSVpath"])

## This loop checks if workflow has already been defined (i.e. running inside workflow)
## if not then not defined (i.e. if just want to run multiqc), so give it a dummy value

try:
    workflow_str
except NameError:
    # Wasn't defined, give a dummy value
    workflow_str = "multiqc_output"
else:
    # Already defined, no need to change
    pass


multiqc_output_folder = os.path.join(get_output_dir(config["project_top_level"], config["multiqc_output_folder"]), workflow_str, "")

rule all:
    input:
        multiqc_output_folder + "multiqc_report.html"


rule multiqc:
    input:
        multiqc_target_files(workflow, SAMPLE_NAMES, UNITS)

    output:
        os.path.join(multiqc_output_folder, "multiqc_report.html")

    params:
        dirs = multiqc_target_dirs(),
        outdir = multiqc_output_folder,
        save_plots = "-p", # Should provide option to alter this in config
        runtime = "--profile-runtime"

    conda:
        "../env/align.yaml"

    shell:
        """
        multiqc \
        {params.save_plots} \
        -o {params.outdir} \
        {params.runtime} \
        {params.dirs}
        """
