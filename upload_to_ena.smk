import os
# a top level folder where the bams reside
directory_with_bams_and_manifests = "/SAN/vyplab/alb_projects/data/ward_bams/bam_files"

bam_suffix = ".pass2Aligned.sortedByCoord.out.bam"

# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------
singularity_path = directory_with_bams_and_manifests.replace("/SAN/vyplab/alb_projects/data", "/home/alb_data")
SAMPLES, = glob_wildcards(directory_with_bams_and_manifests + "{sample}" + bam_suffix)
print(SAMPLES)
print(directory_with_bams_and_manifests)
print(singularity_path)

rule all:
  input:
    #expand(output_dir + "{sample}_namesorted.bam", sample = SAMPLES),
    expand(directory_with_bams_and_manifests + "{sample}_uploaded", sample = SAMPLES)

rule upload_try:
    input:
        manifest = directory_with_bams_and_manifests + "{sample}" + ".manifest"
    output:
       directory_with_bams_and_manifests + "{sample}_uploaded"
    params:
        singularity_path_manifest = singularity_path + "{sample}" + ".manifest"
    shell:
        """
        singularity run --bind /SAN/vyplab/alb_projects/data/:/home/alb_data /SAN/vyplab/alb_projects/tools/webin-cli_latest.sif -context reads -manifest {params.singularity_path_manifest} -userName 'Webin-58069' -password '1HIkW7lgrw' -submit
        touch {output}
        """
