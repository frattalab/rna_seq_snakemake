import os
# a top level folder where the bams reside
directory_with_bams_and_manifests = "/SAN/vyplab/alb_projects/data/buratti_new_shsy5y/STAR_aligned/"

bam_suffix = ".Aligned.sorted.out.bam"

# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------
singularity_path = txt.replace("/SAN/vyplab/alb_projects/data", "/home/alb_data")
SAMPLES, = glob_wildcards(directory_with_bams_and_manifests + "{sample}" + bam_suffix)
print(SAMPLES)

rule all:
  input:
    #expand(output_dir + "{sample}_namesorted.bam", sample = SAMPLES),
    expand(directory_with_bams_and_manifests  + "{sample}_uploaded", sample = SAMPLES)

rule upload_try:
    input:
        manifest = singularity_path + "{sample}" + ".manifest"
    output:
       directory_with_bams_and_manifests + "{sample}_uploaded"
    shell:
        """
        singularity run --bind /SAN/vyplab/alb_projects/data/:/home/alb_data /SAN/vyplab/alb_projects/tools/webin-cli_latest.sif -context reads -manifest {input.manifest} -userName 'Webin-58069' -password '1HIkW7lgrw' -submit
        touch {output}
        """
