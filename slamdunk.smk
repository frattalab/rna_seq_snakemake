import os
output_folder = "/SAN/vyplab/alb_projects/data/4su_tdp_f210i/slam_dunk_analysis/full_run/"
samples = os.listdir("/SAN/vyplab/alb_projects/data/4su_tdp_f210i/merged_fastqs/")
sample_clean = [te[:-16] for te in samples]
rule final:
	input:
		expand(output_folder + "count/" + "{sample}.merged.fastq_slamdunk_mapped_filtered_tcount.tsv", sample = sample_clean)

rule run_slam_dunk:
	input:
		"/SAN/vyplab/alb_projects/data/4su_tdp_f210i/merged_fastqs/{sample}.merged.fastq.gz"
	output:
		output_folder + "count/" + "{sample}.merged.fastq_slamdunk_mapped_filtered_tcount.tsv"
	shell:
		"/SAN/vyplab/alb_projects/tools/slamdunk/bin/slamdunk all -r /SAN/vyplab/vyplab_reference_genomes/sequence/mouse/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa -b /SAN/vyplab/vyplab_reference_genomes/annotation/mouse/GRCm38/Mus_musculus.GRCm38.96.bed -o {output_folder} -rl 75 -5 12 -t 4 -m --skip-sam {input}"