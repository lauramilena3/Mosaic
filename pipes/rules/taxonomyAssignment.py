rule getORFs:
	input:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.fasta",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.fasta"
	output:
		high_coords=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.coords",
		low_coords=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.coords",
		high_aa=dirs_dict["VIRAL_DIR"]+ "/high_confidence_ORFs.{sampling}.fasta",
		low_aa_info=dirs_dict["VIRAL_DIR"]+ "/low_confidence_ORFs.{sampling}.fasta",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		prodigal -i {input.high_contigs} -o {output.high_coords} -a {output.high_aa} -p meta
		prodigal -i {input.high_contigs} -o {output.high_coords} -a {output.high_aa} -p meta
		"""