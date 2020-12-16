rule virSorter:
	input:
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
		virSorter_db=config['virSorter_db'],
	output:
		positive_fasta=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/final-viral-combined.fa",
		table_virsorter=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/final-viral-score.tsv",
		positive_list=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/positive_VS_list_{sampling}.txt",
		viral_boundary=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/final-viral-boundary.tsv",
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		out_folder=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 32
	shell:
		"""
		virsorter run -w {params.out_folder} -i {input.merged_assembly} -j {threads}
		grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list}
		cp {output.positive_fasta} {output.positive_contigs}
		"""

# rule extractViralContigs:
# 	input:
# 		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
# 		positive_list=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/positive_VS_list_{sampling}.txt",
# 	output:
# 		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
# 	message:
# 		"Selecting Viral Contigs"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/vir.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		seqtk subseq {input.merged_assembly} {input.positive_list} > {output.positive_contigs}
# 		"""

rule whatThePhage:
	input:
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
		WTP_dir=config['WTP_dir'],
	output:
		positive_fasta=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/final-viral-combined.fa",
		table_virsorter=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/final-viral-score.tsv",
		positive_list=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/positive_VS_list_{sampling}.txt",
		viral_boundary=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/final-viral-boundary.tsv",
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		out_folder=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 32
	shell:
		"""
		virsorter run -w {params.out_folder} -i {input.merged_assembly} -j {threads}
		grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list}
		cp {output.positive_fasta}
		"""
