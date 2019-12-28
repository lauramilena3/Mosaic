ruleorder: shortReadAsemblySpadesPE > shortReadAsemblySpadesSE


rule shortReadAsemblySpadesPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq",
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/filtered_list.txt"),
		assembly_graph=dirs_dict["ASSEMBLY_DIR"] +"/{sample}_assembly_graph_spades.{sampling}.fastg",
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/scaffolds.fasta",
		assembly_graph=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/assembly_graph.fastg",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}"),
	message:
		"Assembling PE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 16
	shell:
		"""
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		--meta -t {threads} --only-assembler
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		assembly_graph.fastg
		"""
rule shortReadAsemblySpadesSE:
	input:
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq"
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/filtered_list.txt")
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}")
	message:
		"Assembling SE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 16
	shell:
		"""
		spades.py -s {input.unpaired} -o {params.assembly_dir} \
		--sc -t {threads} --only-assembler
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		"""
rule assemblyStatsILLUMINA:
	input:
		quast_dir=directory(config["quast_dir"]),
		scaffolds_spades=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES)
	output:
		quast_report_dir=(dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}"),
		quast_txt=dirs_dict["ASSEMBLY_DIR"] + "/assembly_quast_report.{sampling}.txt"
	message:
		"Creating assembly stats with quast"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		{input.quast_dir}/quast.py {input.scaffolds_spades} -o {output.quast_report_dir} --threads {threads}
		cp {output.quast_report_dir}/report.txt {output.quast_txt}
		"""

rule mergeAssembliesSHORT:
	input:
		scaffolds_spades=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta",sample=SAMPLES)
	output:
		merged_assembly=(dirs_dict["vOUT_DIR"] + "/merged_scaffolds.{sampling}.fasta")
	message:
		"Merging assembled contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		cat {input.scaffolds_spades} > {output.merged_assembly}
		"""
