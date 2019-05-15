ruleorder: mergeAssembliesHIBRID > mergeAssembliesSHORT

rule mergeAssembliesHIBRID:
	input:
		scaffolds_spades=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta",sample=SAMPLES),
		scaffolds_canu=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES),
	output:
		merged_assembly=(dirs_dict["vOUT_DIR"] + "/merged_scaffolds.{sampling}.fasta")
	message:
		"Merging assembled contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		cat {input.scaffolds_canu} {input.scaffolds_spades} > {output.merged_assembly}
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

rule vOUTclustering:
	input:
		merged_assembly=(dirs_dict["vOUT_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
	output:
		clusters=dirs_dict["vOUT_DIR"] + "/merged_scaffolds.{sampling}_95-80.clstr",
		representatives=dirs_dict["vOUT_DIR"] + "/merged_scaffolds.{sampling}_95-80.fna",
		representative_lengths=dirs_dict["vOUT_DIR"] + "/representative_lengths.{sampling}.txt"
	message:
		"Creating vOUTs with stampede"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		./scripts/stampede-Cluster_genomes.pl -f {input.merged_assembly} -c 80 -i 95
		cat {output.representatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
		$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.representative_lengths}
		"""

#rule circularizeContigs:
#	input:
#		seeds=dirs_dict["vOUT_DIR"] + "/{sample}_merged_scaffolds_95-80.fna",
#		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.fastq"),
#		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.fastq"),
#		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.fastq"
#	output:
#		circular_contigs=(dirs_dict["vOUT_DIR"] + "/{sample}_merged_scaffolds.fasta"),
#		circlator_dir=dirs_dict["vOUT_DIR"] + "/{sample}_circlator",
#		merged_reads=dirs_dict["vOUT_DIR"] + "/{sample}_merged_reads.fastq"
#		06.fixstart.fasta
#	message:
#		"Circularizing representative clusters with Circlator"
#	conda:
#		dirs_dict["ENVS_DIR"] + "/env1.yaml"
#	threads: 4
#	shell:
#		"""
#		cat {input.forward_paired} {input.reverse_paired} {input.unpaired} > {output.merged_reads}
#		circlator all --threads {threads} {input.seeds} {output.merged_reads} {output.circlator_dir}
#		"""