ruleorder: mapReadsToContigsPE > mapReadsToContigsSE
rule createContigBowtieDb:
	input:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{type}.fasta",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{type}.fasta"
	output:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{type}.1.bt2",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{type}.1.bt2"
	params:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{type}",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{type}"
	message:
		"Selecting Viral Contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bowtie2-build -f {input.high_contigs} {params.high_contigs}
		bowtie2-build -f {input.low_contigs} {params.low_contigs}
		"""

rule mapReadsToContigsPE:
	input:
		high_bt2=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{type}.1.bt2",
		low_bt2=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{type}.1.bt2",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{type}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{type}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{type}.fastq",
	output:
		high_sam=dirs_dict["VIRAL_DIR"]+ "/{sample}_high_confidence.sam",
		low_sam=dirs_dict["VIRAL_DIR"]+ "/{sample}_low_confidence.sam",
	params:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{type}",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{type}"
	message:
		"Selecting Viral Contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bowtie2 --non-deterministic -x {params.high_contigs} -1 {input.forward_paired} -2 {input.reverse_paired} \
		-U {input.unpaired} -S {output.high_sam} 
		bowtie2 --non-deterministic -x {params.low_contigs} -1 {input.forward_paired} -2 {input.reverse_paired} \
		-U {input.unpaired} -S {output.low_sam} 
		"""
rule mapReadsToContigsSE:
	input:
		high_bt2=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{type}.1.bt2",
		low_bt2=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{type}.1.bt2",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{type}.fastq",
	output:
		high_sam=dirs_dict["VIRAL_DIR"]+ "/{sample}_high_confidence.sam",
		low_sam=dirs_dict["VIRAL_DIR"]+ "/{sample}_low_confidence.sam",
	params:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{type}",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{type}"
	message:
		"Selecting Viral Contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bowtie2 --non-deterministic -x {params.high_contigs} -U {input.unpaired} -S {output.high_sam} 
		bowtie2 --non-deterministic -x {params.low_contigs} -U {input.unpaired} -S {output.low_sam} 
		"""
