rule vOUTclustering:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		clusters=dirs_dict["VIRAL_DIR"] + "/"+ VIRAL_CONTIGS_BASE + ".{sampling}_95-80.clstr",
		representatives_temp=temp(dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}_95-80.fna"),
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		representative_lengths=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_lengths.{sampling}.txt",
	message:
		"Creating vOUTs with stampede"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		./scripts/stampede-Cluster_genomes.pl -f {input.positive_contigs} -c 80 -i 95
		cat {output.representatives_temp} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
		$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.representative_lengths}
		cp {output.representatives_temp} {output.representatives}
		"""
rule estimateGenomeCompletness:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		quality_summary=dirs_dict["VIRAL_DIR"] + "/checkV_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["VIRAL_DIR"] + "/checkV_{sampling}/completeness.tsv",
		contamination=dirs_dict["VIRAL_DIR"] + "/checkV_{sampling}/contamination.tsv",
		repeats=dirs_dict["VIRAL_DIR"] + "/checkV_{sampling}/repeats.tsv",
	params:
		checkv_outdir=dirs_dict["VIRAL_DIR"] + "/checkV_{sampling}",
		checkv_db=dirs_dict["VIRAL_DIR"] + "/checkV_{sampling}",

	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 4
	shell:
		"""
		checkv contamination {input.representatives} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv completeness {input.representatives} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv repeats {input.representatives} {params.checkv_outdir}
		checkv quality_summary {input.representatives} {params.checkv_outdir}
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
