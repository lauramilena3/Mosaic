rule subsampleReadsIllumina_PE:
	input:
		paired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.tot.txt",),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq"),
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_{subsample}_forward_paired_clean.sub.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_{subsample}_reverse_paired_clean.sub.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_{subsample}_unpaired_clean.sub.fastq",
		paired_size=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_{subsample}_paired_clean.sub.txt"),
		unpaired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_{subsample}_unpaired_clean.sub.txt"
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		#paired
		paired=$( cat {input.paired_sizes} )
		p=$(echo "$paired"*{wildcards.subsample}/100 | bc)
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads=$p
		"""
    
rule estimateGenomeCompletness:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		quality_summary=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/completeness.tsv",
		contamination=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/contamination.tsv",
		repeats=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/repeats.tsv",
	params:
		checkv_outdir=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}",
		checkv_db=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}",

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
