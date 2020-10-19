rule subsampleReadsIllumina_PE_test_depth:
	input:
		paired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.tot.txt",),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq"),
	output:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_clean.tot.fastq"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_clean.tot.fastq"),
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
rule normalizeReads_PE_test_depth:
	input:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_clean.tot.fastq"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_clean.tot.fastq"),
	output:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_norm.{sampling}.fastq"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_norm.{sampling}.fastq"),
	message:
		"Normalizing reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	params:
		min_depth=config['min_norm'],
		max_depth=config['max_norm']
	threads: 4
	resources:
		mem_mb=6000
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
		target={params.max_depth} mindepth={params.min_depth} t={threads}
		"""
rule shortReadAsemblySpadesPE__test_depth:
	input:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_norm.tot.fastq"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_norm.tot.fastq"),
	output:
		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_filtered_scaffolds.{sampling}.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}/filtered_list.txt"),
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_{sampling}"),
	message:
		"Assembling PE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired} -o {params.assembly_dir} \
		--meta -t {threads} --only-assembler --memory 350
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		"""

rule virSorter_test_depth:
	input:
		merged_assembly=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_spades_filtered_scaffolds.{sampling}.fasta"),
		virSorter_db=config['virSorter_db'],
	output:
		positive_fasta=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/final-viral-combined.fa",
		table_virsorter=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/final-viral-score.tsv",
		positive_list=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/positive_VS_list_{sampling}.txt",
		viral_boundary=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/final-viral-boundary.tsv",
	params:
		out_folder=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 32
	shell:
		"""
		virsorter run -w {params.out_folder} -i {input.merged_assembly} -j {threads}
		grep ">" {output.positive_fasta} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list}
		"""

rule estimateGenomeCompletness_test_depth:
	input:
		positive_fasta=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_virSorter_{sampling}/final-viral-combined.fa",
	output:
		quality_summary=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_checkV_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_checkV_{sampling}/completeness.tsv",
		contamination=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_checkV_{sampling}/contamination.tsv",
		repeats=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_checkV_{sampling}/repeats.tsv",
	params:
		checkv_outdir=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_checkV_{sampling}",
		checkv_db=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_checkV_{sampling}",

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
