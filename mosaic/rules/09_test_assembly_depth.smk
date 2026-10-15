rule subsampleReadsIllumina_PE_test_depth:
	input:
		paired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}_read_count.txt",),
		unpaired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}_read_count.txt",),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq.gz"),
	output:
		forward_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_unpaired_clean.{sampling}.fastq.gz"),
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/subsampleReadsIllumina_PE_test_depth/sample={sample}__sampling={sampling}__subsample={subsample}.tsv"
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		#paired
		paired=$( cat {input.paired_sizes} )
		p=$(echo "$paired"*{wildcards.subsample}/100 | bc)
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads=$p sampleseed=1
		#unpaired
		unpaired=$( cat {input.unpaired_sizes} )
		up=$(echo "$unpaired"*{wildcards.subsample}/100 | bc)
		reformat.sh in={input.unpaired} out={output.unpaired} reads=$up sampleseed=1
		"""

rule normalizeReads_test_depth:
	input:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_unpaired_clean.{sampling}.fastq.gz"),
	output:
		forward_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_unpaired_norm.{sampling}.fastq.gz"),
	message:
		"Normalizing reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalizeReads_test_depth/sample={sample}__sampling={sampling}__subsample={subsample}.tsv"
	params:
		min_depth=config['min_norm'],
		max_depth=config['max_norm'],
		heap_mb=MEMORY_ECORR
	threads: 4
	resources:
		mem_mb=BBTOOLS_MEM_MB
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{params.heap_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
		target={params.max_depth} mindepth={params.min_depth} t={threads} 
		#unpaired
		bbnorm.sh -Xmx{params.heap_mb}m in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth} t={threads}
		"""

rule metaspadesPE_test_depth:
	input:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_unpaired_norm.{sampling}.fastq.gz"),
	output:
		scaffolds=temp(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_filtered_scaffolds.{sampling}.unrenamed.fasta"),
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_{sampling}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_{sampling}"),
		filtered_list=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_{sampling}/filtered_list.txt"),
		mem_gb=ASSEMBLY_MEM_GB,
	message:
		"Assembling PE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/metaspadesPE_test_depth/sample={sample}__sampling={sampling}__subsample={subsample}.tsv"
	threads: 8
	resources:
		mem_mb=ASSEMBLY_MEM_MB
	shell:
		"""
		rm -rf {params.assembly_dir}
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		--meta -t {threads} --memory {params.mem_gb}
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {params.filtered_list}
		seqtk subseq {params.raw_scaffolds} {params.filtered_list} > {output.scaffolds}
		sed "s/>/>{wildcards.sample}_{wildcards.subsample}_/g" -i {output.scaffolds}
		rm -rf {params.assembly_dir}
		"""

rule assemblyStatsILLUMINA_test_depth:
	input:
		scaffolds_assembly=expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES, subsample=subsample_test),
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/assembly_statistics_quast_{sampling}"),
		quast_txt=dirs_dict["ASSEMBLY_TEST"] + "/assembly_quast_report.{sampling}.txt",
		quast_tsv=dirs_dict["ASSEMBLY_TEST"] + "/assembly_statistics_quast_{sampling}/transposed_report.tsv",
	message:
		"Creating assembly-depth statistics with QUAST"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assemblyStatsILLUMINA_test_depth/sampling={sampling}.tsv"
	log:
		dirs_dict["ASSEMBLY_TEST"] + "/quast.{sampling}.log"
	threads: 4
	resources:
		mem_mb=QUAST_DEPTH_MEM_MB
	shell:
		"""
		quast_fastas=()
		quast_labels=""
		for quast_fasta in {input.scaffolds_assembly:q}; do
			if [ -s "$quast_fasta" ]; then
				quast_fastas+=("$quast_fasta")
				quast_label=$(basename "$quast_fasta" .fasta)
				quast_labels+="${{quast_label}},"
			fi
		done
		if [ "${{#quast_fastas[@]}}" -gt 0 ]; then
			quast.py "${{quast_fastas[@]}}" -o {output.quast_report_dir:q} \
				--labels "${{quast_labels%,}}" --min-contig 1 --threads {threads} > {log:q} 2>&1
			cp {output.quast_report_dir:q}/report.txt {output.quast_txt:q}
		else
			mkdir -p {output.quast_report_dir:q}
			printf 'Assembly\n' > {output.quast_tsv:q}
			printf 'No contigs to assess.\n' > {output.quast_report_dir:q}/report.txt
			cp {output.quast_report_dir:q}/report.txt {output.quast_txt:q}
			printf 'No contigs to assess.\n' > {log:q}
		fi
		"""

rule genomad_viral_id_subassembly:
	input:
		scaffolds=(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_filtered_scaffolds.{sampling}.fasta"),
		genomad_db=(config['genomad_db']),
	output:
		genomad_outdir=directory(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_geNomad_{sampling}/"),
		final_viral_contigs=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_positive_geNomad.{sampling}.fasta",
	params:
		viral_fasta=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_geNomad_{sampling}/{sample}_{subsample}_metaspades_filtered_scaffolds.{sampling}_summary/{sample}_{subsample}_metaspades_filtered_scaffolds.{sampling}_virus.fna",
	message:
		"Identifying viral contigs with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/genomad_viral_id_subassembly/sample={sample}__sampling={sampling}__subsample={subsample}.tsv"
	threads: 8
	shell:
		"""
		if [ -s {input.scaffolds} ]; then
				genomad end-to-end --cleanup --splits 8 -t {threads} {input.scaffolds} {output.genomad_outdir} {input.genomad_db} --relaxed
				seqtk seq -l 0 {params.viral_fasta} | sed "s/|/_/g" > {output.final_viral_contigs}
		else
				echo "The FASTA file {input.scaffolds} is empty"
				mkdir -p {output.genomad_outdir}
				touch {output.final_viral_contigs}
		fi
		"""

rule viralStatsILLUMINA_test_depth:
	input:
		scaffolds_viral=expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_positive_geNomad.{{sampling}}.fasta", sample=SAMPLES, subsample=subsample_test),
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_TEST"] + "/assembly_statistics_viral_contigs_quast_{sampling}"),
		quast_txt=dirs_dict["ASSEMBLY_TEST"] + "/assembly_quast_report_viral.{sampling}.txt",
		quast_tsv=dirs_dict["ASSEMBLY_TEST"] + "/assembly_statistics_viral_contigs_quast_{sampling}/transposed_report.tsv",
	message:
		"Creating viral assembly-depth statistics with QUAST"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/viralStatsILLUMINA_test_depth/sampling={sampling}.tsv"
	log:
		dirs_dict["ASSEMBLY_TEST"] + "/quast_viral.{sampling}.log"
	threads: 4
	resources:
		mem_mb=QUAST_DEPTH_MEM_MB
	shell:
		"""
		quast_fastas=()
		quast_labels=""
		for quast_fasta in {input.scaffolds_viral:q}; do
			if [ -s "$quast_fasta" ]; then
				quast_fastas+=("$quast_fasta")
				quast_label=$(basename "$quast_fasta" .fasta)
				quast_labels+="${{quast_label}},"
			fi
		done
		if [ "${{#quast_fastas[@]}}" -gt 0 ]; then
			quast.py "${{quast_fastas[@]}}" -o {output.quast_report_dir:q} \
				--labels "${{quast_labels%,}}" --min-contig 1 --threads {threads} > {log:q} 2>&1
			cp {output.quast_report_dir:q}/report.txt {output.quast_txt:q}
		else
			mkdir -p {output.quast_report_dir:q}
			printf 'Assembly\n' > {output.quast_tsv:q}
			printf 'No viral contigs to assess.\n' > {output.quast_report_dir:q}/report.txt
			cp {output.quast_report_dir:q}/report.txt {output.quast_txt:q}
			printf 'No viral contigs to assess.\n' > {log:q}
		fi
		"""

rule estimateGenomeCompletness_test_depth:
	input:
		final_viral_contigs=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_positive_{viral_id_tool}.{sampling}.fasta",
	output:
		quality_summary=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_{viral_id_tool}_checkV_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_{viral_id_tool}_checkV_{sampling}/completeness.tsv",
		contamination=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_{viral_id_tool}_checkV_{sampling}/contamination.tsv",
	params:
		checkv_outdir=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_{viral_id_tool}_checkV_{sampling}",
		tmp=dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_{viral_id_tool}_checkV_{sampling}/tmp",
#		checkv_db=dirs_dict["ASSEMBLY_TEST"] + "/95-80_merged_positive_virsorter_checkV_{sampling}",
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/estimateGenomeCompletness_test_depth/sample={sample}__sampling={sampling}__subsample={subsample}__viral_id_tool={viral_id_tool}.tsv"
	threads: 4
	shell:
		"""
		if [ -s {input.final_viral_contigs} ]; then
				rm -rf {params.checkv_outdir} || true
				checkv contamination {input.final_viral_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
				checkv completeness {input.final_viral_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
				checkv complete_genomes {input.final_viral_contigs} {params.checkv_outdir}
				checkv quality_summary {input.final_viral_contigs} {params.checkv_outdir}
				rm -rf {params.tmp}
		else
				echo "The FASTA file {input.final_viral_contigs} is empty"
				mkdir -p {params.checkv_outdir}
				touch {output.quality_summary}
				touch {output.completeness}
				touch {output.contamination}
		fi

		"""
