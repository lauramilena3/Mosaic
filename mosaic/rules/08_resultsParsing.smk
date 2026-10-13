def benchmark_snapshot(wildcards):
	from pathlib import Path
	root=Path(dirs_dict["BENCHMARKS"])
	return [(str(path.relative_to(root)), path.stat().st_size, path.stat().st_mtime_ns)
		for path in sorted(root.rglob("*.tsv"))
		if path.relative_to(root).parts[0] != "benchmark_summary"]

def benchmark_registry(wildcards):
	registry=[]
	for rule in workflow.rules:
		if rule.benchmark is None or rule.name == "benchmark_summary":
			continue
		threads=rule.resources.get("_cores", 1)
		registry.append(dict(name=rule.name, pattern=str(rule.benchmark),
			threads_declared=threads if isinstance(threads, (int, float)) else None,
			inputs=[str(value) for value in rule.input if not callable(value)]))
	return registry

rule sourmash_clean_reads_report:
	input:
		gather=([SOURMASH_CLEAN_DIR + "/POOLED/pooled_gather_sourmash.with-lineages.csv"] +
			expand(SOURMASH_CLEAN_DIR + "/SAMPLES/{sample}_gather_sourmash.with-lineages.csv", sample=SAMPLES)),
		profiles=([SOURMASH_CLEAN_DIR + "/POOLED/pooled_sourmash.summarized.csv"] +
			expand(SOURMASH_CLEAN_DIR + "/SAMPLES/{sample}_sourmash.summarized.csv", sample=SAMPLES)),
	output:
		html=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_reads.tot.html",
		metadata=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_reads.tot.tsv",
		taxonomy=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_taxonomy.tot.tsv",
		status=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_profile_summary.tot.tsv",
		genus_matrix=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_genus.tot.tsv",
		pooled_png=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_pooled.tot.png",
		pooled_svg=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_pooled.tot.svg",
		genus_png=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_genus.tot.png",
		genus_svg=dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_genus.tot.svg",
	params:
		samples=SAMPLES,
		negative_control=str(config.get("negative_control", "")),
		minimum_shared_bp=int(config.get("sourmash_clean_min_shared_bp", 50000)),
		minimum_reference_fraction=float(config.get("sourmash_clean_min_reference_fraction", 0.10)),
		database=config["sourmash_rocksdb"],
	message:
		"Reporting per-sample and pooled GTDB profiles of all final clean reads"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_clean_reads_report/tot.tsv"
	threads: 1
	resources:
		mem_mb=4000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/01_sourmash_clean_reads.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/01_sourmash_clean_reads.py.ipynb"

rule sourmash_contig_metadata:
	input:
		queries=dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_sourmash_queries.tsv",
		gather=dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_gather_sourmash.csv",
		taxonomy=dirs_dict["ANNOTATION"] + "/sourmash_{sourmash_catalogue}.classifications.csv",
	output:
		metadata=dirs_dict["ANNOTATION"] + "/sourmash_{sourmash_catalogue}.metadata.tsv",
	wildcard_constraints:
		sourmash_catalogue="|".join(SOURMASH_CONTIG_FASTAS),
	params:
		phase="sourmash",
		sampling="tot",
		samples=[],
		results_dir=RESULTS_DIR,
		database=config["sourmash_rocksdb"],
		minimum_length=lambda wildcards: 0 if wildcards.sourmash_catalogue == "combined_microbial_derreplicated_tot" else int(config.get("sourmash_contig_catalogue_min_length", 5000)),
	message:
		"Collecting per-representative Sourmash evidence with the existing catalogue parser"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_contig_metadata/catalogue={sourmash_catalogue}.tsv"
	threads: 1
	resources:
		mem_mb=4000,
	log:
		notebook=dirs_dict["BENCHMARKS"] + "/sourmash_contig_metadata/catalogue={sourmash_catalogue}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_isolate_contig_catalogue.py.ipynb"

rule abundance_lifestyle_summary:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		summary=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_summary.tot.csv",
		rpkm=dirs_dict["MAPPING_DIR"] + "/filtered_RPKM_normalised_tot.txt",
		counts=dirs_dict["MAPPING_DIR"] + "/filtered_counts_normalised_tot.txt",
		coverage=dirs_dict["MAPPING_DIR"] + "/breadth_coverage_percent_tot.txt",
		composition=dirs_dict["ANNOTATION"] + "/nucleotide_content_viral_contigs_clustered_with_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.tsv",
		bacphlip=dirs_dict["ANNOTATION"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot_bacphlip.csv",
		phagcn=dirs_dict["ANNOTATION"] + "/PhaGCN_taxonomy_report_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.csv",
		taxmyphage=dirs_dict["ANNOTATION"] + "/taxmyphage_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot",
		iphop=dirs_dict["ANNOTATION"] + "/iphop_hostID_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot_resultsDir",
		refseq=dirs_dict["ANNOTATION"] + "/blast_output_ViralRefSeq_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.csv",
		metavr=[dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.tot/vOTU_summary.tsv"] if METAVR_blast else [],
		crispr=[dirs_dict["ANNOTATION"] + "/spacepharer_minced_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.tsv"] if config["microbial_spacers"] else [],
		crispr_hosts=[config["abundance_crispr_host_metadata"]] if config.get("abundance_crispr_host_metadata", "") else [],
		sourmash=[dirs_dict["ANNOTATION"] + "/sourmash_filtered_votus_tot.metadata.tsv"] if SOURMASH_CONTIG_CATALOGUE else [],
		sample_metadata=[config["abundance_sample_metadata"]] if config.get("abundance_sample_metadata", "") else [],
	output:
		html=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot.html",
		metadata=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/vOTU_metadata_polished.csv",
		samples=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/sample_summary.csv",
		relative=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/relative_abundance.csv",
		presence=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/presence_absence.csv",
		prevalence=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/vOTU_prevalence.csv",
		categories=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/category_summary.csv",
		genomes=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/genome_summary.csv",
		accumulation=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/sample_accumulation.csv",
		groups=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/group_summary.csv",
		group_presence=dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/group_ubiquitous_vOTUs.csv",
		figures=directory(dirs_dict["PLOTS_DIR"] + "/10_abundance_lifestyle_summary.tot/figures"),
	params:
		samples=SAMPLES,
		group_column=config.get("abundance_group_column", ""),
		common_prevalence=float(config.get("abundance_common_prevalence", 0.2)),
		core_prevalence=float(config.get("abundance_core_prevalence", 0.5)),
	message:
		"Summarizing final vOTU abundance, quality, lifestyle, hosts and prevalence"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/abundance_lifestyle_summary/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/10_abundance_lifestyle_summary.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/10_abundance_lifestyle_summary.py.ipynb"

rule microdiversity_summary:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		summaries=expand(MICRO_DIR + "/{sample}/summary.tsv", sample=SAMPLES),
	output:
		html=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot.html",
		combined=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/vOTU_sample_microdiversity.csv",
		samples=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/sample_summary.csv",
		pi=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/pi_callable.csv",
		callable_fraction=dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/callable_fraction.csv",
		figures=directory(dirs_dict["PLOTS_DIR"] + "/11_microdiversity_summary.tot/figures"),
	params:
		samples=SAMPLES,
	message:
		"Summarizing callable-site nucleotide diversity across samples"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/microdiversity_summary/tot.tsv"
	threads: 1
	resources:
		mem_mb=4000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/11_microdiversity_summary.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/11_microdiversity_summary.py.ipynb"

rule benchmark_summary:
	output:
		jobs=dirs_dict["PLOTS_DIR"] + "/09_benchmark_jobs.csv",
		rules=dirs_dict["PLOTS_DIR"] + "/09_benchmark_rules.csv",
		figure=dirs_dict["PLOTS_DIR"] + "/09_benchmark_summary.png",
		html=dirs_dict["PLOTS_DIR"] + "/09_benchmark_summary.html",
	params:
		benchmark_dir=dirs_dict["BENCHMARKS"],
		registry=benchmark_registry,
		snapshot=benchmark_snapshot,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/benchmark_summary/tot.tsv"
	threads: 1
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/09_benchmark_summary.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/09_benchmark_summary.py.ipynb"

rule plot_assemblies:
	input:
		aa="{fasta}_ORFs.{sampling}.fasta",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
		corrected1=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		corrected2=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		corrected3=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_3_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		corrected4=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_4_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		corrected_medaka=dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.faa",
		scaffolds_pilon1=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_1_{sampling}/pilon.faa"),
		scaffolds_pilon2=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_2_{sampling}/pilon.faa"),
		scaffolds_pilon3=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_3_{sampling}/pilon.faa"),
		scaffolds_pilon4=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_4_{sampling}/pilon.faa"),
	output:
		plot=(dirs_dict["CLEAN_DATA_DIR"] + "/protein_lengths_plot.{sampling}.png"),
		svg=(dirs_dict["CLEAN_DATA_DIR"] + "/protein_lengths_plot.{sampling}.svg"),
	message:
		"Plot unique reads with BBtools"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/plot_assemblies/sampling={sampling}.tsv"
	threads: 1
	run:
		import pandas as pd
		import seaborn as sns; sns.set()
		import matplotlib.pyplot as plt

		plt.figure(figsize=(12,12))
		sns.set(font_scale=2)
		sns.set_style("whitegrid")

		read_max=0

		for h in input.histograms:
			df=pd.read_csv(h, sep="\t")
			df.columns=["count", "percent", "c", "d", "e", "f", "g", "h", "i", "j"]
			df=df[["count", "percent"]]
			ax = sns.lineplot(x="count", y="percent", data=df,err_style='band', label=h.split("/")[-1].split("_kmer")[0])
			read_max=max(read_max,df["count"].max())

		ax.set(ylim=(0, 100))
		ax.set(xlim=(0, read_max*1.2))

		ax.set_xlabel("Read count",fontsize=20)
		ax.set_ylabel("New k-mers (%)",fontsize=20)
		ax.figure.savefig(output.plot)
		ax.figure.savefig(output.svg, format="svg")

rule QC_parsing:
	input:
		inputReadsCount,
		fastp_json=expand(dirs_dict["QC_DIR"] + "/{sample}_fastp.json", sample=SAMPLES),
		fastp_html=expand(dirs_dict["QC_DIR"] + "/{sample}_fastp.html", sample=SAMPLES),
		histograms=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{{sampling}}.csv", sample=SAMPLES),
		preqc_txt=dirs_dict["QC_DIR"]+ "/preQC_illumina_report_data/multiqc_fastqc.txt",
		postqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
		read_count_raw_forward=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + "_read_count.txt", sample=SAMPLES),
		read_count_raw_reverse=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + "_read_count.txt", sample=SAMPLES),
		read_count_trimmed_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_read_count.txt", sample=SAMPLES),
		read_count_trimmed_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_read_count.txt", sample=SAMPLES),
		read_count_duk_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot_read_count.txt", sample=SAMPLES),
		read_count_duk_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot_read_count.txt", sample=SAMPLES),
    	read_count_norm_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.tot_read_count.txt", sample=SAMPLES),
		read_count_norm_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.tot_read_count.txt", sample=SAMPLES),
		supper_dedup=expand(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.log", sample=SAMPLES),
	output:
		kmer_png=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_plot.{sampling}.png"),
		kmer_svg=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_plot.{sampling}.svg"),
		kmer_summary=dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefaction_summary.{sampling}.csv",
		kmer_fit_png=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_plot_fitted.{sampling}.png"),
		kmer_fit_svg=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_plot_fitted.{sampling}.svg"),
		kmer_fit_html=(dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefraction_fitted.{sampling}.html"),
		qc_summary_html=(dirs_dict["PLOTS_DIR"] + "/01_post_qc_read_summary.{sampling}.html"),
		fastp_summary=dirs_dict["PLOTS_DIR"] + "/01_fastp_summary.{sampling}.csv",
		percentage_kept_reads_png=(dirs_dict["PLOTS_DIR"] + "/01_percentage_kept_reads.{sampling}.png"),
		percentage_kept_reads_svg=(dirs_dict["PLOTS_DIR"] + "/01_percentage_kept_reads.{sampling}.svg"),
		percentage_kept_Mbp_png=(dirs_dict["PLOTS_DIR"] + "/01_percentage_kept_Mbp.{sampling}.png"),
		percentage_kept_Mbp_svg=(dirs_dict["PLOTS_DIR"] + "/01_percentage_kept_Mbp.{sampling}.svg"),
		step_qc_reads_html=(dirs_dict["PLOTS_DIR"] + "/01_multistep_qc_report.{sampling}.html"),
		steps_qc_reads_png=(dirs_dict["PLOTS_DIR"] + "/01_qc_bystep_counts.{sampling}.png"),
		steps_qc_reads_svg=(dirs_dict["PLOTS_DIR"] + "/01_qc_bystep_counts.{sampling}.svg"),
		steps_qc_percentage_png=(dirs_dict["PLOTS_DIR"] + "/01_qc_bystep_percentage.{sampling}.png"),
		steps_qc_percentage_svg=(dirs_dict["PLOTS_DIR"] + "/01_qc_bystep_percentage.{sampling}.svg"),
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.{sampling}.csv",
		supperdedupper_html=(dirs_dict["PLOTS_DIR"] + "/01_superdedupper_PCR.{sampling}.html"),
		supperdedupper_png=(dirs_dict["PLOTS_DIR"] + "/01_superdedupper_PCR.{sampling}.png"),
		supperdedupper_svg=(dirs_dict["PLOTS_DIR"] + "/01_superdedupper_PCR.{sampling}.svg"),
	params:
		results_dir=RESULTS_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		samples=SAMPLES,
		forward_tag=config['forward_tag'],
		reverse_tag=config['reverse_tag'],
		raw_dir=dirs_dict["RAW_DATA_DIR"],
		qc_dir=dirs_dict["QC_DIR"],
		remove_euk=REMOVE_EUK,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/QC_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/01_QC.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/01_QC.py.ipynb"


rule assembly_parsing_short:
	input:
		quast_report_dir=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}",
	output:
		log_number_contigs_png=(dirs_dict["PLOTS_DIR"] + "/03_log_number_contigs_plot.{sampling}.png"),
		log_number_contigs_svg=(dirs_dict["PLOTS_DIR"] + "/03_log_number_contigs_plot.{sampling}.svg"),
		contig_length_bp_png=(dirs_dict["PLOTS_DIR"] + "/03_contig_length_bp_plot.{sampling}.png"),
		contig_length_bp_svg=(dirs_dict["PLOTS_DIR"] + "/03_contig_length_bp_plot.{sampling}.svg"),
		contig_number_total_png=(dirs_dict["PLOTS_DIR"] + "/03_contig_number_total_plot.{sampling}.png"),
		contig_number_total_svg=(dirs_dict["PLOTS_DIR"] + "/03_contig_number_total_plot.{sampling}.svg"),
		contig_length_total_png=(dirs_dict["PLOTS_DIR"] + "/03_contig_length_total_plot.{sampling}.png"),
		contig_length_total_svg=(dirs_dict["PLOTS_DIR"] + "/03_contig_length_total_plot.{sampling}.svg"),
	params:
		input_quast_report=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
		results_dir=RESULTS_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		samples=SAMPLES,
		forward_tag=config['forward_tag'],
		reverse_tag=config['reverse_tag'],
		raw_dir=dirs_dict["RAW_DATA_DIR"],
		qc_dir=dirs_dict["QC_DIR"],
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assembly_parsing_short/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_short.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_short.py.ipynb"

rule assembly_parsing_short_RNA:
	input:
		assemblies=expand(RNA_DIR + "/{sample}_{assembler}.fasta", sample=SAMPLES, assembler=RNA_ASSEMBLERS),
		combined=expand(RNA_DIR + "/{sample}_combined_derreplicated.fasta", sample=SAMPLES),
		provenance=expand(RNA_DIR + "/{sample}_assembly_provenance.tsv", sample=SAMPLES),
		quast=RNA_DIR + "/statistics_quast_tot/transposed_report.tsv",
	output:
		summary=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_summary.tot.csv",
		fate=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_provenance.tot.csv",
		contig_counts_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_contig_counts.tot.png",
		contig_counts_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_contig_counts.tot.svg",
		total_length_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_total_length.tot.png",
		total_length_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_total_length.tot.svg",
		length_distribution_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_length_distribution.tot.png",
		length_distribution_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_length_distribution.tot.svg",
		provenance_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_provenance.tot.png",
		provenance_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_short_RNA_provenance.tot.svg",
	params:
		samples=SAMPLES,
		assemblers=RNA_ASSEMBLERS,
		min_length=int(config.get("rna_min_contig_length", 500)),
	message:
		"Summarizing RNA assemblies across samples and assemblers"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assembly_parsing_short_RNA/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_short_RNA.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_short_RNA.py.ipynb"

rule assembly_parsing_long:
	input:
		caudovirales=("db/caudovirales_orf_lengths_09_05_2023.txt"),
		hybrid=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds_ORFs_length.tot.txt"),
		canu=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER +"_ORFs_length.tot.txt"),
		medaka=(dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		racon1=(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_1_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		racon2=(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		scaffolds_pilon1_final=(dirs_dict["ASSEMBLY_DIR"] + "/pilon_1_polished_{sample}_contigs_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		scaffolds_pilon2_final=(dirs_dict["ASSEMBLY_DIR"] + "/pilon_2_polished_{sample}_contigs_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		scaffolds_pilon3_final=(dirs_dict["ASSEMBLY_DIR"] + "/pilon_3_polished_{sample}_contigs_"+ LONG_ASSEMBLER + "_ORFs_length.tot.txt"),
		scaffolds_pilon4_final=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon_ORFs_length.tot.txt"),
	output:
		orf_length_png=(dirs_dict["PLOTS_DIR"] + "/03_ORF_length_{sample}.png"),
		orf_length_svg=(dirs_dict["PLOTS_DIR"] + "/03_ORF_length_{sample}.svg"),
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assembly_parsing_long/sample={sample}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_long_{sample}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_long.py.ipynb"

def inputAssemblyContigs(wildcards):
	inputs=[]
	inputs.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES))
	if NANOPORE:
		inputs.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample_nanopore}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{{sampling}}.fasta", sample_nanopore=NANOPORE_SAMPLES))
	if CROSS_ASSEMBLY:
		inputs.append(dirs_dict["ASSEMBLY_DIR"] + "/ALL_spades_filtered_scaffolds.{sampling}.fasta")
	if SUBASSEMBLY:
		inputs.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_metaspades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES, subsample=subsample_test)),
	return inputs

rule viralID_parsing:
	input:
		viral_sequences=input_vOTU_clustering,
		assembled_sequences=inputAssemblyContigs,
	output:
		viral_sequences_count_plot_png=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_count_{sampling}.png"),
		viral_sequences_count_plot_svg=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_count_{sampling}.svg"),
		viral_sequences_count_table_html=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_count_{sampling}.html"),
		viral_sequences_length_plot_png=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_length_{sampling}.png"),
		viral_sequences_length_plot_svg=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_length_{sampling}.svg"),
		viral_sequences_length_table_html=(dirs_dict["PLOTS_DIR"] + "/04_viral_sequences_length_{sampling}.html"),
	params:
		samples=SAMPLES,
		contig_dir=dirs_dict["ASSEMBLY_DIR"],
		viral_dir=dirs_dict['VIRAL_DIR'],
		subassembly=SUBASSEMBLY,
		cross_assembly=CROSS_ASSEMBLY,
		long_assembler=LONG_ASSEMBLER
	benchmark:
		dirs_dict["BENCHMARKS"] + "/viralID_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/04_viral_ID_{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/04_viral_ID.py.ipynb"

def input_assembly_flagstats(wildcards):
	if ISOLATES:
		return []
	inputs=[]
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT)),
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_assembled_contigs.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT)),
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_viral_contigs.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT)),
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_unfiltered_contigs.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT)),
	return inputs

rule mapping_statistics_parsing:
	input:
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.{sampling}.csv",
		assembled_sequences=inputAssemblyContigs,
		assembly_flagstats=input_assembly_flagstats,
		accounting=expand(ISOLATE_MAPPING_DIR + "/{sample}/{stage}.summary.tsv", sample=SAMPLES, stage=ISOLATE_STAGES) if ISOLATES else [],
		all_assembled_flagstats=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_flagstats_filtered_all_assembled_{sample}.tot.txt", sample=SAMPLES) if ISOLATES and MAP_TO_ALL_ASSEMBLED else [],
		all_assembled_mapped_pairs=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_mapped_pairs_filtered_all_assembled_{sample}.tot.txt", sample=SAMPLES) if MAP_TO_ALL_ASSEMBLED and not ISOLATES else [],
		all_assembled_top_summary=[ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_mapping_summary_tot.tsv"] if MAP_TO_ALL_ASSEMBLED else [],
		all_assembled_top_rpkm_png=[dirs_dict["PLOTS_DIR"] + "/07_all_assembled_top_RPKM_clustermap_tot.png"] if MAP_TO_ALL_ASSEMBLED else [],
	output:
		mapping_stats_html=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_{sampling}.html"),
		filtered_viral_png=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_filtered_viral_{sampling}.png") if not ISOLATES else [],
		filtered_viral_svg=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_filtered_viral_{sampling}.svg") if not ISOLATES else [],
		filtered_unfiltered_png=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_filtered_unfiltered_{sampling}.png") if not ISOLATES else [],
		filtered_unfiltered_svg=(dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_filtered_unfiltered_{sampling}.svg") if not ISOLATES else [],
		all_assembled_png=[dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_all_assembled_{sampling}.png"] if MAP_TO_ALL_ASSEMBLED and not ISOLATES else [],
		all_assembled_svg=[dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_all_assembled_{sampling}.svg"] if MAP_TO_ALL_ASSEMBLED and not ISOLATES else [],
		isolate_png=[dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_isolates_{sampling}.png"] if ISOLATES else [],
		isolate_svg=[dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_isolates_{sampling}.svg"] if ISOLATES else [],
		isolate_table=[dirs_dict["PLOTS_DIR"] + "/07_mapping_statistics_isolates_{sampling}.tsv"] if ISOLATES else [],
	params:
		isolates=ISOLATES,
		stages=ISOLATE_STAGES,
		samples=SAMPLES,
		mapping_dir=dirs_dict["MAPPING_DIR"],
		sampling="{sampling}",
		map_to_all_assembled=MAP_TO_ALL_ASSEMBLED,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapping_statistics_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_mapping_statistics_{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + ("/07_mapping_statistics_isolates.py.ipynb" if ISOLATES else "/07_mapping_statistics.py.ipynb")


def input_phage_isolates_assembled_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs.{sampling}_covstats.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_viral_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs.{sampling}_covstats.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_unfiltered_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs.{sampling}_covstats.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_filtered_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_covstats.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_filtered_flagstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.{sampling}.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_assembled_flagstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_assembled_contigs.{sampling}.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_viral_flagstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_viral_contigs.{sampling}.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_unfiltered_flagstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_unfiltered_contigs.{sampling}.txt", sample=SAMPLES, sampling=wildcards.sampling)


def input_phage_isolates_host_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_filtered_{sample}_vs_{host}_covstats.txt", host=HOSTS, sample=SAMPLES)


def input_phage_isolates_host_masked_covstats(wildcards):
	return expand(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_filtered_{sample}_vs_{host}_masked_prophages_covstats.txt", host=HOSTS, sample=SAMPLES)


def input_phage_isolates_host_blast(wildcards):
	return expand(dirs_dict["vOUT_DIR"] + "/blastn_out_assembly_{host}.{sampling}.csv", host=HOSTS, sampling=wildcards.sampling)


def input_phage_isolates_host_fastas(wildcards):
	return expand(dirs_dict["HOST_DIR"] + "/{host}.fasta", host=HOSTS)


def input_isolate_host_assignments(wildcards):
	return HOST_ASSIGNMENT_FILES


rule phage_isolates_catalogue:
	input:
		host_assignments=input_isolate_host_assignments,
		selection=expand(ISOLATE_CONTIG_DIR + "/SELECTION/{sample}.tot.tsv", sample=SAMPLES),
		catalogue_provenance=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_provenance.tot.tsv",
		catalogue_exact_clusters=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_cluster.tot.tsv",
		catalogue_mosaic_clusters=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_rep_seq.tot_95-85.clstr",
		catalogue_genomad=ALL_ASSEMBLED_MAPPING_DIR + "/phage_isolates_contigs_existing_annotations_tot.tsv",
		catalogue_checkv=expand(dirs_dict["ASSEMBLY_DIR"] + "/checkV_isolates_{sample}_tot/quality_summary.tsv", sample=SAMPLES),
		catalogue_contamination=expand(dirs_dict["ASSEMBLY_DIR"] + "/checkV_isolates_{sample}_tot/contamination.tsv", sample=SAMPLES),
		host_blast=input_phage_isolates_host_blast,
		host_viral_blast=expand(dirs_dict["vOUT_DIR"] + "/blastn_out_assembly_{host}_viral_regions.tot.csv", host=HOSTS),
		catalogue_sourmash_tax=[dirs_dict["ANNOTATION"] + "/sourmash_phage_isolates_cluster_representatives_tot.classifications.csv"] if SOURMASH_CONTIG_CATALOGUE else [],
		catalogue_sourmash_gather=[dirs_dict["ANNOTATION"] + "/phage_isolates_cluster_representatives_tot_gather_sourmash.csv"] if SOURMASH_CONTIG_CATALOGUE else [],
		catalogue_sourmash_queries=[dirs_dict["ANNOTATION"] + "/phage_isolates_cluster_representatives_tot_sourmash_queries.tsv"] if SOURMASH_CONTIG_CATALOGUE else [],
	output:
		all_contig_metadata=temp(ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/all_contig_metadata.unflagged.tsv"),
		cluster_metadata=temp(ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/cluster_metadata.unflagged.tsv"),
		sample_hosts=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/sample_host_assignments.tsv",
	params:
		phase="catalogue",
		samples=SAMPLES,
		hosts=HOSTS,
		sampling="{sampling}",
		results_dir=RESULTS_DIR,
		min_depth=float(config.get("isolate_min_depth", 5)),
		min_length=int(config.get("isolate_min_length_bp", 4000)),
		short_depth=float(config.get("isolate_short_min_depth", 10)),
		host_min_identity=float(config.get("isolate_host_min_identity", 90)),
		host_min_query_coverage=float(config.get("isolate_host_min_query_coverage", 90)),
		sourmash_contig_catalogue_database=config["sourmash_rocksdb"],
	benchmark:
		dirs_dict["BENCHMARKS"] + "/phage_isolates_catalogue/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/08_isolate_contig_catalogue.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_isolate_contig_catalogue.py.ipynb"


rule isolate_fragment_candidates:
	input:
		metadata=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/all_contig_metadata.unflagged.tsv",
		clusters=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/cluster_metadata.unflagged.tsv",
		ani=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_rep_seq.tot-aniout.csv",
		fasta=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs.tot.fasta",
	output:
		metadata=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/all_contig_metadata.tsv",
		clusters=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/cluster_metadata.tsv",
		matches=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/fragment_candidate_matches.tsv",
	params:
		min_identity=float(config.get("isolate_fragment_min_identity", 95)),
		min_short_coverage=float(config.get("isolate_fragment_min_short_coverage", 95)),
		max_length_ratio=float(config.get("isolate_fragment_max_length_ratio", 0.80)),
	message:
		"Flagging retained fragment candidates without changing clustering or selection"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/isolate_fragment_candidates/sampling={sampling}.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	shell:
		r"""
		python - {input.metadata:q} {input.clusters:q} {input.ani:q} {input.fasta:q} \
			{output.metadata:q} {output.clusters:q} {output.matches:q} \
			{params.min_identity} {params.min_short_coverage} {params.max_length_ratio} <<-'PYTHON'
		import sys
		import numpy as np
		import pandas as pd
		from Bio import SeqIO

		def flag_fragment_candidates(metadata, ani, sequences, min_identity, min_short_coverage, max_length_ratio):
		    from Bio.Seq import reverse_complement

		    # Only retained originals are compared; cluster membership is never changed.
		    retained=metadata.loc[metadata["retained"], ["sample", "original_id", "exact_rep", "cluster_rep", "length_bp"]].copy()
		    targets=retained.rename(columns={{"sample": "contained_in_sample", "original_id": "contained_in_contig",
		        "exact_rep": "contained_in_exact_rep", "cluster_rep": "contained_in_cluster",
		        "length_bp": "contained_in_length_bp"}})
		    columns=list(retained.columns) + list(targets.columns) + [
		        "containment_identity_percent", "containment_short_coverage_percent", "containment_length_ratio",
		        "containment_evidence_scope"]

		    # ANI percentages apply to exact representatives, not to shorter MMseqs members.
		    # Expand them only to verified full-length identical copies (either strand).
		    same_sequence=[]
		    for row in retained.itertuples():
		        sequence=sequences.get(row.original_id, "")
		        representative=sequences.get(row.exact_rep, "")
		        same_sequence.append(bool(sequence) and (sequence == representative or
		                             reverse_complement(sequence) == representative))
		    aliases=retained.loc[same_sequence]
		    alias_targets=targets[targets["contained_in_contig"].isin(aliases["original_id"])]
		    ani=ani.reindex(columns=["qname", "tname", "pid", "qcov"]).astype({{"qname": object, "tname": object}})
		    ani=ani[ani["pid"].ge(min_identity) & ani["qcov"].ge(min_short_coverage)]
		    matches=ani.merge(aliases, left_on="qname", right_on="exact_rep").merge(
		        alias_targets, left_on="tname", right_on="contained_in_exact_rep")
		    matches["containment_length_ratio"]=matches["length_bp"] / matches["contained_in_length_bp"].replace(0, np.nan)
		    matches=matches[matches["length_bp"].lt(matches["contained_in_length_bp"]) &
		                    matches["containment_length_ratio"].le(max_length_ratio)].copy()
		    matches=matches.rename(columns={{"pid": "containment_identity_percent", "qcov": "containment_short_coverage_percent"}})
		    matches["containment_evidence_scope"]="ANI alignment (exact representatives)"

		    # Exact dereplication can also contain shorter sequences. Verify their actual
		    # sequences within each existing MMseqs group instead of inheriting ANI values.
		    exact_matches=[]
		    if min_identity <= 100 and min_short_coverage <= 100:
		        for exact_rep, members in retained.groupby("exact_rep", sort=False):
		            for query in members.itertuples():
		                sequence=sequences.get(query.original_id, "")
		                if not sequence:
		                    continue
		                reverse=reverse_complement(sequence)
		                longer=members[members["length_bp"].gt(query.length_bp) &
		                               (query.length_bp / members["length_bp"]).le(max_length_ratio)]
		                for target in longer.itertuples():
		                    parent=sequences.get(target.original_id, "")
		                    if sequence in parent or reverse in parent:
		                        exact_matches.append(dict(sample=query.sample, original_id=query.original_id,
		                            exact_rep=exact_rep, cluster_rep=query.cluster_rep, length_bp=query.length_bp,
		                            contained_in_sample=target.sample, contained_in_contig=target.original_id,
		                            contained_in_exact_rep=exact_rep, contained_in_cluster=target.cluster_rep,
		                            contained_in_length_bp=target.length_bp, containment_identity_percent=100.0,
		                            containment_short_coverage_percent=100.0,
		                            containment_length_ratio=query.length_bp / target.length_bp,
		                            containment_evidence_scope="Exact sequence containment (MMseqs members)"))
		    # Prefer the strongest identity, then coverage, longest target and lexical ID.
		    tables=[table for table in [pd.DataFrame(exact_matches, columns=columns), matches.reindex(columns=columns)] if not table.empty]
		    matches=pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=columns)
		    matches=matches.sort_values(["sample", "original_id", "containment_identity_percent",
		        "containment_short_coverage_percent", "contained_in_length_bp", "contained_in_contig"],
		        ascending=[True, True, False, False, False, True], kind="stable").drop_duplicates(
		        ["original_id", "contained_in_contig"])
		    matches["match_rank"]=matches.groupby("original_id", sort=False).cumcount() + 1
		    best=matches.drop_duplicates("original_id").set_index("original_id")
		    metadata=metadata.copy()
		    metadata["fragment_candidate"]=metadata["original_id"].isin(best.index)
		    for column in list(targets.columns) + columns[-4:]:
		        metadata[column]=metadata["original_id"].map(best[column])
		    return metadata, matches

		(metadata_path, cluster_path, ani_path, fasta_path, output_metadata, output_clusters,
		 output_matches, min_identity, min_short_coverage, max_length_ratio)=sys.argv[1:]
		metadata=pd.read_csv(metadata_path, sep="\t")
		clusters=pd.read_csv(cluster_path, sep="\t")
		ani=pd.read_csv(ani_path, sep="\t")
		with open(fasta_path) as handle:
		    sequences={{record.id: str(record.seq).upper() for record in SeqIO.parse(handle, "fasta")}}
		metadata, matches=flag_fragment_candidates(metadata, ani, sequences, float(min_identity),
		    float(min_short_coverage), float(max_length_ratio))
		counts=metadata.groupby("cluster_rep")["fragment_candidate"].sum()
		clusters["number_fragment_candidates"]=clusters["cluster_rep"].map(counts).fillna(0).astype(int)
		for table, path in [(metadata, output_metadata), (clusters, output_clusters), (matches, output_matches)]:
		    table.to_csv(path, sep="\t", index=False, na_rep="not reported")

		PYTHON
		"""


rule phage_isolates_summary:
	input:
		viridic_selection=lambda wc: input_isolate_viridic_selection(wc),
		viridic_members=lambda wc: input_isolate_viridic_members(wc),
		viridic_results=lambda wc: input_isolate_viridic_results(wc),
		host_assignments=input_isolate_host_assignments,
		catalogue_provenance=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_provenance.tot.tsv",
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.{sampling}.csv",
		kmer_summary=dirs_dict["PLOTS_DIR"] + "/01_kmer_rarefaction_summary.{sampling}.csv",
		pcr_duplicates=expand(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.log", sample=SAMPLES),
		kraken_reports=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_tot.csv", sample=SAMPLES),
		quast=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
		metadata=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/all_contig_metadata.tsv",
		clusters=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/cluster_metadata.tsv",
		fragment_matches=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/fragment_candidate_matches.tsv",
		catalogue_blast=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_rep_seq.tot-blastout.csv",
		accounting=expand(ISOLATE_MAPPING_DIR + "/{sample}/{stage}.summary.tsv", sample=SAMPLES, stage=ISOLATE_STAGES),
		reference_counts=expand(ISOLATE_MAPPING_DIR + "/{sample}/{stage}.reference_reads.tsv", sample=SAMPLES, stage=ISOLATE_STAGES),
		unexplained=expand(ISOLATE_MAPPING_DIR + "/{sample}/unexplained_{mate}.fastq.gz", sample=SAMPLES, mate=["R1", "R2", "unpaired"]),
		host_activity=[dirs_dict["PLOTS_DIR"] + "/08_host_prophage_activity.tot.tsv"] if HOST_PROPHAGE_ACTIVITY and HOSTS else [],
		host_covstats=input_phage_isolates_host_covstats if HOST_IDENTIFICATION_TEST else [],
		host_masked_covstats=input_phage_isolates_host_masked_covstats if HOST_IDENTIFICATION_TEST else [],
		host_blast=input_phage_isolates_host_blast,
		viral_refseq_blast=dirs_dict["ANNOTATION"] + "/blast_output_ViralRefSeq_phage_isolates_cluster_representatives.tot.csv",
		metavr_blast=[dirs_dict["ANNOTATION"] + "/blast_output_METAVR_phage_isolates_cluster_representatives.tot.csv"] if METAVR_blast else [],
		clean_sourmash_report_metadata=[dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_reads.tot.tsv"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_report_taxonomy=[dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_taxonomy.tot.tsv"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_report_status=[dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_profile_summary.tot.tsv"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_report_pooled_png=[dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_pooled.tot.png"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_report_pooled_svg=[dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_pooled.tot.svg"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_report_genus_png=[dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_genus.tot.png"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_report_genus_svg=[dirs_dict["PLOTS_DIR"] + "/01_sourmash_clean_genus.tot.svg"] if SOURMASH_CLEAN_READS else [],
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_summary.{sampling}.html",
		summary_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_summary.{sampling}.csv",
		contig_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_contigs.{sampling}.csv",
		closest_relatives_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_closest_relatives.{sampling}.csv",
		host_mapping_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_mapping.{sampling}.csv",
		cluster_members_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_members.{sampling}.csv",
		cluster_summary_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_summary.{sampling}.csv",
		cluster_samples=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_sample_clusters.{sampling}.tsv",
		cluster_coverage_matrix_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_coverage_matrix.{sampling}.csv",
		cluster_counts_matrix_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_counts_matrix.{sampling}.csv",
		cluster_annotations=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_annotations.{sampling}.tsv",
		cluster_coverage_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_coverage.{sampling}.png",
		cluster_coverage_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_cluster_coverage.{sampling}.svg",
		viridic_figures=[directory(dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_viridic.{sampling}")] if ISOLATE_VIRIDIC else [],
		single_contig_samples_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_single_contig_samples.{sampling}.csv",
		host_blast_summary_csv=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_blast_summary.{sampling}.csv",
		read_accounting=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_read_accounting.{sampling}.tsv",
		purity_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_read_accounting.{sampling}.png",
		purity_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_read_accounting.{sampling}.svg",
		mapping_heatmap_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_mapping_heatmap.{sampling}.png",
		mapping_heatmap_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_mapping_heatmap.{sampling}.svg",
		remaining_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_remaining_contigs.{sampling}.png",
		remaining_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_remaining_contigs.{sampling}.svg",
		contig_selection_counts=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_contig_selection.{sampling}.tsv",
		contig_selection_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_contig_selection.{sampling}.png",
		contig_selection_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_contig_selection.{sampling}.svg",
		host_blast_png=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_blast_clustermap.{sampling}.png",
		host_blast_svg=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_host_blast_clustermap.{sampling}.svg",
		clean_sourmash_metadata=[dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_sourmash_clean_reads.{sampling}.tsv"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_taxonomy=[dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_sourmash_clean_taxonomy.{sampling}.tsv"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_pooled_png=[dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_sourmash_clean_pooled.{sampling}.png"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_pooled_svg=[dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_sourmash_clean_pooled.{sampling}.svg"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_genus_png=[dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_sourmash_clean_genus.{sampling}.png"] if SOURMASH_CLEAN_READS else [],
		clean_sourmash_genus_svg=[dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_sourmash_clean_genus.{sampling}.svg"] if SOURMASH_CLEAN_READS else [],
	params:
		phase="summary",
		ani_script=workflow.basedir + "/scripts/anicalc_checkv.py",
		hosts=HOSTS,
		host_min_identity=float(config.get("isolate_host_min_identity", 90)),
		host_min_query_coverage=float(config.get("isolate_host_min_query_coverage", 90)),
		host_identification_test=HOST_IDENTIFICATION_TEST,
		min_mapping_percent=float(config.get("isolate_min_retained_mapping_percent", 70)),
		max_host_percent=float(config.get("isolate_max_host_mapping_percent", 10)),
		max_unexplained_percent=float(config.get("isolate_max_unexplained_percent", 10)),
		host_test_min_reads=int(config.get("isolate_host_test_min_reads", 100)),
		host_test_min_breadth=float(config.get("isolate_host_test_min_breadth_percent", 1)),
		plot_max_host_contigs=int(config.get("isolate_plot_max_host_contigs", 100)),
		stages=ISOLATE_STAGES,
		samples=SAMPLES,
		sampling="{sampling}",
		results_dir=RESULTS_DIR,
		isolates=ISOLATES,
		metagenome=METAGENOME,
		microbial=MICROBIAL,
		remove_euk=REMOVE_EUK,
		sourmash_contig_catalogue=SOURMASH_CONTIG_CATALOGUE,
		sourmash_contig_catalogue_min_length=int(config.get("sourmash_contig_catalogue_min_length", 5000)),
		sourmash_contig_catalogue_database=config["sourmash_rocksdb"],
		sourmash_clean_reads=SOURMASH_CLEAN_READS,
		sourmash_clean_min_shared_bp=int(config.get("sourmash_clean_min_shared_bp", 50000)),
		sourmash_clean_min_reference_fraction=float(config.get("sourmash_clean_min_reference_fraction", 0.10)),
		sourmash_clean_database=config["sourmash_rocksdb"],
	benchmark:
		dirs_dict["BENCHMARKS"] + "/phage_isolates_summary/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/08_phage_isolates_summary.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_phage_isolates_summary.py.ipynb"


rule subsample_reads:
	input:
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.tot.csv",
		flagstats=expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.tot.txt", sample=SAMPLES),
	output:
		viral_subsampling=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sub_sampling_reads.txt", sample=SAMPLES),
	params:
		samples=SAMPLES,
		mapping_dir=dirs_dict["MAPPING_DIR"],
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="tot",
		key_samples=SAMPLES_key
	benchmark:
		dirs_dict["BENCHMARKS"] + "/subsample_reads/tot.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_subsampling.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_subsampling.py.ipynb"

rule normalise_reads:
	input:
		postqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
		covstats=expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{{sampling}}_covstats.txt", sample=SAMPLES),
		covstats_unique=expand(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{{sampling}}_unique_covstats.txt", sample=SAMPLES),
	output:
		raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/RPKM_raw_{sampling}.txt",
		norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/RPKM_normalised_{sampling}.txt",
		raw_count_file=dirs_dict["MAPPING_DIR"] + "/counts_raw_{sampling}.txt",
		norm_count_file=dirs_dict["MAPPING_DIR"] + "/counts_normalised_{sampling}.txt",
		coverage_RPKM_file=dirs_dict["MAPPING_DIR"] + "/breadth_coverage_percent_{sampling}.txt",
		coverage_bases_RPKM_file=dirs_dict["MAPPING_DIR"] + "/breadth_coverage_bases_{sampling}.txt",
		filtered_raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/filtered_RPKM_raw_{sampling}.txt",
		filtered_norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/filtered_RPKM_normalised_{sampling}.txt",
		filtered_raw_count_file=dirs_dict["MAPPING_DIR"] + "/filtered_counts_raw_{sampling}.txt",
		filtered_norm_count_file=dirs_dict["MAPPING_DIR"] + "/filtered_counts_normalised_{sampling}.txt",
		filtered_75_raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/filtered_75_RPKM_raw_{sampling}.txt",
		filtered_75_norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/filtered_75_RPKM_normalised_{sampling}.txt",
	params:
		samples=SAMPLES,
		mapping_dir=dirs_dict["MAPPING_DIR"],
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="{sampling}",
		threshold_bases=200,
		threshold_RPKM=0.1,
		reference="",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalise_reads/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"

rule normalise_reads_reference:
	input:
		postqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
		covstats=expand(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_{{sampling}}_covstats.txt", sample=SAMPLES),
		covstats_unique=expand(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_{{sampling}}_unique_covstats.txt", sample=SAMPLES),
	output:
		raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_RPKM_raw_{sampling}.txt",
		norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_RPKM_normalised_{sampling}.txt",
		raw_count_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_counts_raw_{sampling}.txt",
		norm_count_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_counts_normalised_{sampling}.txt",
		coverage_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_breadth_coverage_percent_{sampling}.txt",
		coverage_bases_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/" + REFERENCE + "_breadth_coverage_bases_{sampling}.txt",
		filtered_raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_" + REFERENCE + "_RPKM_raw_{sampling}.txt",
		filtered_norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_" + REFERENCE + "_RPKM_normalised_{sampling}.txt",
		filtered_raw_count_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_" + REFERENCE + "_counts_raw_{sampling}.txt",
		filtered_norm_count_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_" + REFERENCE + "_counts_normalised_{sampling}.txt",
		filtered_75_raw_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_75_" + REFERENCE + "_RPKM_raw_{sampling}.txt",
		filtered_75_norm_RPKM_file=dirs_dict["MAPPING_DIR"] + "/REFERENCES/filtered_75_" + REFERENCE + "_RPKM_normalised_{sampling}.txt",
	params:
		samples=SAMPLES,
		mapping_dir=dirs_dict["MAPPING_DIR"],
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="{sampling}",
		threshold_bases=200,
		threshold_RPKM=0.1,
		reference=REFERENCE,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalise_reads_reference/" + REFERENCE + "/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise_" + REFERENCE + ".{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"


# Reuse the reference-normalisation notebook and criteria, with separate inputs/outputs.
rule normalise_reads_RefSeq:
	input:
		postqc_txt=dirs_dict["QC_DIR"] + "/postQC_illumina_report_data/multiqc_fastqc.txt",
		covstats=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_covstats.txt", sample=SAMPLES),
		covstats_unique=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_covstats.txt", sample=SAMPLES),
	output:
		raw_RPKM_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_RPKM_raw_tot.txt",
		norm_RPKM_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_RPKM_normalised_tot.txt",
		raw_count_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_counts_raw_tot.txt",
		norm_count_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_counts_normalised_tot.txt",
		coverage_RPKM_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_breadth_coverage_percent_tot.txt",
		coverage_bases_RPKM_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_breadth_coverage_bases_tot.txt",
		mean_coverage_file=REFSEQ_MAPPING_DIR + "/RefSeqViral_mean_depth_tot.txt",
		filtered_raw_RPKM_file=REFSEQ_MAPPING_DIR + "/filtered_RefSeqViral_RPKM_raw_tot.txt",
		filtered_norm_RPKM_file=REFSEQ_MAPPING_DIR + "/filtered_RefSeqViral_RPKM_normalised_tot.txt",
		filtered_raw_count_file=REFSEQ_MAPPING_DIR + "/filtered_RefSeqViral_counts_raw_tot.txt",
		filtered_norm_count_file=REFSEQ_MAPPING_DIR + "/filtered_RefSeqViral_counts_normalised_tot.txt",
		filtered_75_raw_RPKM_file=REFSEQ_MAPPING_DIR + "/filtered_75_RefSeqViral_RPKM_raw_tot.txt",
		filtered_75_norm_RPKM_file=REFSEQ_MAPPING_DIR + "/filtered_75_RefSeqViral_RPKM_normalised_tot.txt",
	params:
		samples=SAMPLES,
		mapping_dir=REFSEQ_MAPPING_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="tot",
		threshold_bases=200,
		threshold_RPKM=0.1,
		reference="RefSeqViral",
		index_label="RefSeq_accession",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalise_reads_RefSeq/tot.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise_RefSeqViral.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"

rule normalise_reads_all_assembled:
	input:
		postqc_txt=dirs_dict["QC_DIR"] + "/postQC_illumina_report_data/multiqc_fastqc.txt",
		covstats=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_all_assembled_{sample}_tot_covstats.txt", sample=SAMPLES),
		covstats_unique=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_all_assembled_{sample}_tot_unique_covstats.txt", sample=SAMPLES),
	output:
		raw_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_RPKM_raw_tot.txt",
		norm_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_RPKM_normalised_tot.txt",
		raw_count_file=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_counts_raw_tot.txt",
		norm_count_file=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_counts_normalised_tot.txt",
		coverage_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_breadth_coverage_percent_tot.txt",
		coverage_bases_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_breadth_coverage_bases_tot.txt",
		mean_coverage_file=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_mean_depth_tot.txt",
		filtered_raw_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_all_assembled_RPKM_raw_tot.txt",
		filtered_norm_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_all_assembled_RPKM_normalised_tot.txt",
		filtered_raw_count_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_all_assembled_counts_raw_tot.txt",
		filtered_norm_count_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_all_assembled_counts_normalised_tot.txt",
		filtered_75_raw_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_75_all_assembled_RPKM_raw_tot.txt",
		filtered_75_norm_RPKM_file=ALL_ASSEMBLED_MAPPING_DIR + "/filtered_75_all_assembled_RPKM_normalised_tot.txt",
	params:
		samples=SAMPLES,
		mapping_dir=ALL_ASSEMBLED_MAPPING_DIR,
		clean_dir=dirs_dict["CLEAN_DATA_DIR"],
		sampling="tot",
		threshold_bases=200,
		threshold_RPKM=0.1,
		reference="all_assembled",
		index_label="assembled_contig",
		plot_max_points=5000,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalise_reads_all_assembled/tot.tsv"
	resources:
		mem_mb=64000
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_Normalise_all_assembled.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_Normalise.py.ipynb"

rule all_assembled_mapping_summary:
	input:
		qc=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.tot.csv",
		mapped_pairs=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_mapped_pairs_filtered_{{assembled_reference}}_{sample}.tot.txt", sample=SAMPLES),
		full_flagstats=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_flagstats_filtered_{{assembled_reference}}_{sample}.tot.txt", sample=SAMPLES) if ISOLATES else [],
		accounting=expand(ISOLATE_MAPPING_DIR + "/{sample}/01_own_retained.summary.tsv", sample=SAMPLES) if ISOLATES else [],
	output:
		tsv=ALL_ASSEMBLED_MAPPING_DIR + "/{assembled_reference}_mapping_summary_tot.tsv",
	params:
		samples=SAMPLES,
		isolates=ISOLATES,
	wildcard_constraints:
		assembled_reference="all_assembled|all_assembled_top",
	message:
		"Summarizing reads mapping to {wildcards.assembled_reference}"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/all_assembled_mapping_summary/reference={assembled_reference}.tsv"
	threads: 1
	run:
		import csv

		with open(output.tsv, "w") as handle:
			writer=csv.writer(handle, delimiter="\t", lineterminator="\n")
			if getattr(params, "isolates", False):
				writer.writerow(["sample", "mapping_input_reads", "mapped_reads", "unmapped_reads", "mapped_percent", "read_scope"])
				for sample, path, ledger in zip(params.samples, input.full_flagstats, input.accounting):
					with open(ledger) as counts:
						total=int(next(csv.DictReader(counts, delimiter="\t"))["original_reads"])
					with open(path) as counts:
						mapped=next((int(line.split()[0]) for line in counts if " primary mapped (" in line), 0)
					writer.writerow([sample, total, mapped, total - mapped,
						round(100 * mapped / total, 2) if total else 0, "all_qc_passed_paired_and_orphan_reads"])
			else:
				with open(input.qc) as counts:
					cleaned_reads={row["sample"]: int(float(row["bbduk"])) for row in csv.DictReader(counts)}
				writer.writerow(["sample", "cleaned_read_pairs", "subsampled_read_pairs", "properly_mapped_pairs", "properly_mapped_percent"])
				for sample, path in zip(params.samples, input.mapped_pairs):
					with open(path) as counts:
						mapped=int(counts.read().strip())
					cleaned=cleaned_reads[sample]
					subsampled=min(2000000, cleaned)
					writer.writerow([sample, cleaned, subsampled, mapped, round(100 * mapped / subsampled, 2) if subsampled else 0])


rule all_assembled_top_abundance_plot:
	input:
		covstats=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_all_assembled_top_{sample}_tot_covstats.txt", sample=SAMPLES),
	output:
		rpkm=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_RPKM_raw_tot.txt",
		counts=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_counts_raw_tot.txt",
		breadth=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_breadth_coverage_percent_tot.txt",
		covered_bases=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_breadth_coverage_bases_tot.txt",
		depth=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_mean_depth_tot.txt",
		png=dirs_dict["PLOTS_DIR"] + "/07_all_assembled_top_RPKM_clustermap_tot.png",
		svg=dirs_dict["PLOTS_DIR"] + "/07_all_assembled_top_RPKM_clustermap_tot.svg",
	params:
		samples=SAMPLES,
		isolates=ISOLATES,
	message:
		"Plotting remapped raw RPKM for the selected cluster representatives"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/all_assembled_top_abundance_plot/tot.tsv"
	threads: 1
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_all_assembled_top_abundance.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/07_all_assembled_top_abundance.py.ipynb"

rule select_all_assembled_top_contigs:
	input:
		fasta=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_derreplicated_rep_seq.tot.fasta",
		clusters=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_derreplicated_cluster.tot.tsv",
		provenance=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_provenance.tot.tsv",
		rpkm=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_RPKM_raw_tot.txt",
		counts=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_counts_raw_tot.txt",
		breadth=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_breadth_coverage_percent_tot.txt",
		depth=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_mean_depth_tot.txt",
		unique=expand(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_all_assembled_{sample}_tot_unique_covstats.txt", sample=SAMPLES),
	output:
		fasta=ALL_ASSEMBLED_SELECTED_PREFIX + ".fasta",
		ranking=ALL_ASSEMBLED_DIR + "/all_assembled_selected_contigs_tot.tsv",
		membership=ALL_ASSEMBLED_DIR + "/all_assembled_selected_contigs_membership_tot.tsv",
	params:
		samples=SAMPLES,
		top_n=int(config.get("all_assembled_cluster_top_n", 1000)),
		top_per_sample=int(config.get("all_assembled_top_per_sample", 100)),
		negative_control=str(config.get("negative_control", "")).strip(),
	message:
		"Selecting the union of sample top contigs and mean raw-RPKM top contigs for clustering"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/select_all_assembled_top_contigs/tot.tsv"
	threads: 1
	resources:
		mem_mb=16000,
	run:
		import pandas as pd

		rpkm=pd.read_csv(input.rpkm, index_col=0).reindex(columns=params.samples).fillna(0)
		ranking_samples=[sample for sample in params.samples if sample != params.negative_control]
		ranked=rpkm[ranking_samples]
		top=pd.DataFrame(index=rpkm.index)
		top["mean_RPKM_raw"]=ranked.mean(axis=1).fillna(0)
		top["max_RPKM_raw"]=ranked.max(axis=1).fillna(0)
		top["highest_abundance_sample"]=ranked.idxmax(axis=1) if ranking_samples else "not reported"
		top["contig_id"]=top.index
		top=top.loc[top["mean_RPKM_raw"] > 0].sort_values(
			["mean_RPKM_raw", "contig_id"], ascending=[False, True], kind="stable")
		# Include each non-NC sample's most abundant contigs, counting shared contigs once.
		sample_top_ids=set()
		for sample in ranking_samples:
			sample_ranking=top[["contig_id"]].copy()
			sample_ranking["RPKM_raw"]=ranked[sample]
			sample_ranking=sample_ranking.loc[sample_ranking["RPKM_raw"] > 0].sort_values(
				["RPKM_raw", "contig_id"], ascending=[False, True], kind="stable")
			sample_top_ids.update(sample_ranking.head(max(0, params.top_per_sample)).index)
		# Union with the overall mean-RPKM selection; neither selection caps the other.
		mean_top_ids=set(top.head(max(0, params.top_n)).index)
		top=top.loc[top.index.isin(sample_top_ids | mean_top_ids)].copy()
		top["selected_by_sample_top"]=top.index.isin(sample_top_ids)
		top["selected_by_mean_top"]=top.index.isin(mean_top_ids)
		top.insert(0, "rank", range(1, len(top) + 1))
		provenance=pd.read_csv(input.provenance, sep="\t", usecols=["contig_id", "sample", "assembler", "length_bp"])
		origins=provenance.set_index("contig_id")
		for field in ["sample", "assembler", "length_bp"]:
			top["assembly_sample" if field == "sample" else field]=origins[field].reindex(top.index)
		clusters=pd.read_csv(input.clusters, sep="\t", header=None, names=["representative_id", "member_id"])
		clusters=clusters.drop_duplicates()
		top["cluster_size"]=clusters.groupby("representative_id").size().reindex(top.index).fillna(1).astype(int)
		for path, prefix in [(input.rpkm, "RPKM_raw"), (input.counts, "mapped_reads"),
				(input.breadth, "breadth_percent"), (input.depth, "mean_depth")]:
			matrix=pd.read_csv(path, index_col=0).reindex(index=top.index, columns=params.samples).fillna(0)
			for sample in params.samples:
				top[prefix + "_" + sample]=matrix[sample]
		for sample, path in zip(params.samples, input.unique):
			unique=pd.read_csv(path, sep="\t", usecols=[0, 4]).set_index("Contig")
			top["unique_mapped_reads_" + sample]=unique.iloc[:, 0].reindex(top.index).fillna(0)
			top["unique_mapping_ratio_" + sample]=(top["unique_mapped_reads_" + sample] /
				top["mapped_reads_" + sample].replace(0, float("nan"))).fillna(0)
		membership=clusters.loc[clusters["representative_id"].isin(top.index)].merge(
			provenance.rename(columns={"contig_id": "member_id", "sample": "assembly_sample"}), on="member_id", how="left")
		membership["rank"]=membership["representative_id"].map(top["rank"])
		membership=membership.sort_values(["rank", "member_id"], kind="stable")
		membership.to_csv(output.membership, sep="\t", index=False)

		def records(path):
			name, chunks=None, []
			with open(path) as handle:
				for line in handle:
					if line.startswith(">"):
						if name is not None:
							yield name, "".join(chunks)
						name, chunks=line[1:].split()[0], []
					else:
						chunks.append(line.strip())
			if name is not None:
				yield name, "".join(chunks)

		selected={name: sequence for name, sequence in records(input.fasta) if name in top.index}
		with open(output.fasta, "w") as handle:
			for name in top.index:
				sequence=selected[name]
				handle.write(f">{name}\n{sequence}\n")
				top.loc[name, "gc_percent"]=100 * sum(sequence.upper().count(base) for base in "GC") / len(sequence)
		top.to_csv(output.ranking, sep="\t", index=False)

rule collect_all_assembled_top_evidence:
	input:
		ranking=lambda wc: ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_provenance.tot.tsv" if wc.catalogue == "phage_isolates_contigs" else ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_contigs_tot.tsv",
		genomad=lambda wc: expand(dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_tot", sample=SAMPLES) + (expand(RNA_DIR + "/{sample}/{assembler}_genomad", sample=SAMPLES, assembler=RNA_ASSEMBLERS) if RNA_MODE and wc.catalogue == "all_assembled_top_contigs" else []),
		circularity=lambda wc: expand(dirs_dict["VIRAL_DIR"] + "/{sample}_{assembler}_circularity.tot.tsv", sample=SAMPLES, assembler=["spades"] + (RNA_ASSEMBLERS if RNA_MODE and wc.catalogue == "all_assembled_top_contigs" else [])),
	output:
		tsv=ALL_ASSEMBLED_MAPPING_DIR + "/{catalogue}_existing_annotations_tot.tsv",
	params:
		prefix=lambda wc: [sample + "_spades_filtered_scaffolds.tot" for sample in SAMPLES] + ([sample + "_" + assembler for sample in SAMPLES for assembler in RNA_ASSEMBLERS] if RNA_MODE and wc.catalogue == "all_assembled_top_contigs" else []),
	wildcard_constraints:
		catalogue="all_assembled_top_contigs|phage_isolates_contigs",
	message:
		"Collecting existing per-assembly geNomad and terminal-repeat evidence using contig identifiers"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/collect_all_assembled_top_evidence/catalogue={catalogue}.tsv"
	threads: 1
	run:
		import csv
		from pathlib import Path

		def rows(path):
			if path.is_file() and path.stat().st_size:
				with path.open() as handle:
					yield from csv.DictReader(handle, delimiter="\t")

		circles={}
		for path in input.circularity:
			for row in rows(Path(path)):
				circles[row["contig_id"]]=row
		circle_fields=["dtr_bp", "itr_bp", "dtr_sequence", "dtr_left_start", "dtr_left_end", "dtr_right_start", "dtr_right_end",
			"itr_left_sequence", "itr_right_sequence", "itr_left_start", "itr_left_end", "itr_right_start", "itr_right_end",
			"terminal_repeat_type", "repeat_warning"]
		fields=["contig_id", "genomad_assessed", "genomad_classification", "genomad_virus_score", "genomad_plasmid_score",
			"genomad_taxonomy", "genomad_topology", "genomad_evidence_scope", "genomad_reported_sequence", "genomad_source",
			"genomad_virus_call", "genomad_plasmid_call", "genomad_contig_virus_score", "genomad_contig_plasmid_score",
			"genomad_contig_chromosome_score", "genomad_score_source", "genomad_reported_virus_scores", "genomad_reported_plasmid_scores",
			"genomad_provirus_sequences", "genomad_provirus_coordinates", "genomad_provirus_scores",
			"genomad_plasmid_reported_sequence", "genomad_plasmid_topology",
			"genomad_hallmarks_assessed",
			"genomad_virus_hallmark_count", "genomad_plasmid_hallmark_count",
			"genomad_virus_hallmark_genes", "genomad_plasmid_hallmark_genes",
			"genomad_virus_hallmark_annotations", "genomad_plasmid_hallmark_annotations", "phage_plasmid_candidate"] + circle_fields
		origins=list(rows(Path(input.ranking)))
		wanted={row["contig_id"] for row in origins}
		origin_samples={row["contig_id"]: row.get("sample") or row.get("assembly_sample", "") for row in origins}
		origin_assemblers={row["contig_id"]: row.get("assembler", "") for row in origins}
		folders=[input.genomad] if isinstance(input.genomad, str) else list(input.genomad)
		prefixes=[params.prefix] if isinstance(params.prefix, str) else list(params.prefix)
		virus_hits, plasmid_hits, scores, sources, hallmarks={}, {}, {}, {}, {}
		for folder, prefix in zip(folders, prefixes):
			folder=Path(folder)
			dna=prefix.endswith("_spades_filtered_scaffolds.tot")
			sample=prefix.removesuffix("_spades_filtered_scaffolds.tot") if dna else prefix.rsplit("_", 1)[0]
			assemblers={"spades", "metaspades"} if dna else {prefix.rsplit("_", 1)[1]}
			summary=folder / (prefix + "_summary")
			virus=summary / (prefix + "_virus_summary.tsv")
			plasmid=summary / (prefix + "_plasmid_summary.tsv")
			for row in rows(virus):
				name=row["seq_name"].split("|provirus_", 1)[0]
				if name in wanted:
					virus_hits.setdefault(name, []).append(row)
			for row in rows(plasmid):
				if row["seq_name"] in wanted:
					plasmid_hits.setdefault(row["seq_name"], []).append(row)
			score_path=folder / (prefix + "_aggregated_classification") / (prefix + "_aggregated_classification.tsv")
			for row in rows(score_path):
				if row["seq_name"] in wanted:
					scores[row["seq_name"]]=dict(row, source=str(score_path))
			# These are geNomad's own hallmark assignments, not new score cutoffs.
			# The whole-contig annotation table also contains genes within provirus regions.
			gene_path=folder / (prefix + "_annotate") / (prefix + "_genes.tsv")
			for row in rows(gene_path):
				name=row["gene"].rsplit("_", 1)[0]
				if name in wanted:
					for biology in ["virus", "plasmid"]:
						if row.get(biology + "_hallmark") == "1":
							hallmarks.setdefault((name, biology), {})[row["gene"]]=row.get("annotation_description", "")
			for name in wanted:
				if origin_samples[name] == sample and origin_assemblers[name] in assemblers:
					sources[name]=(str(folder), virus.is_file(), gene_path.is_file())
		with open(output.tsv, "w") as handle:
			writer=csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
			writer.writeheader()
			for origin in origins:
				name=origin["contig_id"]
				exact=[row for row in virus_hits.get(name, []) if row["seq_name"] == name]
				regional=[row for row in virus_hits.get(name, []) if row["seq_name"] != name]
				plasmids=plasmid_hits.get(name, [])
				hits=exact or regional or plasmids
				result=dict.fromkeys(fields, "not reported")
				folder, assessed, genes_assessed=sources.get(name, ("not reported", False, False))
				result.update(contig_id=name, genomad_assessed="yes" if assessed else "not reported", genomad_source=folder,
					genomad_virus_call=bool(exact or regional), genomad_plasmid_call=bool(plasmids),
					genomad_hallmarks_assessed="yes" if genes_assessed else "not reported")
				if hits:
					result["genomad_classification"]="virus" if exact or regional else "plasmid"
					result["genomad_evidence_scope"]="provirus region" if regional and not exact else "whole contig"
					for column, source in [("genomad_virus_score", "virus_score"), ("genomad_plasmid_score", "plasmid_score"),
							("genomad_taxonomy", "taxonomy"), ("genomad_topology", "topology"), ("genomad_reported_sequence", "seq_name")]:
						values=[row[source] for row in hits if row.get(source)]
						result[column]=" | ".join(values) if values else "not reported"
				for column, source in [("genomad_contig_virus_score", "virus_score"), ("genomad_contig_plasmid_score", "plasmid_score"),
						("genomad_contig_chromosome_score", "chromosome_score"), ("genomad_score_source", "source")]:
					result[column]=scores.get(name, {}).get(source) or "not reported"
				for biology, reported in [("virus", exact or regional), ("plasmid", plasmids)]:
					values=[row[biology + "_score"] for row in reported if row.get(biology + "_score")]
					result["genomad_reported_" + biology + "_scores"]=" | ".join(values) if values else "not reported"
					# Keep existing reported scores, and provide whole-contig scores for non-calls.
					result["genomad_" + biology + "_score"]=" | ".join(values) if values else result["genomad_contig_" + biology + "_score"]
					genes=hallmarks.get((name, biology), {})
					result["genomad_" + biology + "_hallmark_count"]=len(genes) if genes_assessed else "not reported"
					result["genomad_" + biology + "_hallmark_genes"]=" | ".join(genes) or "not reported"
					result["genomad_" + biology + "_hallmark_annotations"]=" | ".join(dict.fromkeys(
						value for value in genes.values() if value and value != "NA")) or "not reported"
				for column, reported, source in [("genomad_provirus_sequences", regional, "seq_name"),
						("genomad_provirus_coordinates", regional, "coordinates"), ("genomad_provirus_scores", regional, "virus_score"),
						("genomad_plasmid_reported_sequence", plasmids, "seq_name"), ("genomad_plasmid_topology", plasmids, "topology")]:
					result[column]=" | ".join(row[source] for row in reported if row.get(source)) or "not reported"
				result["phage_plasmid_candidate"]=(bool(exact or regional) or bool(hallmarks.get((name, "virus")))) and (bool(plasmids) or bool(hallmarks.get((name, "plasmid"))))
				circle=circles.get(name, {})
				for column in circle_fields:
					result[column]=circle.get(column) or "not reported"
				writer.writerow(result)

rule all_assembled_top_metadata:
	input:
		ranking=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_contigs_tot.tsv",
		existing=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_contigs_existing_annotations_tot.tsv",
		virsorter=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + ALL_ASSEMBLED_TOP_NAME + "_tot/final-viral-score.tsv",
		vibrant_quality=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + ALL_ASSEMBLED_TOP_NAME + "_positive_quality.tot.csv",
		vibrant_summary=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + ALL_ASSEMBLED_TOP_NAME + "_summary_results.tot.csv",
		checkv=ALL_ASSEMBLED_TOP_PREFIX + "_checkV/quality_summary.tsv",
		refseq=dirs_dict["ANNOTATION"] + "/blast_output_ViralRefSeq_" + ALL_ASSEMBLED_TOP_NAME + ".tot.csv",
		metavr=dirs_dict["ANNOTATION"] + "/blast_output_METAVR_" + ALL_ASSEMBLED_TOP_NAME + ".tot.csv",
		metavr_metadata=dirs_dict["PLOTS_DIR"] + "/" + ALL_ASSEMBLED_TOP_NAME + "_METAVR.tot/METAVR_main_table_for_hits.tsv",
		pharokka=ALL_ASSEMBLED_TOP_PREFIX + "_pharokka",
		cenote=[ALL_ASSEMBLED_TOP_CENOTE] if MAP_TO_ALL_ASSEMBLED and RUN_CENOTE else [],
		sourmash=[dirs_dict["ANNOTATION"] + "/sourmash_all_assembled_top_contigs_tot.metadata.tsv"] if SOURMASH_CONTIG_CATALOGUE else [],
	output:
		metadata=ALL_ASSEMBLED_DIR + "/all_assembled_top_contigs_metadata_tot.tsv",
		pharokka_cds=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_contigs_pharokka_cds_tot.tsv",
	message:
		"Combining original representative abundance, classification and annotation of every selected cluster"
	params:
		cenote_enabled=MAP_TO_ALL_ASSEMBLED and RUN_CENOTE,
		cenote_min_length=int(config.get("cenote_min_contig_length", 1000)),
	benchmark:
		dirs_dict["BENCHMARKS"] + "/all_assembled_top_metadata/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	run:
		import pandas as pd
		from pathlib import Path

		def table(path, **kwargs):
			path=Path(path)
			return pd.read_csv(path, sep="\t", **kwargs) if path.is_file() and path.stat().st_size else pd.DataFrame()

		metadata=table(input.ranking).set_index("contig_id")
		existing=table(input.existing).set_index("contig_id")
		metadata=metadata.join(existing)
		if getattr(input, "sourmash", []):
			sourmash=table(input.sourmash[0]).set_index("cluster_rep")
			metadata=metadata.join(sourmash.reindex(metadata.index))
		for path, key, prefix in [(input.virsorter, "seqname", "VirSorter2"), (input.checkv, "contig_id", "CheckV")]:
			frame=table(path)
			if key in frame:
				frame=frame.drop_duplicates(key).set_index(key).add_prefix(prefix + "_")
				metadata=metadata.join(frame.reindex(metadata.index))

		if params.cenote_enabled:
			metadata["Cenote_assessed"]="yes"
			metadata.loc[metadata["length_bp"] < params.cenote_min_length, "Cenote_assessed"]="below length cutoff"
			metadata["Cenote_classification"]="not reported"
			cenote=table(Path(input.cenote[0]) / "mosaic_ct3_virus_summary.tsv")
			if "input_name" in cenote:
				# Cenote's internal names are not the assembly catalogue identifiers.
				cenote["contig_id"]=cenote["input_name"].str.split().str[0]
				cenote=cenote.loc[cenote["contig_id"].isin(metadata.index)].copy()
				cenote["evidence_scope"]="whole contig"
				processed=pd.to_numeric(cenote["virus_seq_length"], errors="coerce") != cenote["contig_id"].map(metadata["length_bp"])
				cenote.loc[processed, "evidence_scope"]="processed contig"
				cenote.loc[cenote["contig"].str.contains("@", regex=False), "evidence_scope"]="region"
				metadata.loc[metadata.index.isin(cenote["contig_id"]), "Cenote_classification"]="virus"
				cenote=cenote.drop(columns="input_name").groupby("contig_id", sort=False).agg(lambda values: " | ".join(dict.fromkeys(values.dropna().astype(str))))
				metadata=metadata.join(cenote.add_prefix("Cenote_").reindex(metadata.index))

		def original_name(name):
			if name in metadata.index:
				return name
			for separator in ["_fragment_", "|provirus_", "_provirus_"]:
				if name.split(separator, 1)[0] in metadata.index:
					return name.split(separator, 1)[0]
			return name

		for path, prefix in [(input.vibrant_quality, "VIBRANT"), (input.vibrant_summary, "VIBRANT_summary")]:
			frame=table(path)
			if "scaffold" in frame:
				frame["contig_id"]=frame["scaffold"].map(original_name)
				frame["evidence_scope"]=frame.apply(lambda row: "whole contig" if row["scaffold"] == row["contig_id"] else "region", axis=1)
				frame=frame.groupby("contig_id", sort=False).agg(lambda values: " | ".join(dict.fromkeys(values.dropna().astype(str))))
				frame.columns=[prefix + "_" + column.replace(" ", "_") for column in frame.columns]
				metadata=metadata.join(frame.reindex(metadata.index))

		blast_columns=["qseqid", "sseqid", "description", "qstart", "qend", "qlen", "slen", "qcovs", "evalue", "alignment_length", "pident"]
		for path, prefix in [(input.refseq, "RefSeq"), (input.metavr, "METAVR")]:
			frame=table(path, header=None, names=blast_columns)
			frame=frame.loc[frame["qseqid"] != "qseqid"] if "qseqid" in frame else pd.DataFrame(columns=blast_columns)
			for field in blast_columns[3:]:
				frame[field]=pd.to_numeric(frame[field], errors="coerce")
			frame=frame.sort_values(["evalue", "qcovs", "pident", "alignment_length", "sseqid"], ascending=[True, False, False, False, True], kind="stable")
			frame=frame.drop_duplicates("qseqid").set_index("qseqid").add_prefix(prefix + "_")
			metadata=metadata.join(frame.reindex(metadata.index))
		metavr=table(input.metavr_metadata)
		if "uvig" in metavr:
			metavr=metavr.drop_duplicates("uvig").set_index("uvig")
			for field in metavr.columns:
				metadata["METAVR_" + field]=metadata["METAVR_sseqid"].fillna("").str.split("|").str[0].map(metavr[field])

		cds=table(Path(input.pharokka) / "pharokka_cds_final_merged_output.tsv")
		if "contig" in cds:
			cds=cds.loc[cds["contig"].isin(metadata.index)].copy()
			metadata["Pharokka_CDS_count"]=cds.groupby("contig").size().reindex(metadata.index, fill_value=0)
			if "annot" in cds:
				known=cds.loc[cds["annot"].notna() & ~cds["annot"].str.contains("hypothetical|unknown", case=False, na=False)]
				metadata["Pharokka_annotated_CDS_count"]=known.groupby("contig").size().reindex(metadata.index, fill_value=0)
			if "category" in cds:
				metadata["Pharokka_function_categories"]=cds.groupby("contig")["category"].agg(lambda values: "; ".join(sorted(set(values.dropna()))))
		else:
			cds=pd.DataFrame(columns=["contig", "gene", "annot", "category"])
			metadata["Pharokka_CDS_count"]="not reported"
			metadata["Pharokka_annotated_CDS_count"]="not reported"
			metadata["Pharokka_function_categories"]="not reported"
		cds.to_csv(output.pharokka_cds, sep="\t", index=False)
		metadata=metadata.replace("", "not reported").fillna("not reported")
		metadata.reset_index().to_csv(output.metadata, sep="\t", index=False)


if MAP_TO_REFSEQ:
	rule refseq_detection:
		input:
			fasta=config["RefSeqViral_db"],
			covstats=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_covstats.txt", sample=SAMPLES),
			unique_covstats=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_covstats.txt", sample=SAMPLES),
			basecov=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_basecov.txt.gz", sample=SAMPLES) if PER_BASE_COVERAGE else [],
			unique_basecov=expand(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_basecov.txt.gz", sample=SAMPLES) if PER_BASE_COVERAGE else [],
			rpkm=REFSEQ_MAPPING_DIR + "/RefSeqViral_RPKM_raw_tot.txt",
			counts=REFSEQ_MAPPING_DIR + "/RefSeqViral_counts_raw_tot.txt",
			breadth=REFSEQ_MAPPING_DIR + "/RefSeqViral_breadth_coverage_percent_tot.txt",
			depth=REFSEQ_MAPPING_DIR + "/RefSeqViral_mean_depth_tot.txt",
		output:
			long=REFSEQ_MAPPING_DIR + "/RefSeqViral_detection_long_tot.tsv",
			matrix=REFSEQ_MAPPING_DIR + "/RefSeqViral_detection_matrix_tot.tsv",
			metadata=REFSEQ_MAPPING_DIR + "/RefSeqViral_detection_metadata_tot.tsv",
			high=REFSEQ_MAPPING_DIR + "/RefSeqViral_detected_high_confidence_tot.tsv",
			candidates=REFSEQ_MAPPING_DIR + "/RefSeqViral_detected_candidates_tot.tsv",
			contaminants=REFSEQ_MAPPING_DIR + "/RefSeqViral_detected_contaminants_tot.tsv",
			figures=directory(dirs_dict["PLOTS_DIR"] + "/07_RefSeq_detection"),
		params:
			samples=SAMPLES,
			metadata_cache=REFSEQ_MAPPING_DIR + "/RefSeqViral_NCBI_metadata_cache.sqlite",
			metadata_refresh=config_bool("refseq_metadata_refresh", False),
			negative_control=str(config.get("negative_control", "")).strip(),
			high_min_length_bp=int(config.get("refseq_high_min_length_bp", 6667)),
			high_min_covered_bases=int(config.get("refseq_high_min_covered_bases", 5000)),
			high_min_breadth_percent=float(config.get("refseq_high_min_breadth_percent", 75)),
			high_min_unique_reads=int(config.get("refseq_high_min_unique_reads", 5)),
			candidate_min_unique_reads=int(config.get("refseq_candidate_min_unique_reads", 10)),
			candidate_min_breadth_percent=float(config.get("refseq_candidate_min_breadth_percent", 25)),
			ambiguous_min_reads=int(config.get("refseq_ambiguous_min_reads", 5)),
			ambiguous_max_unique_ratio=float(config.get("refseq_ambiguous_max_unique_ratio", 0.2)),
			negative_control_min_reads=int(config.get("refseq_negative_control_min_reads", 5)),
			negative_control_max_enrichment=float(config.get("refseq_negative_control_max_enrichment", 2)),
			enrichment_pseudocount=float(config.get("refseq_enrichment_pseudocount", 0.01)),
			plot_max_accessions=int(config.get("refseq_plot_max_accessions", 60)),
		message:
			"Classifying RefSeq viral read-mapping evidence"
		benchmark:
			dirs_dict["BENCHMARKS"] + "/refseq_detection/tot.tsv"
		threads: 1
		resources:
			mem_mb=8000,
		log:
			notebook=dirs_dict["NOTEBOOKS_DIR"] + "/07_RefSeq_detection.py.ipynb"
		notebook:
			dirs_dict["RAW_NOTEBOOKS"] + "/07_RefSeq_detection.py.ipynb"

def input_QC_long_only_nanopore_pre(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanostats_preQC.html", sample_nanopore=NANOPORE_SAMPLES))
	return(input_list)

def input_QC_long_only_nanopore_post(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanostats_postQC.html", sample_nanopore=NANOPORE_SAMPLES))
	return(input_list)

def input_QC_long_only_pacbio_pre(wildcards):
	input_list=[]
	if PACBIO:
		input_list.extend(expand(dirs_dict["QC_DIR"] + "/{sample_pacbio}_pacbio_nanostats_preQC.html", sample_pacbio=PACBIO_SAMPLES))
	return(input_list)

def input_QC_long_only_pacbio_post(wildcards):
	input_list=[]
	if PACBIO:
		input_list.extend(expand(dirs_dict["QC_DIR"] + "/{sample_pacbio}_pacbio_nanostats_postQC_{sampling}.html", sample_pacbio=PACBIO_SAMPLES, sampling=wildcards.sampling))
	return(input_list)


rule QC_long_only_parsing:
	input:
		nanopore_pre=input_QC_long_only_nanopore_pre,
		nanopore_post=input_QC_long_only_nanopore_post,
		pacbio_pre=input_QC_long_only_pacbio_pre,
		pacbio_post=input_QC_long_only_pacbio_post
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_summary.{sampling}.html",
		read_count_png=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_count.{sampling}.png",
		read_count_svg=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_count.{sampling}.svg",
		total_bases_png=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_total_bases.{sampling}.png",
		total_bases_svg=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_total_bases.{sampling}.svg",
		read_length_png=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_length.{sampling}.png",
		read_length_svg=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_length.{sampling}.svg",
		read_quality_png=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_quality.{sampling}.png",
		read_quality_svg=dirs_dict["PLOTS_DIR"] + "/01_QC_long_only_read_quality.{sampling}.svg"
	params:
		sampling="{sampling}",
		samples_nanopore=NANOPORE_SAMPLES,
		samples_pacbio=PACBIO_SAMPLES
	benchmark:
		dirs_dict["BENCHMARKS"] + "/QC_long_only_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/01_QC_long_only.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/01_QC_long_only.py.ipynb"

rule assembly_long_only_parsing:
	input:
		quast_report_dir=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}"
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_summary.{sampling}.html",
		contig_number_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_contig_number.{sampling}.png",
		contig_number_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_contig_number.{sampling}.svg",
		total_length_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_total_length.{sampling}.png",
		total_length_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_total_length.{sampling}.svg",
		n50_png=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_N50.{sampling}.png",
		n50_svg=dirs_dict["PLOTS_DIR"] + "/03_assembly_long_only_N50.{sampling}.svg"
	params:
		input_quast_report=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
		sampling="{sampling}"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assembly_long_only_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/03_assembly_long_only.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/03_assembly_long_only.py.ipynb"

def input_bacterial_results_checkm(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if ISOLATES:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}", sample=SAMPLES, sampling=wildcards.sampling))
	return(input_list)

def input_bacterial_results_sourmash(wildcards):
	input_list=[]
	if NANOPORE & NANOPORE_ONLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore.classifications.csv", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if NANOPORE & (not NANOPORE_ONLY) & PAIRED:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_hybrid.classifications.csv", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO & PACBIO_ONLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio.classifications.csv", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if PACBIO & PACBIO_HYBRID:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_hybrid.classifications.csv", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if ISOLATES:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}.classifications.csv", sample=SAMPLES, sampling=wildcards.sampling))
	return(input_list)

def input_bacterial_results_coverage(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["MAPPING_DIR"] + "/{sample}_{sampling}_long_read_contig_coverage.tsv", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO:
		input_list.extend(expand(dirs_dict["MAPPING_DIR"] + "/{sample}_{sampling}_long_read_contig_coverage.tsv", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	return(input_list)


def input_bacterial_results_genomad(wildcards):
	input_list=[]
	if NANOPORE:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"] + "/{sample}_long_geNomad_{sampling}", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"] + "/{sample}_long_geNomad_{sampling}", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if ISOLATES:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}", sample=SAMPLES, sampling=wildcards.sampling))
	return(input_list)


rule selectMetaVRMetadata:
	input:
		blast=lambda wc: dirs_dict["ANNOTATION"] + "/blast_output_METAVR_" + ("filtered_" + REPRESENTATIVE_CONTIGS_BASE if wc.report == "08_METAVR_analysis" else wc.report.removesuffix("_METAVR")) + "." + wc.sampling + ".csv",
		metadata=os.path.join(config["METAVR_db"], "METAVR_main_table.parquet"),
	output:
		selected=dirs_dict["PLOTS_DIR"] + "/{report}.{sampling}/METAVR_main_table_for_hits.tsv",
	message:
		"Selecting MetaVR metadata for the matched viral genomes"
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/selectMetaVRMetadata/report={report}__sampling={sampling}.tsv"
	threads: 1
	wildcard_constraints:
		report="[^/]+",
		sampling="tot|sub"
	resources:
		mem_mb=8000,
	shell:
		"""
		set -euo pipefail
		mkdir -p $(dirname {output.selected:q})
		python - {input.blast:q} {input.metadata:q} {output.selected:q} <<'PYCODE'
import csv
import sys
from fastparquet import ParquetFile

columns = ["uvig", "taxon_oid", "quality", "completeness", "genome_type",
           "ictv_taxonomy", "ictv_taxonomy_method", "host_taxonomy", "host_taxonomy_method"]
with open(sys.argv[1], newline="") as handle:
    subjects = set()
    for row in csv.reader(handle, delimiter="\t"):
        if len(row) > 1 and row[1] != "sseqid":
            subjects.add(row[1].split("|")[0])
with open(sys.argv[3], "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(columns)
    if subjects:
        parquet = ParquetFile(sys.argv[2])
        for group in parquet.iter_row_groups(columns=columns):
            selected = group.loc[group["uvig"].isin(subjects), columns]
            writer.writerows(selected.itertuples(index=False, name=None))
PYCODE
		"""


rule METAVR_analysis:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		blast=dirs_dict["ANNOTATION"] + "/blast_output_METAVR_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.csv",
		checkv=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		uvig_metadata=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/METAVR_main_table_for_hits.tsv",
		source_metadata=os.path.join(config["METAVR_db"], "IMG_full_metadata.tsv.gz"),
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}.html",
		hit_pairs=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/hit_pairs.tsv",
		query_summary=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/vOTU_summary.tsv",
		sample_summary=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/sample_summary.tsv",
		environment_counts=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/environment_counts.tsv",
		viral_taxonomy_counts=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/viral_taxonomy_counts.tsv",
		host_taxonomy_counts=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/host_taxonomy_counts.tsv",
		host_method_counts=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/host_method_counts.tsv",
		metadata_subset=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/METAVR_metadata_for_hits.tsv",
		metadata_manifest=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/METAVR_metadata_for_hits.manifest.json",
		run_summary=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/run_summary.tsv",
		provenance=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/provenance.json",
		composition_png=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/coverage_composition.png",
		composition_svg=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/coverage_composition.svg",
		samples_png=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/sample_composition.png",
		samples_svg=dirs_dict["PLOTS_DIR"] + "/08_METAVR_analysis.{sampling}/sample_composition.svg",
	params:
		samples=SAMPLES,
		sampling="{sampling}",
		metadata_chunksize=250000,
	message:
		"Summarizing MetaVR hit environments, viral taxonomy and recorded host taxonomy"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/METAVR_analysis/sampling={sampling}.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/08_METAVR_analysis.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_METAVR_analysis.py.ipynb"


rule bacterial_results_parsing:
	input:
		checkm=input_bacterial_results_checkm,
		sourmash=input_bacterial_results_sourmash,
		gtdbtk=dirs_dict["ASSEMBLY_DIR"] + "/assembly_bacteria_GTDB-Tk_{sampling}",
		quast=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
		coverage=input_bacterial_results_coverage,
		genomad=input_bacterial_results_genomad,
	output:
		summary_html=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_summary.{sampling}.html",
		summary_csv=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_summary.{sampling}.csv",
		checkm_png=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_checkm.{sampling}.png",
		checkm_svg=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_checkm.{sampling}.svg",
		taxonomy_png=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_taxonomy.{sampling}.png",
		taxonomy_svg=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_taxonomy.{sampling}.svg",
		quality_taxonomy_png=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_quality_taxonomy.{sampling}.png",
		quality_taxonomy_svg=dirs_dict["PLOTS_DIR"] + "/06_bacterial_results_quality_taxonomy.{sampling}.svg"
	params:
		sampling="{sampling}",
		long_assembler=LONG_ASSEMBLER,
		long_assembler_pacbio=LONG_ASSEMBLER_PACBIO,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/bacterial_results_parsing/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/06_bacterial_results.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/06_bacterial_results.py.ipynb"
