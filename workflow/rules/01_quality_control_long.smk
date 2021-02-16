rule qualityCheckNanopore:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"]+"/{sample_nanopore}_nanopore.fastq",
	output:
		nanoqc_dir=temp(directory(dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanoplot")),
		nanoqc=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_preQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
#	threads: 1
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.raw_fastq}
		mv {output.nanoqc_dir}/nanoQC.html {output.nanoqc}
		"""

rule remove_adapters_quality_nanopore:
	input:
		raw_data=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore.fastq",
	output:
		trimmed_data=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq",
		porechopped=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_porechopped.fastq"),
	params:
		headcrop=50,
		tailcrop=50,
		quality=config['nanopore_quality'],
	message:
		"Trimming Nanopore Adapters with Porechop"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env3.yaml"
	threads: 2
	shell:
		"""
		porechop -i {input.raw_data} -o {output.porechopped} --threads {threads}
		NanoFilt -q {params.quality} -l 1000 --headcrop {params.headcrop} --tailcrop {params.tailcrop} {output.porechopped} > {output.trimmed_data}
		"""

rule remove_contaminants_nanopore:
	input:
		trimmed_data=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq",
		contaminants_fasta=expand(dirs_dict["CONTAMINANTS_DIR"] +"/{contaminants}.fasta",contaminants=CONTAMINANTS),
	output:
		fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.fastq"),
		size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.txt",
		phix_contaminants_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{sample_nanopore}_nanopore_contaminants.fasta",
	message:
		"Remove contamination with Porechop"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 2
	shell:
		"""
		cat {input.contaminants_fasta} > {output.phix_contaminants_fasta}
		minimap2 -ax map-ont {output.phix_contaminants_fasta} {input.trimmed_data} | samtools fastq -n -f 4 - > {output.fastq}
		grep -c "^@" {output.fastq} > {output.size}
		"""


rule postQualityCheckNanopore:
	input:
		fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.fastq"),
	output:
		nanoqc_dir=temp(directory(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanoqc_post")),
		nanoplot_dir=temp(directory(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanoplot_post")),
		nanoqc=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_postQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
#	threads: 1
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.fastq}
		NanoPlot --fastq {input.fastq} -o {output.nanoplot_dir}
		mv {output.nanoqc_dir}/nanoQC.html {output.nanoqc}
		"""
rule subsampleReadsNanopore:
	input:
		nano_sizes=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.txt", sample_nanopore=NANOPORE_SAMPLES),
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.tot.fastq",
	output:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.sub.fastq",
		size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.sub.txt",
	message:
		"Subsampling Nanopore reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	params:
		sizes=dirs_dict["CLEAN_DATA_DIR"] + "/*_nanopore_clean.tot.txt"
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		nanopore=$( cat {params.sizes} | sort -n | head -1 )
		reformat.sh in={input.nanopore} out={output.nanopore} reads=$nanopore ignorebadquality
		grep -c "^@" {output.nanopore} > {output.size}
		"""
