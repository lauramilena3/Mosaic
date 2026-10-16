#ruleorder: shortReadAsemblySpadesPE > shortReadAsemblySpadesSE
rule rename_assembled_contigs:
	input:
		fasta=RESULTS_DIR + "/{assembly_path}.unrenamed.fasta",
	output:
		fasta=RESULTS_DIR + "/{assembly_path}.fasta",
		ids=RESULTS_DIR + "/{assembly_path}.ids.tsv",
	params:
		spades_assembler="metaspades" if METAGENOME_FLAG == "--meta" else "spades",
		min_length=lambda wc: int(config.get("rna_min_contig_length", 500)) if wc.assembly_path.startswith("03_CONTIGS/RNA/") else 0,
	message:
		"Naming assembled contigs and recording their original identifiers"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rename_assembled_contigs/assembly_path={assembly_path}.tsv"
	threads: 1
	wildcard_constraints:
		assembly_path=r"03_CONTIGS/(?:RNA/[^/]+_(?:rnaviralspades|megahit|trinity)|[^/]+_spades_filtered_scaffolds\.(?:tot|sub)|[^/]+_(?:contigs_[^/]+|corrected_scaffolds_pilon)\.(?:tot|sub))|08_ASSEMBLY_TEST/[^/]+_metaspades_filtered_scaffolds\.(?:tot|sub)"
	run:
		import csv
		import re
		from Bio import SeqIO

		path=wildcards.assembly_path
		if path.startswith("03_CONTIGS/RNA/"):
			sample, assembler=os.path.basename(path).rsplit("_", 1)
			prefix=sample
		else:
			stem, sampling=os.path.basename(path).rsplit(".", 1)
			if path.startswith("08_ASSEMBLY_TEST/"):
				sample, percentage=stem.removesuffix("_metaspades_filtered_scaffolds").rsplit("_", 1)
				assembler="metaspades"
				prefix=sample + "_" + percentage + "pct"
			elif stem.endswith("_spades_filtered_scaffolds"):
				sample=stem.removesuffix("_spades_filtered_scaffolds")
				assembler=params.spades_assembler
				prefix=sample
			elif stem.endswith("_corrected_scaffolds_pilon"):
				sample, assembler=stem.removesuffix("_corrected_scaffolds_pilon").rsplit("_", 1)
				prefix=sample + "_pilon"
			else:
				sample, assembler=stem.rsplit("_contigs_", 1)
				stage=re.match(r"^(medaka_polished|polypolish|racon|pilon_[1-4]_polished)_(.+)$", sample)
				if stage:
					sample=stage.group(2)
					if stage.group(1) == "racon":
						iteration, assembler=assembler.split("_", 1)
						prefix=sample + "_racon" + iteration
					else:
						prefix=sample + "_" + stage.group(1)
				else:
					prefix=sample
			if sampling == "sub":
				prefix += "_sub"

		with open(input.fasta) as source:
			if params.min_length:
				total_contigs=sum(1 for record in SeqIO.parse(source, "fasta") if len(record.seq) >= params.min_length)
			else:
				total_contigs=sum(1 for line in source if line.startswith(">"))
		width=max(5, len(str(total_contigs)))

		with open(input.fasta) as source, open(output.fasta, "w") as fasta, open(output.ids, "w") as handle:
			writer=csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["contig_id", "sample", "assembler", "original_id", "length_bp", "original_description"])
			number=0
			for record in SeqIO.parse(source, "fasta"):
				if len(record.seq) < params.min_length:
					continue
				number += 1
				original_id, description=record.id, record.description
				contig_id=f"{prefix}_{assembler}_{number:0{width}d}_len_{len(record.seq)}"
				writer.writerow([contig_id, sample, assembler, original_id, len(record.seq), description])
				record.id=record.name=contig_id
				record.description=""
				SeqIO.write(record, fasta, "fasta-2line")

def input_error_correction(wildcards):
	params_ecc=""
	if wildcards.sample=="ALL":
		params_ecc="--only-assembler"
	return params_ecc

def input_threads_assembler(wildcards):
	use_threads=12
	if wildcards.sample=="ALL":
		use_threads=64
	return use_threads

rule shortReadAsemblySpadesPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq.gz",
	output:
		scaffolds=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.unrenamed.fasta"),
		assembly_graph=dirs_dict["ASSEMBLY_DIR"] +"/{sample}_assembly_graph_spades.{sampling}.fastg",
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/scaffolds.fasta",
		assembly_graph=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/assembly_graph.fastg",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}"),
		metagenomic_flag=METAGENOME_FLAG,
		error_correction=input_error_correction,
		mem_gb=ASSEMBLY_MEM_GB,
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/filtered_list.txt"),
	message:
		"Assembling PE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/shortReadAsemblySpadesPE/sample={sample}__sampling={sampling}.tsv"
	log:
		dirs_dict["ASSEMBLY_DIR"] + "/LOGS/{sample}_spades.{sampling}.log"
	threads: input_threads_assembler
	resources:
		mem_mb=ASSEMBLY_MEM_MB
	priority: 1
	wildcard_constraints:
		sampling="tot|sub"  
	shell:
		"""
		mkdir -p "$(dirname {log:q})"
		rm -rf {params.assembly_dir}
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		{params.metagenomic_flag} -t {threads} --memory {params.mem_gb} {params.error_correction} 2>&1 | tee {log:q}
		cp {params.assembly_dir:q}/spades.log {log:q}
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {params.filtered_list}
		seqtk subseq {params.raw_scaffolds} {params.filtered_list} > {output.scaffolds}
		cp {params.assembly_graph} {output.assembly_graph}
		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds}
		rm -rf {params.assembly_dir}
		"""

def input_Quast(wildcards):
	input_list=[]
	if NANOPORE & NANOPORE_ONLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_" + LONG_ASSEMBLER + ".{sampling}.fasta", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if NANOPORE & (not NANOPORE_ONLY) & PAIRED:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_" + LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}.fasta", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO & PACBIO_ONLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_" + LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if PACBIO & PACBIO_HYBRID:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/polypolish_{sample}_contigs_" + LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["ASSEMBLY_DIR"] + "/ALL_spades_filtered_scaffolds." + wildcards.sampling + ".fasta")
	if PAIRED & (not CROSS_ASSEMBLY):
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta", sample=SAMPLES, sampling=wildcards.sampling))
	return(input_list)

rule assemblyStats:
	input:
		scaffolds=input_Quast,
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}"),
		quast_txt=dirs_dict["ASSEMBLY_DIR"] + "/assembly_quast_report.{sampling}.txt",
		quast_tsv=dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}/transposed_report.tsv",
	message:
		"Creating assembly stats with quast"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assemblyStats/sampling={sampling}.tsv"
	threads: 1
	resources:
		mem_mb=QUAST_MEM_MB
	shell:
		"""
		quast.py {input.scaffolds} -o {output.quast_report_dir}
		cp {output.quast_report_dir}/report.txt {output.quast_txt}
		"""
