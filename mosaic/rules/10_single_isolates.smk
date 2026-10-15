rule estimateGenomeCompletnessIsolates:
	input:
		assembly_fasta=(dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_spades_filtered_scaffolds.{sampling}.fasta"),
		checkv_db=(config['checkv_db']),
	output:
		quality_summary=dirs_dict["ASSEMBLY_DIR"] + "/checkV_isolates_{sample}_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["ASSEMBLY_DIR"] + "/checkV_isolates_{sample}_{sampling}/completeness.tsv",
		contamination=dirs_dict["ASSEMBLY_DIR"] + "/checkV_isolates_{sample}_{sampling}/contamination.tsv",
	params:
		checkv_outdir=dirs_dict["ASSEMBLY_DIR"] + "/checkV_isolates_{sample}_{sampling}",
		checkv_db=dirs_dict["ASSEMBLY_DIR"] + "/checkV_isolates_{sample}_{sampling}",
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/estimateGenomeCompletnessIsolates/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		rm -rf {params.checkv_outdir}
		if [ -s {input.assembly_fasta:q} ]; then
		checkv contamination {input.assembly_fasta} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv completeness {input.assembly_fasta} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv complete_genomes {input.assembly_fasta} {params.checkv_outdir}
		checkv quality_summary {input.assembly_fasta} {params.checkv_outdir}
		else
			mkdir -p {params.checkv_outdir:q}
			printf 'contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n' > {output.quality_summary:q}
			printf 'contig_id\\n' > {output.completeness:q}
			printf 'contig_id\\n' > {output.contamination:q}
		fi
		"""

rule filter_isolates:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		quality_summary=dirs_dict["vOUT_DIR"] + "/checkV_isolates_{sampling}/quality_summary.tsv",
	output:
		filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_list.{sampling}.txt",
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates.{sampling}.fasta",
	message:
		"Filtering vOTUs "
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/filter_isolates/sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		grep "High-quality" {input.quality_summary} | cut -f1 > {output.filtered_list}
		seqtk subseq {input.positive_contigs} {output.filtered_list} > {output.filtered_positive_contigs}
		"""

rule makeblastdb_isolates:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
	output:
		orfs_blast_db=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta.pdb",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/makeblastdb_isolates/sampling={sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 8
	shell:
		"""
		makeblastdb -in {input.filtered_positive_contigs} -dbtype prot
		"""

rule blastall_isolates:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
		blast=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta.pdb",
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/isolates_blastall.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/blastall_isolates/sampling={sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 8
	shell:
		"""
		blastp -num_threads {threads} -db {input.filtered_positive_contigs} -query {input.filtered_positive_contigs} \
			-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs qcovhsp length pident evalue positive gaps" > {output.blast_output}
		"""

def inputDatabase(wildcards):
	if wildcards.db=="RefSeqViral":
		return config["RefSeqViral_protein_db"]
	elif wildcards.db=="GenBankViral":
		return config["GenBankViral_protein_db"]
	elif wildcards.db=="METAVR":
		return config["METAVR_protein_db"]

rule blastp_Reference_db:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
		blast_database=inputDatabase
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/isolates_blast_{db}.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/blastp_Reference_db/db={db}__sampling={sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 32
	shell:
		"""
		blastp -num_threads {threads} -db {input.blast_database} -query {input.filtered_positive_contigs} \
			-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs qcovhsp length pident evalue positive gaps" > {output.blast_output}
		"""

rule blastp_database_lengths:
	input:
		blast="{fasta}.faa",
	output:
		length_temp=temp("{fasta}_length_temp.txt"),
		names_temp=temp("{fasta}_names_temp.txt"),
		length="{fasta}_lengths_and_names.txt",
		cummulative_length="{fasta}_coding_lengths.txt",
	message:
		"Getting database ORF lengths"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/blastp_database_lengths/fasta={fasta}.tsv"
	threads: 1
	shell:
		"""
		cat {input.blast} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
			$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.length_temp}
		cat {output.length_temp} | cut -f1 | cut -f1 -d' ' | rev | cut -d_ -f2- | rev > {output.names_temp}
		paste {output.length_temp} {output.names_temp} > {output.length}
		cat {output.length} | cut -f3,2 | awk '{{arr[$2]+=$1}} END {{for (i in arr) {{print i,arr[i]}}}}' > {output.cummulative_length}
		"""

def inputDatabaseLen(wildcards):
	if wildcards.db=="RefSeqViral":
		return (config["RefSeqViral_protein_db"].split(".faa")[0] + "_lengths_and_names.txt")
	elif wildcards.db=="GenBankViral":
		return (config["GenBankViral_protein_db"].split(".faa")[0] + "_lengths_and_names.txt")
	elif wildcards.db=="METAVR":
		return (config["METAVR_protein_db"].split(".faa")[0] + "_lengths_and_names.txt")

def inputDatabaseCoding(wildcards):
	if wildcards.db=="RefSeqViral":
		return (config["RefSeqViral_protein_db"].split(".faa")[0] + "_coding_lengths.txt")
	elif wildcards.db=="GenBankViral":
		return (config["GenBankViral_protein_db"].split(".faa")[0] + "_coding_lengths.txt")
	elif wildcards.db=="METAVR":
		return (config["METAVR_protein_db"].split(".faa")[0] + "_coding_lengths.txt")

rule get_relatives_list:
	input:
		blast_output=dirs_dict["ANNOTATION"] + "/isolates_blast_{db}.{sampling}.csv",
		isolate_len_file=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs_length.{sampling}.txt",
		isolate_coding=dirs_dict["ANNOTATION"]+ "/filtered_isolates_ORFs_coding_lengths.{sampling}.txt",
		reference_len_file=inputDatabaseLen,
		reference_coding=inputDatabaseCoding,
	output:
		blast_relatives_phages=(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_phages_{db}.{sampling}.txt"),
		blast_relatives_proteins=(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_ORFs_{db}.{sampling}.txt"),
	message:
		"Extracting relative sequences NAMES"
	params:
		similarity_cutoff=20,
		length_cutoff=20,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_relatives_list/db={db}__sampling={sampling}.tsv"
	threads: 1
	run:
		# ALL TOGETHER
		import pandas as pd
		from Bio import SeqIO
		import seaborn as sns
		import matplotlib.pyplot as plt
		import time
		#------------------------
		# DEFINE CONSTANTS
		#------------------------
		#input
		blast_out_file=input.blast_output
		isolate_len_file=input.isolate_len_file
		isolates_coding_len=input.isolate_coding
		reference_len_file=input.reference_len_file
		reference_coding_len=input.reference_coding
		#output
		positive_phages=output.blast_relatives_phages
		positive_orfs=output.blast_relatives_proteins

		similarity_cutoff=params.similarity_cutoff
		len_cutoff=params.length_cutoff
		start = time.time()

		#------------------------
		# READ CODING LENGTH
		#------------------------
		print("Calculate coding length")
		print(time.time() - start)

		isolates_len_df=pd.read_csv(isolate_len_file, sep="\t", names=["orf", "aa_length", "phage"]).dropna()
		isolates_coding_len_df=pd.read_csv(isolates_coding_len, sep=" ", names=["phage", "aa_length"]).dropna()

		Reference_len_df=pd.read_csv(reference_len_file, sep="\t", names=["orf", "aa_length", "phage"]).dropna()
		Reference_coding_len_df=pd.read_csv(reference_coding_len, sep=" ", names=["phage", "aa_length"]).dropna()

		isolates_coding_len_df

		#-------------------------
		# PARSE BLAST RESULTS
		#-------------------------
		print("PARSE BLAST RESULTS")
		print(time.time() - start)

		#BLASTALL
		blast_results_reference=pd.read_csv(blast_out_file, names=["qseqid", "sseqid", "qstart", "qend", "qlen", "slen", "qcovs", "qcovhsp", "length", "pident", "evalue", "positive", "gaps"], sep="\t")
		blast_results_reference["phage_A"]=blast_results_reference["qseqid"].str.rsplit("_",1).str[0]
		blast_results_reference["phage_B"]=blast_results_reference["sseqid"].str.rsplit("_",1).str[0]
		blast_results_reference["match_positions"]=blast_results_reference["pident"]*blast_results_reference["length"]/100

		#-------------------------
		# CALCULATE BLAST IDENTITY
		#-------------------------
		print("CALCULATE BLAST IDENTITY")
		print(time.time() - start)

		grouped = blast_results_reference.groupby(by=["phage_A", "phage_B"])
		phage_A=[]
		phage_B=[]
		positive=[]
		for name, group in grouped:
			 best_hit=group.groupby(by=["qseqid"]).first()
			 phage_A.append(name[0])
			 phage_B.append(name[1])
			 positive.append(best_hit["match_positions"].sum())

		similarity_df=pd.DataFrame()
		similarity_df["Phage_A"]=phage_A
		similarity_df["Phage_B"]=phage_B
		similarity_df["Positive"]=positive


		similarity_df1=similarity_df.merge(isolates_coding_len_df, left_on="Phage_A", right_on="phage").drop(columns=["phage"])
		similarity_df1 = similarity_df1.rename(columns={'aa_length': 'lenA'})
		similarity_df2=similarity_df1.merge(Reference_coding_len_df, left_on="Phage_B", right_on="phage").drop(columns=["phage"])
		similarity_df2 = similarity_df2.rename(columns={'aa_length': 'lenB'})
		similarity_df2["simAB"]=similarity_df2["Positive"]*100/similarity_df2["lenA"]
		similarity_df2["ratioAB"]=similarity_df2["lenA"]*100/similarity_df2["lenB"]

		similarity_df2
		#-------------------------
		#FILTER POSITIVE HITS
		#-------------------------
		print("FILTER POSITIVE HITS")
		print(time.time() - start)

		similarity_df2=similarity_df2[similarity_df2["simAB"]>similarity_cutoff]
		similarity_df2=similarity_df2[(similarity_df2["ratioAB"]>(100-len_cutoff)) & (similarity_df2["ratioAB"]<(100*100/(100-len_cutoff)))]

		positive_hits_nucleotide=list(set(similarity_df2["Phage_B"].to_list()))
		positive_hits_nucleotide
		positive_hits_protein=list(set(Reference_len_df[Reference_len_df["phage"].isin(positive_hits_nucleotide)]["orf"].to_list()))
		positive_hits_protein

		#-------------------------
		#PRINT TO FILE
		#-------------------------
		print("PRINTING RESULTS TO FILE")
		print(time.time() - start)

		with open(positive_phages, 'w') as f:
			 for item in positive_hits_nucleotide:
				  f.write("%s\n" % item)

		with open(positive_orfs, 'w') as f:
			 for item in positive_hits_protein:
				  f.write("%s\n" % item)

def inputDatabaseExtract(wildcards):
	if wildcards.type=="ORFs":
		if wildcards.db=="RefSeqViral":
			return config["RefSeqViral_protein_db"]
		elif wildcards.db=="GenBankViral":
			return config["GenBankViral_protein_db"]
		elif wildcards.db=="METAVR":
			return config["METAVR_protein_db"]
	elif wildcards.type=="phages":
		if wildcards.db=="RefSeqViral":
			return config["RefSeqViral_db"]
		elif wildcards.db=="GenBankViral":
			return config["GenBankViral_db"]
		elif wildcards.db=="METAVR":
			return config["METAVR_reference_fasta"]

rule get_relatives_fasta:
	input:
		blast_relatives=(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_{type}_{db}.{sampling}.txt"),
		blast_database=inputDatabaseExtract
	output:
		blast_relatives=(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_{type}_{db}.{sampling}.fasta"),
	message:
		"Extracting relative sequences FASTA"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_relatives_fasta/db={db}__sampling={sampling}__type={type}.tsv"
	threads: 1
	shell:
		"""
		seqtk subseq {input.blast_database} {input.blast_relatives} > {output.blast_relatives}
		"""

rule makeblastdb_relatives_ORFs:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates_ORFs.{sampling}.fasta",
		blast_relatives=expand(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_ORFs_{db}.{{sampling}}.fasta", db=REFERENCE_DATABASES),
	output:
		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs.{sampling}.fasta",
		orfs_blast_db=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs.{sampling}.fasta.pdb",
		length_temp=temp(dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs_length_temp.{sampling}.txt"),
		names_temp=temp(dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs_names_temp.{sampling}.txt"),
		length=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs_lengths_and_names.{sampling}.txt",
		cummulative_length=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs_coding_lengths.{sampling}.txt",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/makeblastdb_relatives_ORFs/sampling={sampling}.tsv"
	message:
		"Creating ORFs relatives blast database"
	threads: 8
	shell:
		"""
		cat {input.filtered_positive_contigs} {input.blast_relatives} > {output.cat_isolates_relatives}
		cat {output.cat_isolates_relatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
			$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.length_temp}
		cat {output.length_temp} | cut -f1 | cut -f1 -d' ' | rev | cut -d_ -f2- | rev > {output.names_temp}
		paste {output.length_temp} {output.names_temp} > {output.length}
		cat {output.length} | cut -f3,2 | awk '{{arr[$2]+=$1}} END {{for (i in arr) {{print i,arr[i]}}}}' > {output.cummulative_length}
		makeblastdb -in {output.cat_isolates_relatives} -dbtype prot
		"""

rule blastall_relatives_ORFs:
	input:
		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs.{sampling}.fasta",
		orfs_blast_db=dirs_dict["ANNOTATION"]+ "/isolates_relatives_ORFs.{sampling}.fasta.pdb",
	output:
		isolates_relatives_blastall=(dirs_dict["ANNOTATION"] + "/isolates_relatives_blastall_ORFs.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/blastall_relatives_ORFs/sampling={sampling}.tsv"
	message:
		"Annotating contigs with BLAST"
	threads: 32
	shell:
		"""
		blastp -num_threads {threads} -db {input.cat_isolates_relatives} -query {input.cat_isolates_relatives} \
			-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs qcovhsp length pident evalue positive gaps" > {output.isolates_relatives_blastall}
		"""

rule cat_relatives_phages:
	input:
		filtered_positive_contigs=dirs_dict["vOUT_DIR"]+ "/filtered_isolates.{sampling}.fasta",
		blast_relatives=expand(dirs_dict["ANNOTATION"] + "/isolates_blast_relatives_phages_{db}.{{sampling}}.fasta", db=REFERENCE_DATABASES),
	output:
		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages.{sampling}.fasta",
		length_temp=temp(dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages_length_temp.{sampling}.txt"),
		names_temp=temp(dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages_names_temp.{sampling}.txt"),
		length=dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages_lengths_and_names.{sampling}.txt",
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/cat_relatives_phages/sampling={sampling}.tsv"
	message:
		"Creating ORFs relatives blast database"
	threads: 8
	shell:
		"""
		cat {input.filtered_positive_contigs} {input.blast_relatives} > {output.cat_isolates_relatives}
		cat {output.cat_isolates_relatives} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
			$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.length_temp}
		cat {output.length_temp} | cut -f1 | cut -f1 -d' ' | rev | cut -d_ -f2- | rev > {output.names_temp}
		paste {output.length_temp} {output.names_temp} > {output.length}
		"""

rule viridic_relatives_phages:
	input:
		cat_isolates_relatives=dirs_dict["ANNOTATION"]+ "/isolates_relatives_phages.{sampling}.fasta",
		viridic_singularity_folder=config["viridic_folder"]
	output:
		viridic_out=directory(dirs_dict["ANNOTATION"] + "/VIRIDIC_isolates_relatives_phages.{sampling}/"),
	message:
		"Calculating intergenomic similarities with VIRIDIC"
	params:
		steps="ALL",
		run_log="/dev/stdout",
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/viridic_relatives_phages/sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		mkdir -p $(dirname {output.viridic_out:q})
		# A legacy Conda Singularity can shadow the installed SIF-compatible runtime.
		viridic_runtime=""
		while IFS= read -r viridic_candidate; do
			viridic_version=$("$viridic_candidate" --version 2>&1) || continue
			if [[ "$viridic_version" =~ ([0-9]+)[.]([0-9]+) ]]; then
				if (( BASH_REMATCH[1] > 3 || (BASH_REMATCH[1] == 3 && BASH_REMATCH[2] >= 5) )); then
					viridic_runtime="$viridic_candidate"
					break
				fi
			fi
		done < <(type -a -P singularity)
		if [ -z "$viridic_runtime" ]; then
			echo "VIRIDIC requires Singularity >=3.5 on PATH; no compatible runtime was found." >&2
			exit 1
		fi
		export PATH="$(dirname "$viridic_runtime"):$PATH"
		# VIRIDIC's R future workers otherwise ignore ncor and use every server core.
		export SINGULARITYENV_R_FUTURE_AVAILABLECORES_SYSTEM={threads}
		export SINGULARITYENV_R_PARALLELLY_AVAILABLECORES_SYSTEM={threads}
		export SINGULARITYENV_MC_CORES={threads}
		export SINGULARITYENV_OMP_NUM_THREADS=1
		{input.viridic_singularity_folder:q}/viridic.bash projdir={output.viridic_out:q} \
			in={input.cat_isolates_relatives:q} ncor={threads} steps={params.steps} > {params.run_log:q} 2>&1
		test -s {output.viridic_out:q}/04_VIRIDIC_out/sim_MA_genCol.csv
		"""

rule genomad_host:
	input:
		host_fasta=dirs_dict["HOST_DIR"] + "/{host}.fasta",
		genomad_db=config["genomad_db"],
	output:
		genomad_outdir=directory(dirs_dict["HOST_DIR"] + "/{host}_geNomad"),
		positive_contigs=dirs_dict["HOST_DIR"] + "/prophages/{host}_prophages.fasta",
		evidence=dirs_dict["HOST_DIR"] + "/prophages/{host}_viral_evidence.tsv",
	message:
		"Identifying embedded prophages and whole-contig viral candidates in host {wildcards.host}"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/genomad_host/host={host}.tsv"
	threads: 16
	shell:
		"""
		genomad end-to-end --restart --cleanup -t {threads} {input.host_fasta:q} {output.genomad_outdir:q} {input.genomad_db:q}
		python - {output.genomad_outdir:q} {wildcards.host:q} {output.positive_contigs:q} {output.evidence:q} <<-'PYTHON'
		import csv
		import sys
		from pathlib import Path
		from Bio import SeqIO
		folder, host, fasta, evidence = sys.argv[1:]
		summary = Path(folder) / (host + "_summary")
		with (summary / (host + "_virus_summary.tsv")).open() as handle:
		    calls = {{row["seq_name"]: row for row in csv.DictReader(handle, delimiter="\\t")}}
		fields = ["host", "viral_id", "genomad_sequence", "parent_contig", "evidence_scope", "start", "end",
		          "length_bp", "virus_score", "taxonomy", "topology", "n_hallmarks", "genetic_code"]
		with open(fasta, "w") as output, open(evidence, "w") as table:
		    writer = csv.DictWriter(table, fieldnames=fields, delimiter="\\t", lineterminator="\\n")
		    writer.writeheader()
		    for record in SeqIO.parse(summary / (host + "_virus.fna"), "fasta"):
		        original = record.id
		        row = calls[original]
		        embedded = "|provirus_" in original
		        coordinates = row.get("coordinates", "").split("-") if embedded else ["1", str(len(record.seq))]
		        record.id = host + "_" + original.replace("|", "_")
		        record.description = ""
		        SeqIO.write(record, output, "fasta-2line")
		        writer.writerow(dict(host=host, viral_id=record.id, genomad_sequence=original,
		            parent_contig=original.split("|provirus_")[0],
		            evidence_scope="embedded_prophage" if embedded else "whole_contig_viral_candidate",
		            start=coordinates[0], end=coordinates[-1], length_bp=len(record.seq),
		            **{{key: row.get(key, "not reported") for key in fields[8:]}}))
		PYTHON
		"""

rule estimateGenomeCompletness_prophages:
	input:
		positive_contigs=dirs_dict["HOST_DIR"] + "/prophages/{host}_prophages.fasta",
		checkv_db=(config['checkv_db']),
	output:
		quality_summary=dirs_dict["HOST_DIR"] + "/prophages/{host}_checkV/quality_summary.tsv",
		completeness=dirs_dict["HOST_DIR"] + "/prophages/{host}_checkV/completeness.tsv",
		contamination=dirs_dict["HOST_DIR"] + "/prophages/{host}_checkV/contamination.tsv",
	params:
		checkv_outdir=dirs_dict["HOST_DIR"] + "/prophages/{host}_checkV",
		tmp=dirs_dict["HOST_DIR"] + "/prophages/{host}_checkV/tmp",
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/estimateGenomeCompletness_prophages/host={host}.tsv"
	threads: 32
	shell:
		"""
		rm -rf {params.checkv_outdir} || true
		if [ -s {input.positive_contigs} ]; then
			 							checkv contamination {input.positive_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
			 							checkv completeness {input.positive_contigs} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
			 							checkv complete_genomes {input.positive_contigs} {params.checkv_outdir}
			 							checkv quality_summary {input.positive_contigs} {params.checkv_outdir}
										rm -rf {params.tmp}
		else
			 							echo "The FASTA file {input.positive_contigs} is empty"
			 							mkdir -p {params.checkv_outdir}
			 							touch {output.quality_summary}
			 							touch {output.completeness}
			 							touch {output.contamination}
		fi
		"""

rule mask_prophages:
	input:
		host_fasta=dirs_dict["HOST_DIR"] + "/{host}.fasta",
		evidence=dirs_dict["HOST_DIR"] + "/prophages/{host}_viral_evidence.tsv",
	output:
		masked_prophages=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.fasta",
		mask_regions=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_mask_regions.bed",
	params:
		mask_additional_bases=0 if ISOLATES else 500,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mask_prophages/host={host}.tsv"
	threads: 1
	run:
		import csv
		from Bio import SeqIO

		with open(input.evidence) as handle:
			regions=list(csv.DictReader(handle, delimiter="\t"))
		with open(output.mask_regions, "w") as bed, open(output.masked_prophages, "w") as fasta:
			for record in SeqIO.parse(input.host_fasta, "fasta"):
				sequence=list(str(record.seq))
				whole_viral=False
				for row in regions:
					if row["parent_contig"] != record.id:
						continue
					start=max(0, int(row["start"]) - 1 - params.mask_additional_bases)
					end=min(len(sequence), int(row["end"]) + params.mask_additional_bases)
					bed.write(f"{record.id}\t{start}\t{end}\n")
					whole_viral=whole_viral or row["evidence_scope"] == "whole_contig_viral_candidate"
					sequence[start:end]=["N"] * (end - start)
				# A standalone viral contig is not host chromosome sequence.
				if not whole_viral:
					fasta.write(f">{record.id}\n{''.join(sequence)}\n")

rule buildBowtieDB_host:
	input:
		host_fasta = dirs_dict["HOST_DIR"] + "/{host}.fasta",
	output:
		contigs_bt2_1=temp(dirs_dict["HOST_DIR"]+ "/{host}.1.bt2"),
		contigs_bt2_2=temp(dirs_dict["HOST_DIR"]+ "/{host}.2.bt2"),
		contigs_bt2_3=temp(dirs_dict["HOST_DIR"]+ "/{host}.3.bt2"),
		contigs_bt2_4=temp(dirs_dict["HOST_DIR"]+ "/{host}.4.bt2"),
		contigs_bt2_rev_1=temp(dirs_dict["HOST_DIR"]+ "/{host}.rev.1.bt2"),
		contigs_bt2_rev_2=temp(dirs_dict["HOST_DIR"]+ "/{host}.rev.2.bt2"),
	params:
		prefix=dirs_dict["HOST_DIR"]+ "/{host}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_host/host={host}.tsv"
	threads: 8
	shell:
		"""
		if [ -s {input.host_fasta:q} ]; then
			bowtie2-build {input.host_fasta:q} {params.prefix:q} --threads {threads}
		else
			touch {output:q}
		fi
		"""

rule run_BLASTn_host:
	input:
		host_fasta=lambda wc: dirs_dict["HOST_DIR"] + ("/prophages/" + wc.host + "_prophages.fasta" if wc.host_scope else "/" + wc.host + ".fasta"),
		assembly_fasta=lambda wc: [dirs_dict["ASSEMBLY_DIR"] + "/" + wc.assembly_query.split("/")[1] + "_spades_filtered_scaffolds.tot.fasta"] if wc.assembly_query else ([ALL_ASSEMBLED_DIR + "/phage_isolates_contigs.tot.fasta"] if ISOLATES else expand(dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_spades_filtered_scaffolds.tot.fasta",sample=SAMPLES)),
	output:
		temp_fasta=temp(dirs_dict["ASSEMBLY_DIR"]+ "/{assembly_query}assembly_contigs_{host}{host_scope}.fasta"),
		blast_output=(dirs_dict["vOUT_DIR"] + "/{assembly_query}blastn_out_assembly_{host}{host_scope}.tot.csv"),
	params:
		outfmt="6 qseqid sseqid salltitles qstart qend qlen slen qcovs evalue length pident" + (" sstart send bitscore btop" if ISOLATES else ""),
		reference=lambda wc, input: ("-subject " if wc.assembly_query else "-db ") + shlex.quote(input.host_fasta),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/run_BLASTn_host/{assembly_query}host={host}__scope={host_scope}.tsv"
	wildcard_constraints:
		assembly_query="|ISOLATES/[^/]+/",
		host="|".join(re.escape(host) for host in HOSTS) or "(?!)",
		host_scope="|_viral_regions",
	message:
		"Annotating contigs with BLAST"
	threads: 8
	shell:
		"""
		cat {input.assembly_fasta:q} > {output.temp_fasta:q}
		if [ -s {input.host_fasta:q} ] && [ -s {output.temp_fasta:q} ]; then
			if [ -z {wildcards.assembly_query:q} ]; then
				makeblastdb -in {input.host_fasta:q} -dbtype nucl
			fi
			blastn -num_threads {threads} {params.reference} -query {output.temp_fasta:q}\
				-outfmt {params.outfmt:q} > {output.blast_output:q}
		else
			printf '' > {output.blast_output:q}
		fi
		"""

# Both references are mapped sequentially in one job per sample/host.
rule map_to_host:
	input:
		host_fasta=dirs_dict["HOST_DIR"] + "/{host}.fasta",
		masked_fasta=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.fasta",
		contigs_bt2_1=(dirs_dict["HOST_DIR"] + "/{host}.1.bt2"),
		contigs_bt2_2=(dirs_dict["HOST_DIR"] + "/{host}.2.bt2"),
		contigs_bt2_3=(dirs_dict["HOST_DIR"] + "/{host}.3.bt2"),
		contigs_bt2_4=(dirs_dict["HOST_DIR"] + "/{host}.4.bt2"),
		contigs_bt2_rev_1=(dirs_dict["HOST_DIR"] + "/{host}.rev.1.bt2"),
		contigs_bt2_rev_2=(dirs_dict["HOST_DIR"] + "/{host}.rev.2.bt2"),
		masked_bt2_1=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.1.bt2",
		masked_bt2_2=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.2.bt2",
		masked_bt2_3=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.3.bt2",
		masked_bt2_4=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.4.bt2",
		masked_bt2_rev_1=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.rev.1.bt2",
		masked_bt2_rev_2=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.rev.2.bt2",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_flagstats_{sample}_vs_{host}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_flagstats_filtered_{sample}_vs_{host}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_covstats.txt",
		covstats_filtered=dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_filtered_{sample}_vs_{host}_covstats.txt",
		basecov=(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_{sample}_vs_{host}_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		basecov_filtered=(dirs_dict["MAPPING_DIR"]+ "/HOST/bowtie2_filtered_{sample}_vs_{host}_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		masked_sam=temp(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_{sample}_vs_{host}_masked_prophages.sam"),
		masked_bam=temp(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_{sample}_vs_{host}_masked_prophages.bam"),
		masked_sorted_bam=temp(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_{sample}_vs_{host}_masked_prophages_sorted.bam"),
		masked_sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_{sample}_vs_{host}_masked_prophages_sorted.bam.bai"),
		masked_filtered_bam=temp(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_{sample}_vs_{host}_masked_prophages_filtered.bam"),
		masked_flagstats=dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_flagstats_{sample}_vs_{host}_masked_prophages.txt",
		masked_flagstats_filtered=dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_flagstats_filtered_{sample}_vs_{host}_masked_prophages.txt",
		masked_covstats=dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_{sample}_vs_{host}_masked_prophages_covstats.txt",
		masked_covstats_filtered=dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_filtered_{sample}_vs_{host}_masked_prophages_covstats.txt",
		masked_basecov=(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_{sample}_vs_{host}_masked_prophages_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		masked_basecov_filtered=(dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_filtered_{sample}_vs_{host}_masked_prophages_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unpaired=mapping_orphan_option,
		prefix=dirs_dict["HOST_DIR"]+ "/{host}",
		masked_prefix=dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages",
	message:
		"Mapping reads to the host and prophage-masked host"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/map_to_host/host={host}__sample={sample}.tsv"
	threads: 8
	shell:
		"""
		# Unmasked host
		if [ -s {input.host_fasta:q} ]; then
			bowtie2 -X {params.max_insert_bp} -x {params.prefix:q} -1 {input.forward_paired:q} -2 {input.reverse_paired:q} {params.unpaired} -S {output.sam:q} --threads {threads} --no-unal --all --very-sensitive
		else
			printf '@HD\\tVN:1.6\\tSO:coordinate\\n' > {output.sam:q}
		fi
		samtools view -@ {threads} -bS {output.sam:q} > {output.bam:q}
		samtools sort -@ {threads} {output.bam:q} -o {output.sorted_bam:q}
		samtools index {output.sorted_bam:q}
		
		if [ -s {input.host_fasta:q} ]; then
			coverm filter -b {output.sorted_bam:q} -o {output.filtered_bam:q} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		else
			cp {output.sorted_bam:q} {output.filtered_bam:q}
		fi
		
		samtools flagstat {output.sorted_bam:q} > {output.flagstats:q}
		samtools flagstat {output.filtered_bam:q} > {output.flagstats_filtered:q}

		if [ -s {input.host_fasta:q} ]; then
			coverm contig -b {output.sorted_bam:q} -m mean length covered_bases count variance trimmed_mean rpkm -o {output.covstats:q}
			coverm contig -b {output.filtered_bam:q} -m mean length covered_bases count variance trimmed_mean rpkm -o {output.covstats_filtered:q}
			for basecov in {output.basecov:q}; do
				bedtools genomecov -dz -ibam {output.sorted_bam:q} | gzip -c > "$basecov"
			done
			for basecov in {output.basecov_filtered:q}; do
				bedtools genomecov -dz -ibam {output.filtered_bam:q} | gzip -c > "$basecov"
			done
		else
			printf 'Contig\\tMean\\tLength\\tCovered Bases\\tRead Count\\tVariance\\tTrimmed Mean\\tRPKM\\n' > {output.covstats:q}
			cp {output.covstats:q} {output.covstats_filtered:q}
			for basecov in {output.basecov:q}; do
				printf '' | gzip -c > "$basecov"
			done
			for basecov in {output.basecov_filtered:q}; do
				printf '' | gzip -c > "$basecov"
			done
		fi

		# Prophage-masked host
		if [ -s {input.masked_fasta:q} ]; then
			bowtie2 -X {params.max_insert_bp} -x {params.masked_prefix:q} -1 {input.forward_paired:q} -2 {input.reverse_paired:q} {params.unpaired} -S {output.masked_sam:q} --threads {threads} --no-unal --all --very-sensitive
		else
			printf '@HD\\tVN:1.6\\tSO:coordinate\\n' > {output.masked_sam:q}
		fi
		samtools view -@ {threads} -bS {output.masked_sam:q} > {output.masked_bam:q}
		samtools sort -@ {threads} {output.masked_bam:q} -o {output.masked_sorted_bam:q}
		samtools index {output.masked_sorted_bam:q}
		samtools flagstat {output.masked_sorted_bam:q} > {output.masked_flagstats:q}
		if [ -s {input.masked_fasta:q} ]; then
			coverm filter -b {output.masked_sorted_bam:q} -o {output.masked_filtered_bam:q} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		else
			cp {output.masked_sorted_bam:q} {output.masked_filtered_bam:q}
		fi
		samtools flagstat {output.masked_filtered_bam:q} > {output.masked_flagstats_filtered:q}

		if [ -s {input.masked_fasta:q} ]; then
			coverm contig -b {output.masked_sorted_bam:q} -m mean length covered_bases count variance trimmed_mean rpkm -o {output.masked_covstats:q}
			coverm contig -b {output.masked_filtered_bam:q} -m mean length covered_bases count variance trimmed_mean rpkm -o {output.masked_covstats_filtered:q}
			for basecov in {output.masked_basecov:q}; do
				bedtools genomecov -dz -ibam {output.masked_sorted_bam:q} | gzip -c > "$basecov"
			done
			for basecov in {output.masked_basecov_filtered:q}; do
				bedtools genomecov -dz -ibam {output.masked_filtered_bam:q} | gzip -c > "$basecov"
			done
		else
			printf 'Contig\\tMean\\tLength\\tCovered Bases\\tRead Count\\tVariance\\tTrimmed Mean\\tRPKM\\n' > {output.masked_covstats:q}
			cp {output.masked_covstats:q} {output.masked_covstats_filtered:q}
			for basecov in {output.masked_basecov:q}; do
				printf '' | gzip -c > "$basecov"
			done
			for basecov in {output.masked_basecov_filtered:q}; do
				printf '' | gzip -c > "$basecov"
			done
		fi
		"""

rule clustering_isolates:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/combined_positive_viral_contigs.tot.fasta",
	output:
		clusters=dirs_dict["vOUT_DIR"] + "/clustering_combined_positive_viral_contigs_{identity}-{coverage}.clstr",
		blastout=temp(dirs_dict["vOUT_DIR"] + "/clustering_combined_positive_viral_contigs_{identity}-{coverage}_blastout.csv"),
		aniout=temp(dirs_dict["vOUT_DIR"] + "/clustering_combined_positive_viral_contigs_{identity}-{coverage}.csv"),
	message:
		"Creating vOUTs with CheckV aniclust"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/clustering_isolates/coverage={coverage}__identity={identity}.tsv"
	#  benchmark:
	#		dirs_dict['BENCHMARKS']+ "/vOUTclustering/{sequence}.tsv",
	threads: 144
	wildcard_constraints:
		identity="95",
		coverage="85",
	shell:
		"""
		makeblastdb -in {input.fasta} -dbtype nucl -out {input.fasta}
		blastn -query {input.fasta} -db {input.fasta} -outfmt '6 std qlen slen' \
				-max_target_seqs 10000000 -out {output.blastout} -num_threads {threads}
		python scripts/anicalc_checkv.py  -i {output.blastout} -o {output.aniout}
		aniclust --fna {input.fasta} --ani {output.aniout} --out {output.clusters} --min_ani 95 --min_tcov 85 --min_qcov 0
		"""

rule estimateBacterialGenomeCompletness_host:
	input:
		host_fasta = dirs_dict["HOST_DIR"] + "/{host}.fasta",
		checkm_db=(config['checkm_db']),
	output:
		checkMoutdir_temp=temp(directory(dirs_dict["HOST_DIR"] + "/{host}_checkM_temp")),
		checkMoutdir=temp(directory(dirs_dict["HOST_DIR"] + "/{host}_checkM")),
	params:
		checkv_db=dirs_dict["HOST_DIR"] + "/{host}_checkV",
	log:
		checkMoutdir=temp(dirs_dict["HOST_DIR"] + "/{host}_checkM.log"),
	message:
		"Estimating genome completeness with CheckM "
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/estimateBacterialGenomeCompletness_host/host={host}.tsv"
	threads: 4
	shell:
		"""
		mkdir -p {output.checkMoutdir_temp}
		cp {input.host_fasta} {output.checkMoutdir_temp}
		cd {output.checkMoutdir_temp}
		checkm lineage_wf -t {threads} -x fasta {output.checkMoutdir_temp} {output.checkMoutdir} 1> {log}
		"""

rule combine_logs_to_csv_hosts:
	input:
		logs = expand((dirs_dict["HOST_DIR"] + "/{host}_checkM.log"), host=HOSTS)
	output:
		csv = dirs_dict["HOST_DIR"] + "/checkM_summary.csv"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/combine_logs_to_csv_hosts/tot.tsv"
	run:
		import pandas as pd
		from io import StringIO
		import glob
		dfs = []
		for log_file in input.logs:
			with open(log_file, "r") as f:
					lines = f.readlines()

			# Find the table start and end
			start_idx = None
			end_idx = None
			for i, line in enumerate(lines):
					if line.strip().startswith("Bin Id"):
						start_idx = i
					elif start_idx and line.strip().startswith("[") and "INFO" in line:
						end_idx = i
						break

			if start_idx is not None:
					# Extract the table lines
					table_lines = lines[start_idx:end_idx]
					# Remove separator lines (-----)
					table_lines = [l for l in table_lines if not l.strip().startswith('---')]
					# Join lines and read with pandas
					table_str = "".join(table_lines)
					df = pd.read_csv(StringIO(table_str), sep=r'\s{2,}', engine='python')
					dfs.append(df)

		# Combine all tables
		if len(dfs) == 0:
			# Create an empty CSV with no rows
			final_df = pd.DataFrame()
		elif len(dfs) == 1:
			final_df = dfs[0]
		else:
			final_df = pd.concat(dfs, ignore_index=True)

		final_df.to_csv(output.csv, index=False)
