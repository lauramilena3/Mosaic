# ruleorder: vOUTclustering_references>vOUTclustering
# ruleorder: filter_vOTUs_references>filter_vOTUs

def input_vOTU_clustering(wildcards):
	input_list=[]
	if (not NANOPORE_ONLY) & (not ISOLATES):
		input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta",sample=SAMPLES))
	if NANOPORE:
		input_list.extend(expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_"+ LONG_ASSEMBLER + "_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.fasta", sample=NANOPORE_SAMPLES))
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["VIRAL_DIR"]+ "/ALL_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta")
	if SUBASSEMBLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_positive_geNomad.{{sampling}}.fasta", sample=SAMPLES, subsample=subsample_test))
	if len(config['additional_reference_contigs'])>0:
		input_list.append(config['additional_reference_contigs'])
	if ISOLATES:
		input_list.extend(expand(dirs_dict["HOST_DIR"] + "/prophages/{host}_prophages.fasta", host=HOSTS))
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_spades_filtered_scaffolds.tot.fasta",sample=SAMPLES))
	if RNA_MODE and wildcards.sampling == "tot":
		input_list.extend(expand(RNA_DIR + "/{sample}_virsorter2_RNA_viral_positive.fasta", sample=SAMPLES))
	return input_list

# if len(config['additional_reference_contigs'])==0:

rule derreplicate_assembly:
	input:
		positive_contigs=input_vOTU_clustering
	output:
		combined_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
		derreplicated_clusters=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_cluster.tsv",
		derreplicated_tmp=directory(dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_tmp"),
	params:
		rep_name="combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated",
		rep_name_full=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_rep_seq.fasta",
		rep_temp="combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_tmp",
		dir_votu=dirs_dict["vOUT_DIR"],
		rna_enabled=lambda wc: RNA_MODE and wc.sampling == "tot",
	message:
		"Derreplicating assembled contigs with mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/derreplicate_assembly/sampling={sampling}.tsv"
	threads: 16
	shell:
		"""
		cat {input.positive_contigs} > {output.combined_positive_contigs}
		cd {params.dir_votu}
		mmseqs easy-cluster --threads {threads} --createdb-mode 1 --min-seq-id 1 -c 1 --cov-mode 1 {output.combined_positive_contigs} {params.rep_name} {params.rep_temp} 
		mv {params.rep_name_full} {output.derreplicated_positive_contigs}
		rm -f -- {params.rep_name:q}_all_seqs.fasta
		"""

rule combine_all_assembled_contigs:
	input:
		dna=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.tot.fasta", sample=SAMPLES),
		rna=lambda wc: expand(RNA_DIR + "/{sample}_{assembler}.fasta", sample=SAMPLES, assembler=RNA_ASSEMBLERS) if RNA_MODE and wc.catalogue == "all_assembled" else [],
	output:
		fasta=ALL_ASSEMBLED_DIR + "/{catalogue}_contigs.tot.fasta",
		provenance=ALL_ASSEMBLED_DIR + "/{catalogue}_contigs_provenance.tot.tsv",
	wildcard_constraints:
		catalogue="all_assembled|phage_isolates" if RNA_MODE else "all_assembled",
	message:
		"Combining retained assemblies into the {wildcards.catalogue} catalogue"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/combine_all_assembled_contigs/catalogue={catalogue}.tsv"
	threads: 1
	run:
		import csv

		def records(path):
			name, chunks = None, []
			with open(path) as handle:
				for line in handle:
					if line.startswith(">"):
						if name is not None:
							yield name, "".join(chunks)
						name, chunks = line[1:].split()[0], []
					else:
						chunks.append(line.strip())
			if name is not None:
				yield name, "".join(chunks)

		os.makedirs(ALL_ASSEMBLED_DIR, exist_ok=True)
		with open(output.fasta, "w") as fasta, open(output.provenance, "w") as handle:
			writer=csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["contig_id", "sample", "assembler", "length_bp"])
			for path in list(input.dna) + list(input.rna):
				if path.endswith("_spades_filtered_scaffolds.tot.fasta"):
					sample=os.path.basename(path).removesuffix("_spades_filtered_scaffolds.tot.fasta")
				else:
					sample=os.path.basename(path).removesuffix(".fasta").rsplit("_", 1)[0]
				for contig_id, sequence in records(path):
					assembler=contig_id.rsplit("_", 4)[1]
					fasta.write(f">{contig_id}\n{sequence}\n")
					writer.writerow([contig_id, sample, assembler, len(sequence)])

rule derreplicate_all_assembled_contigs:
	input:
		fasta=ALL_ASSEMBLED_DIR + "/{catalogue}_contigs.tot.fasta",
	output:
		fasta=ALL_ASSEMBLED_DIR + "/{catalogue}_contigs_derreplicated_rep_seq.tot.fasta",
		clusters=ALL_ASSEMBLED_DIR + "/{catalogue}_contigs_derreplicated_cluster.tot.tsv",
		tmp=directory(ALL_ASSEMBLED_DIR + "/{catalogue}_contigs_derreplicated_tmp"),
	params:
		prefix=ALL_ASSEMBLED_DIR + "/{catalogue}_contigs_derreplicated",
	wildcard_constraints:
		catalogue="all_assembled|phage_isolates" if RNA_MODE else "all_assembled",
	message:
		"Derreplicating all assembled contigs with mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/derreplicate_all_assembled_contigs/catalogue={catalogue}.tsv"
	threads: 16
	shell:
		"""
		if [ -s {input.fasta:q} ]; then
			mmseqs easy-cluster --threads {threads} --createdb-mode 1 --min-seq-id 1 -c 1 --cov-mode 1 \
				{input.fasta:q} {params.prefix:q} {output.tmp:q}
			mv {params.prefix:q}_rep_seq.fasta {output.fasta:q}
			mv {params.prefix:q}_cluster.tsv {output.clusters:q}
		else
			mkdir -p {output.tmp:q}
			printf '' > {output.fasta:q}
			printf '' > {output.clusters:q}
		fi
		rm -f -- {params.prefix:q}_all_seqs.fasta
		"""

if not RNA_MODE:
	rule reuse_dna_contigs_for_isolates:
		input:
			fasta=ALL_ASSEMBLED_DIR + "/all_assembled_contigs.tot.fasta",
			provenance=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_provenance.tot.tsv",
		output:
			fasta=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs.tot.fasta",
			provenance=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_provenance.tot.tsv",
		params:
			fasta=lambda wc, input: os.path.basename(input.fasta),
			provenance=lambda wc, input: os.path.basename(input.provenance),
		message:
			"Reusing the DNA catalogue and provenance for phage isolates"
		benchmark:
			dirs_dict["BENCHMARKS"] + "/reuse_dna_contigs_for_isolates/tot.tsv"
		threads: 1
		shell:
			"""
			ln -sfn -- {params.fasta:q} {output.fasta:q}
			ln -sfn -- {params.provenance:q} {output.provenance:q}
			"""

	rule reuse_dna_dereplication_for_isolates:
		input:
			fasta=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_derreplicated_rep_seq.tot.fasta",
			clusters=ALL_ASSEMBLED_DIR + "/all_assembled_contigs_derreplicated_cluster.tot.tsv",
		output:
			fasta=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_rep_seq.tot.fasta",
			clusters=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_cluster.tot.tsv",
		params:
			fasta=lambda wc, input: os.path.basename(input.fasta),
			clusters=lambda wc, input: os.path.basename(input.clusters),
		message:
			"Reusing exact DNA dereplication for phage isolates"
		benchmark:
			dirs_dict["BENCHMARKS"] + "/reuse_dna_dereplication_for_isolates/tot.tsv"
		threads: 1
		shell:
			"""
			ln -sfn -- {params.fasta:q} {output.fasta:q}
			ln -sfn -- {params.clusters:q} {output.clusters:q}
			"""

rule vOUTclustering:
	input:
		fasta="{basedir}/{sequence}.fasta",
	output:
		clusters="{basedir}/{sequence}_95-85.clstr",
		blastout="{basedir}/{sequence}-blastout.csv",
		aniout="{basedir}/{sequence}-aniout.csv",
	message:
		"Creating vOUTs with CheckV aniclust"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/vOUTclustering/basedir={basedir}__sequence={sequence}.tsv"
	#  benchmark:
	#		dirs_dict['BENCHMARKS']+ "/vOUTclustering/{sequence}.tsv",
	threads: 144
	wildcard_constraints:
		sequence="[^/]+"  # The 'sequence' wildcard cannot contain a slash
	shell:
		"""
		if [ -s {input.fasta:q} ]; then
			makeblastdb -in {input.fasta:q} -dbtype nucl -out {input.fasta:q}
			blastn -query {input.fasta:q} -db {input.fasta:q} -outfmt '6 std qlen slen' \
				-max_target_seqs 10000000 -out {output.blastout:q} -num_threads {threads}
			if [ -s {output.blastout:q} ]; then
				python scripts/anicalc_checkv.py -i {output.blastout:q} -o {output.aniout:q}
			else
				printf 'qname\ttname\tnum_alns\tpid\tqcov\ttcov\n' > {output.aniout:q}
			fi
			aniclust --fna {input.fasta:q} --ani {output.aniout:q} --out {output.clusters:q} --min_ani 95 --min_tcov 85 --min_qcov 0
		else
			printf '' > {output.blastout:q}
			printf 'qname\ttname\tnum_alns\tpid\tqcov\ttcov\n' > {output.aniout:q}
			printf '' > {output.clusters:q}
		fi
		"""

rule collect_all_assembled_cluster_representatives:
	input:
		fasta=ALL_ASSEMBLED_SELECTED_PREFIX + ".fasta",
		clusters=ALL_ASSEMBLED_SELECTED_PREFIX + "_95-85.clstr",
		ranking=ALL_ASSEMBLED_DIR + "/all_assembled_selected_contigs_tot.tsv",
		membership=ALL_ASSEMBLED_DIR + "/all_assembled_selected_contigs_membership_tot.tsv",
	output:
		fasta=ALL_ASSEMBLED_TOP_PREFIX + ".fasta",
		ranking=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_contigs_tot.tsv",
		membership=ALL_ASSEMBLED_MAPPING_DIR + "/all_assembled_top_contigs_membership_tot.tsv",
	message:
		"Collecting every selected cluster representative and its original full-catalogue mapping measurements"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/collect_all_assembled_cluster_representatives/tot.tsv"
	threads: 1
	resources:
		mem_mb=16000,
	run:
		import pandas as pd
		from Bio import SeqIO

		clusters=pd.read_csv(input.clusters, sep="\t", header=None, names=["representative_id", "dereplicated_representative_id"])
		clusters["dereplicated_representative_id"]=clusters["dereplicated_representative_id"].str.split(",")
		clusters=clusters.explode("dereplicated_representative_id").drop_duplicates()
		ranking=pd.read_csv(input.ranking, sep="\t").set_index("contig_id", drop=False)
		top=ranking.loc[ranking.index.isin(clusters["representative_id"])].copy()
		top=top.rename(columns={"rank": "selection_rank"})
		top.insert(0, "rank", range(1, len(top) + 1))

		# Expand ANI membership through the existing MMseqs membership, without old-ID lookups.
		members=pd.read_csv(input.membership, sep="\t").rename(
			columns={"representative_id": "dereplicated_representative_id", "rank": "selection_rank"})
		membership=clusters.merge(members, on="dereplicated_representative_id", how="left")
		membership["rank"]=membership["representative_id"].map(top["rank"])
		membership=membership.sort_values(["rank", "member_id"], kind="stable")
		top["cluster_size"]=membership.groupby("representative_id")["member_id"].nunique().reindex(top.index).astype(int)
		top["selected_contigs_in_cluster"]=clusters.groupby("representative_id").size().reindex(top.index).astype(int)
		top["abundance_scope"]="representative_contig_full_catalogue_mapping"
		membership.to_csv(output.membership, sep="\t", index=False)
		top.to_csv(output.ranking, sep="\t", index=False)

		# Preserve the longest-first centroids chosen by the existing ANI clustering helper.
		with open(input.fasta) as handle:
			sequences={record.id: record for record in SeqIO.parse(handle, "fasta") if record.id in top.index}
		with open(output.fasta, "w") as handle:
			SeqIO.write((sequences[name] for name in top.index), handle, "fasta")

def input_getHighQuality(wildcards):
	input_list=[]
	if (not NANOPORE_ONLY) & (not ISOLATES):
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{{sampling}}/quality_summary.tsv",sample=SAMPLES)),
	if NANOPORE:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{{sampling}}/quality_summary.tsv", sample=NANOPORE_SAMPLES)),
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["vOUT_DIR"] + "/ALL_checkV_{sampling}/quality_summary.tsv"),
	if SUBASSEMBLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_geNomad_checkV_{{sampling}}/quality_summary.tsv", sample=SAMPLES, subsample=subsample_test)),
	if len(config['additional_reference_contigs'])>0:
		input_list.append(dirs_dict["vOUT_DIR"] + "/user_reference_contigs_checkV/quality_summary.tsv"),
	if ISOLATES:
		input_list.extend(expand(dirs_dict["HOST_DIR"] + "/prophages/{host}_checkV/quality_summary.tsv", host=HOSTS))
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/checkV_isolates_{sample}_tot/quality_summary.tsv",sample=SAMPLES)),
	if RNA_MODE and wildcards.sampling == "tot":
		input_list.extend(expand(RNA_DIR + "/{sample}/checkv/quality_summary.tsv", sample=SAMPLES))
	return input_list

rule getHighQuality:
	input:
		input_getHighQuality,
	output:
		quality_summary_concat=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		high_qualty_list=dirs_dict["vOUT_DIR"] + "/checkV_high_quality.{sampling}.txt",
	params:
		rna_enabled=lambda wc: RNA_MODE and wc.sampling == "tot",
	message:
		"Getting list high-quality vOTUs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/getHighQuality/sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		awk 'FNR>1' {input} > {output.quality_summary_concat}
		awk -F '\t' '/High-quality/ {{print $1}}' {output.quality_summary_concat} > {output.high_qualty_list}
		"""

checkpoint getHighQuality_clusters_fasta:
	input:
		new_clusters = dirs_dict["vOUT_DIR"] + "/new_references_clusters.{sampling}.csv",
		high_quality_list = dirs_dict["vOUT_DIR"] + "/checkV_high_quality.{sampling}.txt",
		combined_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		complete_clusters = dirs_dict["vOUT_DIR"] + "/new_references_complete_clusters.{sampling}.csv",
		fasta_dir = directory(dirs_dict["vOUT_DIR"] + "/high_quality_fastas.{sampling}")
	message:
		"Filtering clusters and extracting high-quality vOTUs FASTA sequences"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/getHighQuality_clusters_fasta/sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		mkdir -p {output.fasta_dir}
		awk 'NR==FNR {{a[$1]; next}} $1 in a' {input.high_quality_list} {input.new_clusters} > {output.complete_clusters}
		awk '{{print $2 >> "{output.fasta_dir}/" $1 ".list"}}' {output.complete_clusters}

		for listfile in {output.fasta_dir}/*.list; do
			rep=$(basename "$listfile" .list)
			n_lines=$(wc -l < "$listfile" | tr -d '[:space:]')  # make sure to strip spaces
			if [ "$n_lines" -gt 1 ]; then
				seqtk subseq {input.combined_positive_contigs} "$listfile" > {output.fasta_dir}/"$rep".fasta
			fi
			rm "$listfile"
		done
		"""

rule combine_with_taxmyphage:
	input:
		ref_fasta = "{contigs}.fasta",
		results_dir=(dirs_dict["ANNOTATION"] + "/taxmyphage_combined_positive_viral_contigs.tot"),
	output:
		blastout = temp("{contigs}_references.blastout"),
		aniout = temp("{contigs}_references.aniout"),
		clusters = temp("{contigs}_references.clusters"),
		cluster_rep = temp("{contigs}_references_cluster_rep.txt"),
		tax_fasta_rep = temp("{contigs}_references_cluster_rep.fasta"),
		combined = "{contigs}_with_references.fasta"
	params:
		tax_fasta = lambda wildcards: (dirs_dict["ANNOTATION"] + f"/taxmyphage_combined_positive_viral_contigs.tot/Results_per_genome/{os.path.basename(wildcards.contigs)}/known_taxa.fa")
	# wildcard_constraints:
	# 	contigs="(?!.*with_reference)[^/]+"
	message:
			"Combining {input.ref_fasta} with taxmyphage result: {params.tax_fasta}"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/combine_with_taxmyphage/contigs={contigs}.tsv"
	threads: 4
	shell:
		"""
		if [ -s {params.tax_fasta} ]; then
			makeblastdb -in {params.tax_fasta} -dbtype nucl -out {params.tax_fasta}
			blastn -query {params.tax_fasta} -db {params.tax_fasta} -outfmt '6 std qlen slen' \
					-max_target_seqs 10000000 -out {output.blastout} -num_threads {threads}
			python scripts/anicalc_checkv.py  -i {output.blastout} -o {output.aniout}
			aniclust --fna {params.tax_fasta} --ani {output.aniout} --out {output.clusters} --min_ani 95 --min_tcov 85 --min_qcov 0
			cut -f1 {output.clusters} | sort | uniq > {output.cluster_rep}
			seqtk subseq {params.tax_fasta} {output.cluster_rep} > {output.tax_fasta_rep}
		else
			touch {output.blastout} {output.aniout} {output.clusters} {output.cluster_rep} {output.tax_fasta_rep}> {output.combined}
		fi
		cat {input.ref_fasta} {output.tax_fasta_rep}> {output.combined}
		"""


rule select_vOTU_representative:
	input:
		merged_summary=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		cluster_file=dirs_dict["vOUT_DIR"] + "/combined_"+ VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}_95-85.clstr",
		derreplicated_clusters=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + ".{sampling}_derreplicated_cluster.tsv",
	output:
		representatives=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list.{sampling}.csv",
		checkv_categories=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list_checkv_per_category.{sampling}.csv",
		new_clusters=dirs_dict["vOUT_DIR"]+ "/new_references_clusters.{sampling}.csv"
	params:
		samples=SAMPLES,
		contig_dir=dirs_dict["ASSEMBLY_DIR"],
		viral_dir=dirs_dict['VIRAL_DIR'],
		subassembly=SUBASSEMBLY,
		cross_assembly=CROSS_ASSEMBLY,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/select_vOTU_representative/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/05_vOTU_representative.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/05_vOTU_representative.py.ipynb"

rule vOUTclustering_get_new_references:
	input:
		combined_positive_contigs=lambda wildcards: ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_rep_seq.tot.fasta" if wildcards.reference_catalogue == "phage_isolates_cluster_representatives" else dirs_dict["vOUT_DIR"] + "/combined_" + VIRAL_CONTIGS_BASE + "." + wildcards.sampling + ".fasta",
		representative_list=lambda wildcards: ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_rep_seq.tot_95-85.clstr" if wildcards.reference_catalogue == "phage_isolates_cluster_representatives" else dirs_dict["vOUT_DIR"] + "/vOTU_clustering_rep_list." + wildcards.sampling + ".csv",
	output:
		representatives="{basedir}/{reference_catalogue}.{sampling}.fasta",
		representative_lengths="{basedir}/{reference_catalogue}_lengths.{sampling}.txt",
	wildcard_constraints:
		basedir="(?:" + re.escape(dirs_dict["vOUT_DIR"]) + "|" + re.escape(ALL_ASSEMBLED_DIR) + ")",
		reference_catalogue="(?:" + re.escape(REPRESENTATIVE_CONTIGS_BASE) + "|phage_isolates_cluster_representatives)",
	params:
		cluster_centroids=lambda wildcards: wildcards.reference_catalogue == "phage_isolates_cluster_representatives",
	message:
		"Selecting new representatives with seqtk"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/vOUTclustering_get_new_references/basedir={basedir}__catalogue={reference_catalogue}__sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		if [ "{params.cluster_centroids}" = "True" ]; then
			seqtk subseq {input.combined_positive_contigs:q} <(cut -f1 {input.representative_list:q}) > {output.representatives:q}
			seqtk comp {output.representatives:q} | cut -f1,2 > {output.representative_lengths:q}
		else
			seqtk subseq {input.combined_positive_contigs:q} {input.representative_list:q} > {output.representatives:q}
			cat {output.representatives:q} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
				$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.representative_lengths:q}
		fi
		"""
		
rule get_list_filtered_vOTUs:
	input:
		df_counts_paired=dirs_dict["PLOTS_DIR"] + "/01_qc_read_counts_paired.{sampling}.csv",
		vOTUs_prefiltered=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		merged_summary=dirs_dict["vOUT_DIR"] + "/checkV_merged_quality_summary.{sampling}.txt",
		vibrant_circular=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_circular.{sampling}.csv",
		vibrant_positive=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_list.{sampling}.csv",
		vibrant_quality=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_quality.{sampling}.csv",
		vibrant_summary=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_summary_results.{sampling}.csv",
		virsorter_table=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/final-viral-score.tsv",
		virsorter_positive_list=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/positive_VS_list_{sampling}.txt",	
		genomad_virus_summary=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_virus_summary.tsv",
		genomad_plasmid_summary=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_plasmid_summary.tsv",
		genomad_viral_fasta=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/formatted_viral_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",										
		genomad_viral_fasta_conservative=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/formatted_viral_" + REPRESENTATIVE_CONTIGS_BASE + "_conservative.{sampling}.fasta",
		map_unfiltered=expand(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_unfiltered_contigs.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT),
		covstats_unfiltered=expand(dirs_dict["MAPPING_DIR"] + "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs.{sampling}_covstats.txt", sample=SAMPLES, sampling=SAMPLING_TYPE_TOT),

	output:
		summary=dirs_dict["vOUT_DIR"] + "/vOTU_clustering_summary.{sampling}.csv",
		filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
	params:
		samples=SAMPLES,
		contig_dir=dirs_dict["ASSEMBLY_DIR"],
		viral_dir=dirs_dict['VIRAL_DIR'],
		mapping_dir=dirs_dict['MAPPING_DIR'],
		subassembly=SUBASSEMBLY,
		cross_assembly=CROSS_ASSEMBLY,
		min_votu_len=config['min_votu_length'],
		key_samples=SAMPLES_key,
		rna_enabled=lambda wc: RNA_MODE and wc.sampling == "tot",
		rna_min_length=int(config.get("rna_min_contig_length", 500)),
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_list_filtered_vOTUs/sampling={sampling}.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/05_vOTU_filtering.{sampling}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/05_vOTU_filtering.py.ipynb"

rule filter_vOTUs:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
	output:
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		min_votu_len=config['min_votu_length']
	message:
		"Filtering vOTUs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/filter_vOTUs/sampling={sampling}.tsv"
	threads: 2
	shell:
		"""
		seqtk subseq {input.representatives} {input.filtered_list} > {output.filtered_representatives}
		"""

rule clustered_with_filter_vOTUs:
	input:
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
		new_clusters=dirs_dict["vOUT_DIR"]+ "/new_references_clusters.{sampling}.csv",
		filtered_list=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
	output:
		cluster_filtered_representatives_list=dirs_dict["vOUT_DIR"]+ "/viral_contigs_clustered_with_filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_list.{sampling}.txt",
		cluster_filtered_representatives_fasta=dirs_dict["vOUT_DIR"]+ "/viral_contigs_clustered_with_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		min_votu_len=config['min_votu_length']
	message:
		"Filtering vOTUs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/clustered_with_filter_vOTUs/sampling={sampling}.tsv"
	threads: 2
	shell:
		"""
		grep -f {input.filtered_list} {input.new_clusters} | cut -f2 > {output.cluster_filtered_representatives_list}
		seqtk subseq {input.derreplicated_positive_contigs} {output.cluster_filtered_representatives_list} > {output.cluster_filtered_representatives_fasta}
		"""
