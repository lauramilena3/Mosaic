# ruleorder: gbk_to_faa>getORFs_assembly
def input_representative(wildcards):
	if METAGENOME:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{{sampling}}/quality_summary.tsv",sample=SAMPLES)),
	if NANOPORE:
		input_list.extend(expand(dirs_dict["vOUT_DIR"] + "/nanopore_{sample}_" + LONG_ASSEMBLER + "_checkV_{{sampling}}/quality_summary.tsv", sample=NANOPORE_SAMPLES)),
	if CROSS_ASSEMBLY:
		input_list.append(dirs_dict["vOUT_DIR"] + "/ALL_checkV_{sampling}/quality_summary.tsv"),
	if SUBASSEMBLY:
		input_list.extend(expand(dirs_dict["ASSEMBLY_TEST"] + "/{sample}_{subsample}_geNomad_checkV_{{sampling}}/quality_summary.tsv", sample=SAMPLES, subsample=subsample_test)),
	if len(config['additional_reference_contigs'])>0:
		input_list.append(dirs_dict["vOUT_DIR"] + "/user_reference_contigs_checkV/quality_summary.tsv"),
	return input_list


rule getORFs_prodigal_gv:
	input:
		nuc_fasta="{fasta}.{sampling}.fasta",
	output:
		coords="{fasta}_ORFs.{sampling}.coords",
		aa="{fasta}_ORFs.{sampling}.faa",
	message:
		"Calling ORFs with prodigal-gv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_orf.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/getORFs_prodigal_gv/fasta={fasta}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		python ./scripts/parallel-prodigal-gv.py -q -i {input.nuc_fasta} -o {output.coords} -a {output.aa} -t {threads}
		"""

rule getORFs_coding_length:
	input:
		aa="{fasta}_ORFs.{sampling}.faa",
	output:
		length_temp=temp("{fasta}_ORFs_length_temp.{sampling}.txt"),
		names_temp=temp("{fasta}_ORFs_names_temp.{sampling}.txt"),
		length="{fasta}_ORFs_length.{sampling}.txt",
		cummulative_length="{fasta}_ORFs_coding_lengths.{sampling}.txt",
	message:
		"Calculating ORFs length"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/getORFs_coding_length/fasta={fasta}__sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		awk '$0 ~ ">" {{if (NR > 1) print c; c = 0; printf substr($0,2,100) "\t"}} $0 !~ ">" {{c += length($0)}} END {{print c}}' {input.aa} > {output.length_temp}
		cat {output.length_temp} | cut -f1 | cut -f1 -d' ' | rev | cut -d_ -f2- | rev > {output.names_temp}
		paste {output.length_temp} {output.names_temp} > {output.length}
		cat {output.length} | cut -f3,2 | awk '{{arr[$2]+=$1}} END {{for (i in arr) {{print i,arr[i]}}}}' > {output.cummulative_length}
		"""

rule clusterTaxonomy:
	input:
		aa=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs.{sampling}.faa",
		clusterONE_dir=config["clusterONE_dir"],
		gene2genome_millard=("db/vcontact2/GenomesDB_Aug_2026_vConTACT2_gene_to_genome.csv"),
		vcontact_aa_millard=("db/vcontact2/GenomesDB_Aug_2026_vConTACT2_proteins.faa"),
	output:
		gene2genome=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/gene2genome.csv",
		merged_gene2genome=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/gene2genome_merged.csv",
		merged_ORFs=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs_merged.{sampling}.fasta",
		genome_file=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/genome_by_genome_overview.csv",
		viral_cluster_overview=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/viral_cluster_overview.csv",
	params:
		out_dir=directory(dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}"),
		reference_genomes=config["reference_genomes_vcontact"],
		#reference_genomes='ProkaryoticViralRefSeq94-Merged'
	message:
		"Clustering viral genomes with vContact2"
	conda:
		dirs_dict["ENVS_DIR"] + "/vcontact.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/clusterTaxonomy/sampling={sampling}.tsv"
	# conda:
	# 	dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	threads: 64
	shell:
		"""
		rm -rf {params.out_dir}
		mkdir {params.out_dir}
		vcontact2_gene2genome -p {input.aa} -s Prodigal-FAA -o {output.gene2genome}
		cat {output.gene2genome} {input.gene2genome_millard}  > {output.merged_gene2genome}
		cat {input.aa} {input.vcontact_aa_millard} > {output.merged_ORFs}
		dos2unix {output.merged_gene2genome}
		vcontact2 --raw-proteins {output.merged_ORFs} --rel-mode 'Diamond' --proteins-fp {output.merged_gene2genome} \
		--db {params.reference_genomes} --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {input.clusterONE_dir}/cluster_one-1.0.jar \
		--output-dir {params.out_dir} --threads {threads} -f
		#|| true
		"""

rule parseVcontact:
	input:
		viral_cluster_overview=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/viral_cluster_overview.csv",
		formatting_taxonomy_affiliations=config["taxonomy_file"],
	output:
		taxonomy_results=dirs_dict["ANNOTATION"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + "_vcontact2_taxonomy.{sampling}.csv",
	message:
		"Assigning viral taxonomy with vContact2 results"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/parseVcontact/sampling={sampling}.tsv"
	threads: 1
	run:
		import pandas as pd
		def is_unique(s):
			a = s.to_numpy()
			return (a[0] == a).all()

		taxonomy_df=pd.read_csv(input.formatting_taxonomy_affiliations)
		df=pd.read_csv(input.viral_cluster_overview, index_col=0)

		df[['cluster','subcluster']] = df.VC.str.rsplit('_', 1, expand=True)
		grouped_df=df.groupby('cluster')
		grouped_results_df=pd.DataFrame()
		members=[]
		vc=[]
		for name, group in grouped_df:
			members.append(','.join([str(elem) for elem in group.Members.tolist()]))
			vc.append(name)
		grouped_results_df["Members"]=members
		grouped_results_df["VC"]=vc


		df=grouped_results_df[grouped_results_df['Members'].str.contains("NODE|tig0")]
		df["Members"]=df["Members"].str.split(",")
		accessions=[]
		nodes=[]
		taxonomies=[]

		with open(output.taxonomy_results, 'w') as f:
			for index, row in df.iterrows():
				accession1=([x for x in row['Members'] if not 'NODE' in x])
				accession=[item for item in accession1 if not item.startswith("tig00")]
				accession=([x for x in accession if not '~' in x])
				accessions.append(accession)
				#node=[x for x in row['Members'] if 'NODE' in x]
				node1=([x for x in row['Members'] if 'NODE' in x])
				node2=([x for x in row['Members'] if 'tig00' in x])
				node=node1+node2
				print(node)
				nodes.append(node)
				taxonomy=[]
				for acc in accession:
					print(acc)
					print(taxonomy_df[taxonomy_df["acc"]==acc]["lineage"].values[0].split(";"))
					taxonomy.append(taxonomy_df[taxonomy_df["acc"]==acc]["lineage"].values[0].split(";"))
				if taxonomy:
					tax_df=pd.DataFrame(taxonomy)
					tax_df.columns=["kindom", "phylum", "class", "order", "family", "genus", "species"]
					#print(tax_df)
					tax_df=tax_df.drop(columns="species")
					consensus_tax=""
					for (columnName, columnData) in tax_df.iteritems():
						#print('Colunm Name : ', columnName)
						if (is_unique(tax_df[columnName])):
							if not tax_df[columnName].to_numpy()[0] =="__":
								consensus_tax=(columnName, tax_df[columnName].to_numpy()[0])
					#print(tax_df)
					for n in node:
						print(n + "\t" + consensus_tax[1] +" [" + consensus_tax[0] + "]", file=f)

rule PhaGCNTaxonomy:
	input:
		phabox_db=config["phabox_db"],
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
	output:
		taxonomy_table=dirs_dict["ANNOTATION"] + "/PhaGCN_taxonomy_report_{sequence}.csv",
		raw_taxonomy=dirs_dict["ANNOTATION"] + "/PhaBOX_{sequence}/phagcn_prediction.tsv",
	message:
		"Taxonomy assignment with PhaBOX 2 PhaGCN (taxonomy only)"
	conda:
		dirs_dict["ENVS_DIR"] + "/phabox.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/PhaGCNTaxonomy/sequence={sequence}.tsv"
	params:
		minimum_length=int(config.get("phabox_min_length", 2000)),
	log:
		dirs_dict["ANNOTATION"] + "/PhaBOX_{sequence}.log",
	threads: 32
	wildcard_constraints:
		  sequence="[^/]+"  # The 'sequence' wildcard cannot contain a slash
	shell:
		"""
		mkdir -p "$(dirname {output.raw_taxonomy:q})"
		phabox_work=$(mktemp -d "$(dirname {output.raw_taxonomy:q})/phabox_work.XXXXXX")
		trap 'rm -rf -- "$phabox_work"' EXIT
		if [ -s {input.fasta:q} ]; then
			phabox2 --task phagcn --contigs {input.fasta:q} --dbdir {input.phabox_db:q} \
				--outpth "$phabox_work" --len {params.minimum_length} --threads {threads} > {log:q} 2>&1
		else
			mkdir -p "$phabox_work/final_prediction"
			printf 'Accession\\tLength\\tLineage\\tPhaGCNScore\\tGenus\\tGenusCluster\\n' > "$phabox_work/final_prediction/phagcn_prediction.tsv"
			printf 'Empty catalogue: no PhaBOX taxonomy search needed.\\n' > {log:q}
		fi
		python - "$phabox_work" {output.taxonomy_table:q} {input.phabox_db:q} {params.minimum_length} <<-'PYTHON'
		from pathlib import Path
		import sys
		import pandas as pd
		work, output_path, database, minimum = sys.argv[1:]
		work = Path(work)
		frame = pd.read_csv(work / "final_prediction/phagcn_prediction.tsv", sep="\\t", dtype=str).fillna("")
		classification_path = work / "midfolder/phagcn_classification.tsv"
		original = pd.read_csv(classification_path, sep="\\t", dtype=str).fillna("").set_index("Accession") if classification_path.is_file() else frame.set_index("Accession")
		ranks = ["superkingdom", "realm", "subrealm", "kingdom", "subkingdom", "phylum", "subphylum", "class", "subclass", "order", "suborder", "family", "subfamily", "genus", "subgenus", "species"]
		def rank_values(lineage):
		    return dict(part.strip().split(":", 1) for part in str(lineage).split(";") if ":" in part)
		rows = []
		for _, row in frame.iterrows():
		    accession = row["Accession"]
		    lineage = row.get("Lineage", "")
		    raw_score = row.get("PhaGCNScore", "")
		    values = rank_values(lineage)
		    source = original.loc[accession] if accession in original.index else row
		    original_lineage = source.get("Lineage", "")
		    original_score = source.get("Score", source.get("PhaGCNScore", ""))
		    original_values = rank_values(original_lineage)
		    score_tokens = str(original_score).split(";")
		    scores = dict(zip(original_values, score_tokens)) if len(original_values) == len(score_tokens) else dict()
		    family_score = pd.to_numeric(scores.get("family"), errors="coerce") if values.get("family") and values.get("family") == original_values.get("family") else None
		    status = "below_length_cutoff" if int(row["Length"]) < int(minimum) else ("classified" if values else (lineage if lineage not in ["", "-"] else "not reported"))
		    evidence = dict(contig_name=accession, prediction=values.get("family"), score=family_score,
		        PhaBOX_length_bp=row["Length"], PhaBOX_lineage=lineage, PhaBOX_score_raw=raw_score,
		        PhaBOX_original_lineage=original_lineage, PhaBOX_original_score_raw=original_score,
		        PhaBOX_genus_cluster=row.get("GenusCluster"), PhaBOX_prokaryotic_virus=row.get("Prokaryotic virus (Bacteriophages and Archaeal virus)"),
		        PhaBOX_classification_status=status, PhaBOX_software_version="2.2", PhaBOX_database=database,
		        PhaBOX_database_version="2.2", PhaBOX_min_length_bp=int(minimum))
		    for rank in ranks:
		        evidence["PhaBOX_" + rank] = values.get(rank)
		    rows.append(evidence)
		columns = ["contig_name", "prediction", "score", "PhaBOX_length_bp", "PhaBOX_lineage", "PhaBOX_score_raw", "PhaBOX_original_lineage", "PhaBOX_original_score_raw", "PhaBOX_genus_cluster", "PhaBOX_prokaryotic_virus", "PhaBOX_classification_status", "PhaBOX_software_version", "PhaBOX_database", "PhaBOX_database_version", "PhaBOX_min_length_bp"] + ["PhaBOX_" + rank for rank in ranks]
		pd.DataFrame(rows, columns=columns).to_csv(output_path, index=False, na_rep="not reported")
		PYTHON
		cp "$phabox_work/final_prediction/phagcn_prediction.tsv" {output.raw_taxonomy:q}
	 	"""

rule hostID_iphop:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/{sequence}.fasta",
		iphop_db=(config['iphop_db']),
	output:
		results_dir=directory(dirs_dict["ANNOTATION"] + "/iphop_hostID_{sequence}_resultsDir"),
	message:
		"Host finding with iphop"
	conda:
		dirs_dict["ENVS_DIR"] + "/env2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/hostID_iphop/sequence={sequence}.tsv"
	threads: 64
	wildcard_constraints:
		  sequence="[^/]+"  # The 'sequence' wildcard cannot contain a slash
	shell:
		"""
		iphop predict --fa_file {input.fasta} --db_dir {input.iphop_db} --out_dir {output.results_dir} --num_threads {threads}
		rm -rf {output.results_dir}/Wdir
		"""

rule single_fasta_filtered:
	input:
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		filtered_representatives_checkpoint=temp(dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/filtered_representatives_dir_checkpoint.txt"),
	params:
		filtered_representatives_dir=((dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}")),
	message:
		"formating filtered vOTUs into single fasta"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/single_fasta_filtered/sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		seqkit split --quiet -w 0 -i {input.filtered_representatives} --out-dir {params.filtered_representatives_dir}
		touch {output.filtered_representatives_checkpoint}
	 	"""

rule match_spacers:
	input:
		spacers=config['microbial_spacers'],
		filtered_representatives_checkpoint=(dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/filtered_representatives_dir_checkpoint.txt"),
	output:
		spacer_match=dirs_dict["ANNOTATION"] + "/spacepharer_minced_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tsv",
	message:
		"Matching microbial spacers"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_spacepharer.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/match_spacers/sampling={sampling}.tsv"
	params:
		filtered_representatives_dir=((dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}")),
		viralTargetDB=temp(directory(dirs_dict["ANNOTATION"] + "/viralTargetDB.{sampling}")),
		viralTargetDB_rev=temp(directory(dirs_dict["ANNOTATION"] + "/viralTargetDB_rev.{sampling}")),
		spacers_mincedSetDB=temp(directory(dirs_dict["ANNOTATION"] + "/spacers_mincedSetDB.{sampling}")),
		tmpFolder=temp(directory(dirs_dict["ANNOTATION"] + "/tmpFolder.{sampling}")),	
	threads: 1
	shell:
		"""
		ulimit -S -s unlimited
		rm -rf {params.spacers_mincedSetDB}* {params.viralTargetDB}* {params.viralTargetDB_rev}* {params.tmpFolder}*
		spacepharer createsetdb {params.filtered_representatives_dir}/*fasta {params.viralTargetDB} {params.tmpFolder}
		spacepharer createsetdb {params.filtered_representatives_dir}/*fasta {params.viralTargetDB_rev} {params.tmpFolder} --reverse-fragments 1
		spacepharer createsetdb {input.spacers} {params.spacers_mincedSetDB} {params.tmpFolder} --extractorf-spacer 1
		spacepharer predictmatch {params.spacers_mincedSetDB} {params.viralTargetDB} {params.viralTargetDB_rev} {output.spacer_match} {params.tmpFolder} -s 7.5 
		rm -rf {params.spacers_mincedSetDB}* {params.viralTargetDB}* {params.viralTargetDB_rev}* {params.tmpFolder}*
	 	"""

rule match_spacers_dion:
	input:
		spacers_dion_db=(os.path.join(workflow.basedir, config['dion_db'])),
		filtered_representatives_dir=dirs_dict["vOUT_DIR"]+ "/single_filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	output:
		spacer_match=dirs_dict["ANNOTATION"] + "/spacepharer_dion_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tsv",
		spacer_match_header=dirs_dict["ANNOTATION"] + "/spacepharer_dion_" + REPRESENTATIVE_CONTIGS_BASE + "_header.{sampling}.tsv",
	message:
		"Matching microbial spacers with the DION database"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_spacepharer.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/match_spacers_dion/sampling={sampling}.tsv"
	params:
		viralTargetDB=temp(directory(dirs_dict["ANNOTATION"] + "/viralTargetDB_dion.{sampling}")),
		viralTargetDB_rev=temp(directory(dirs_dict["ANNOTATION"] + "/viralTargetDB_dion_rev.{sampling}")),
		# spacers_dionSetDB=temp(directory(dirs_dict["ANNOTATION"] + "/spacers_dionSetDB.{sampling}")),
		tmpFolder=temp(directory(dirs_dict["ANNOTATION"] + "/tmpFolder_dion.{sampling}")),	
	threads: 1
	shell:
		"""
		rm -rf {params.viralTargetDB}* {params.viralTargetDB_rev}* {params.tmpFolder}*	 	
		spacepharer createsetdb {input.filtered_representatives_dir}/*fasta {params.viralTargetDB} {params.tmpFolder}
		spacepharer createsetdb {input.filtered_representatives_dir}/*fasta {params.viralTargetDB_rev} {params.tmpFolder} --reverse-fragments 1
		spacepharer predictmatch {input.spacers_dion_db}/dionSetDB {params.viralTargetDB} {params.viralTargetDB_rev} {output.spacer_match} {params.tmpFolder} -s 7.5 
		grep "^#" {output.spacer_match} > {output.spacer_match_header}
		rm -rf {params.viralTargetDB}* {params.viralTargetDB_rev}* {params.tmpFolder}*	 	
		"""


rule taxmyphage:
	input:
		fasta=lambda wc: annotation_fasta_path(wc.sequence),
		taxmyphage_db=(config['taxmyphage_db']),
	output:
		results_dir=directory(dirs_dict["ANNOTATION"] + "/taxmyphage_{sequence}"),
	params:
		results_dir=dirs_dict["ANNOTATION"] + "/taxmyphage_{sequence}",
		annotation_dir=directory(dirs_dict["ANNOTATION"]),
	message:
		"Assigning taxonomy with taxmyphage"
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/taxmyphage/sequence={sequence}.tsv"
	threads: 32
	shell:
		r"""
		if [ -s {input.fasta:q} ]; then
			taxmyphage run -i {input.fasta:q} -t {threads} -db {input.taxmyphage_db:q} -o {params.results_dir:q} --force
		else
			mkdir -p {output.results_dir:q}
			printf 'Genome\n' > {output.results_dir:q}/Summary_taxonomy.tsv
		fi
		"""
