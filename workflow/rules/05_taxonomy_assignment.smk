rule getORFs:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		coords=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.coords",
		aa=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs.{sampling}.fasta",
	message:
		"Calling ORFs with prodigal"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		if [ -s {input.representatives} ]
		then
			prodigal -i {input.representatives} -o {output.coords} -a {output.aa} -p meta
		else
			echo "Empty contigs file, no ORFs to detect"
			touch {output.coords} {output.aa}
		fi
		"""
rule clusterTaxonomy:
	input:
		aa=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs.{sampling}.fasta",
		clusterONE_dir=config["clusterONE_dir"],
		gene2genome_format_csv=(os.path.join(workflow.basedir,"db/vcontact2/gene-to-genome.30May2020.csv")),
		vcontact_format_aa=(os.path.join(workflow.basedir,"db/vcontact2/vcontact_format_30May2020.faa")),
	output:
		gene2genome=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/gene2genome.csv",
		merged_gene2genome=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/gene2genome_merged.csv",
		merged_ORFs=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_ORFs_merged.{sampling}.fasta",
		genome_file=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/genome_by_genome_overview.csv",
		viral_cluster_overview=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/viral_cluster_overview.csv",
	params:
		vcontact_dir=config["vcontact_dir"],
		out_dir=directory(dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}"),
	message:
		"Clustering viral genomes with vContact2"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 64
	shell:
		"""
		vcontact2_gene2genome -p {input.aa} -s Prodigal-FAA -o {output.gene2genome}
		cat {output.gene2genome} {input.gene2genome_format_csv}  > {output.merged_gene2genome}
		cat {input.aa} {input.vcontact_format_aa} > {output.merged_ORFs}
		dos2unix {output.merged_gene2genome}
		vcontact2 --raw-proteins {output.merged_ORFs} --rel-mode 'Diamond' --proteins-fp {output.merged_gene2genome} \
		--db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {input.clusterONE_dir}/cluster_one-1.0.jar \
		--output-dir {params.out_dir} --threads {threads} || true
		"""
rule parseVcontact:
	input:
		viral_cluster_overview=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/viral_cluster_overview.csv",
		formatting_taxonomy_affiliations=config["taxonomy_file"],
	output:
		taxonomy_results=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vcontact2_taxonomy.{sampling}.csv",
	message:
		"Assigning viral taxonomy with vContact2 results"
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

rule mmseqsTaxonomy:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		mmseqs_dir=directory(os.path.join(workflow.basedir, config['mmseqs_dir'])),
		refseq=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna")),
		refseq_taxid=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna.taxidmapping")),
	output:
		mmseqsdir=directory(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/"),
		html=(dirs_dict["vOUT_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.html"),
		tsv=(dirs_dict["vOUT_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tsv"),
		table=(dirs_dict["vOUT_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tbl"),
	message:
		"Taxonomy Assignment with MMseqs2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	params:
		positive_contigsDB=directory(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/positive_contigsDB"),
		taxonomyResultDB=directory(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/taxonomyResultDB"),
		tmp=directory(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/tmp"),
		taxdump=(os.path.join(workflow.basedir,"db/ncbi-taxdump/")),
		refDB=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fnaDB")),
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 8
	shell:
		"""
		#analyse
		#mkdir {output.mmseqsdir}
		mmseqs createdb {input.representatives} {params.positive_contigsDB}
		mmseqs taxonomy --threads {threads} {params.positive_contigsDB} {params.refDB} \
			{params.taxonomyResultDB} {params.tmp} --search-type 3 --lca-mode 2 -c 0.3 --cov-mode 2
		#results
		mmseqs createtsv {params.positive_contigsDB} {params.taxonomyResultDB} {output.tsv}
		mmseqs taxonomyreport {params.refDB} {params.taxonomyResultDB} {output.table}
		mmseqs taxonomyreport {params.refDB} {params.taxonomyResultDB} {output.html} --report-mode 1
	 	"""
