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
	params:
		vcontact_dir=config["vcontact_dir"],
		out_dir=directory(dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}"),
	message:
		"Clustering viral genomes with vContact2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 8
	shell:
		"""
		vcontact2_gene2genome -p {input.aa} -s Prodigal-FAA -o {output.gene2genome}
		cat {output.gene2genome} {input.gene2genome_format_csv}  > {output.merged_gene2genome}
		cat {input.aa} {input.vcontact_format_aa} > {output.merged_ORFs}
		dos2unix {output.merged_gene2genome}
		vcontact2 --raw-proteins {output.merged_ORFs} --rel-mode 'Diamond' --proteins-fp {output.merged_gene2genome} \
		--db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {input.clusterONE_dir}/cluster_one-1.0.jar \
		--output-dir {params.out_dir} --threads {threads}
		"""
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
	params:
		positive_contigsDB=directory(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/positive_contigsDB"),
		taxonomyResultDB=directory(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/taxonomyResultDB"),
		tmp=directory(dirs_dict["vOUT_DIR"] + "/taxonomy_mmseqs_"+ REPRESENTATIVE_CONTIGS_BASE + ".{sampling}/tmp"),
		taxdump=(os.path.join(workflow.basedir,"db/ncbi-taxdump/")),
		refDB=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fnaDB")),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 8
	shell:
		"""
		#analyse
		#mkdir {output.mmseqsdir}
		{input.mmseqs_dir}/build/bin/mmseqs createdb {input.representatives} {params.positive_contigsDB}
		{input.mmseqs_dir}/build/bin/mmseqs taxonomy --threads {threads} {params.positive_contigsDB} {params.refDB} \
			{params.taxonomyResultDB} {params.tmp} --search-type 3 --lca-mode 2
		#results
		{input.mmseqs_dir}/build/bin/mmseqs createtsv {params.positive_contigsDB} {params.taxonomyResultDB} {output.tsv}
		{input.mmseqs_dir}/build/bin/mmseqs taxonomyreport {params.refDB} {params.taxonomyResultDB} {output.table}
		{input.mmseqs_dir}/build/bin/mmseqs taxonomyreport {params.refDB} {params.taxonomyResultDB} {output.html} --report-mode 1
	 	"""
