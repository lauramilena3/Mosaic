rule getORFs:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + REFERENCE_CONTIGS_BASE + ".tot.fasta",
	output:
		coords=dirs_dict["VIRAL_DIR"]+ "/" + REFERENCE_CONTIGS_BASE + ".{sampling}.coords",
		aa=dirs_dict["VIRAL_DIR"]+ "/" + REFERENCE_CONTIGS_BASE + "_ORFs.{sampling}.fasta",
	message:
		"Calling ORFs with prodigal"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		if [ ! -s {input.positive_contigs} ]
		then
			prodigal -i {input.positive_contigs} -o {output.coords} -a {output.aa} -p meta
		else
			echo "Empty contigs file, no ORFs to detect"
			touch {output.coords} {output.aa}
		fi
		"""
rule clusterTaxonomy:
	input:
		aa=dirs_dict["VIRAL_DIR"]+ "/" + REFERENCE_CONTIGS_BASE + "_ORFs.{sampling}.fasta",
	output:
		genome_file=dirs_dict["VIRAL_DIR"]+ "/" + REFERENCE_CONTIGS_BASE + "_vContact.{sampling}/genome_by_genome_overview.csv",
	params:
		clusterONE_dir=config["clusterONE_dir"],
		vcontact_dir=config["vcontact_dir"],
		out_dir=directory(dirs_dict["VIRAL_DIR"]+ "/" + REFERENCE_CONTIGS_BASE + "_vContact.{sampling}"),
	message:
		"Clustering viral genomes with vContact2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 8
	shell:
		"""
		# if [ ! -d {params.vcontact_dir} ]
		# the
		# 	git clone https://bitbucket.org/MAVERICLab/vcontact2/
		# 	mv vcontact2 tools
		# 	envir=$( which vcontact | rev | cut -d/ -f3 | rev)
		# 	cp {params.vcontact_dir}/vcontact/data/ViralRefSeq-* .snakemake/conda/$envir/lib/python3.7/site-packages/vcontact/data/
		# 	cp scripts/matrices.py .snakemake/conda/$envir/lib/python3.7/site-packages/vcontact/matrices.py
		# 	cp scripts/vcontact .snakemake/conda/$envir/bin/vcontact
		# 	cp scripts/summaries.py .snakemake/conda/$envir/lib/python3.7/site-packages/vcontact/exports/summaries.py
		# fi
		#three changes in code 1) int 2,3) summary remove excluded
		vcontact2_gene2genome -p {input.aa} -s Prodigal-FAA -o {output.genome_file}
		vcontact --raw-proteins {input.aa} --rel-mode 'Diamond' --proteins-fp {output.genome_file} \
		--db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {params.clusterONE_dir}/cluster_one-1.0.jar \
		--output-dir {params.out_dir} --threads {threads}
		"""

rule mmseqsTaxonomy:
	input:
		representatives=dirs_dict["VIRAL_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		mmseqs_dir=(config['mmseqs_dir']),
		refseq=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna")),
		refseq_taxid=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna.taxidmapping")),
	output:
		mmseqsdir=directory(dirs_dict["VIRAL_DIR"] + "/taxonomy_mmseqs"),
		html=(dirs_dict["VIRAL_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.html"),
		tsv=(dirs_dict["VIRAL_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tsv"),
		table=(dirs_dict["VIRAL_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tbl"),
	message:
		"Taxonomy Assignment with MMseqs2"
	params:
		taxdump=(os.path.join(workflow.basedir,"db/ncbi-taxdump/")),
		refDB=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fnaDB")),
		current_dir=dirs_dict["VIRAL_DIR"],
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 1
	shell:
		"""
		#analyse
		{input.mmseqs_dir}/build/bin/mmseqs createdb {input.representatives} {params.current_dir}/positive_contigsDB
		{input.mmseqs_dir}/build/bin/mmseqs taxonomy {params.current_dir}/positive_contigsDB {params.refDB} \
			{params.current_dir}/taxonomyResult {params.current_dir}/tmp --search-type 2
		#results
		{input.mmseqs_dir}/build/bin/mmseqs createtsv {params.current_dir}/positive_contigsDB {params.current_dir}/taxonomyResult {output.tsv}
		{input.mmseqs_dir}/build/bin/mmseqs taxonomyreport {params.refDB} {params.current_dir}/taxonomyResult {output.table}
		{input.mmseqs_dir}/build/bin/mmseqs taxonomyreport {params.refDB} {params.current_dir}/taxonomyResult {output.html} --report-mode 1
	 	"""
