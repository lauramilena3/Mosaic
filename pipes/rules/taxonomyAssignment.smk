rule getORFs:
	input:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.fasta",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.fasta"
	output:
		high_coords=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.coords",
		low_coords=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.coords",
		high_aa=dirs_dict["VIRAL_DIR"]+ "/high_confidence_ORFs.{sampling}.fasta",
		low_aa=dirs_dict["VIRAL_DIR"]+ "/low_confidence_ORFs.{sampling}.fasta",
		#high_aa_temp=dirs_dict["VIRAL_DIR"]+ "/high_confidence_ORFs.{sampling}.fasta_temp",
		#low_aa_temp=dirs_dict["VIRAL_DIR"]+ "/low_confidence_ORFs.{sampling}.fasta_temp",
		#high_genome_file=dirs_dict["VIRAL_DIR"]+ "/high_confidence_genome_file.{sampling}.csv",
		#low_genome_file=dirs_dict["VIRAL_DIR"]+ "/low_confidence_genome_file.{sampling}.csv",
	message:
		"Calling ORFs with prodigal"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		prodigal -i {input.high_contigs} -o {output.high_coords} -a {output.high_aa} -p meta
		prodigal -i {input.low_contigs} -o {output.low_coords} -a {output.low_aa} -p meta
		"""
rule clusterTaxonomy:
	input:
		high_aa=dirs_dict["VIRAL_DIR"]+ "/high_confidence_ORFs.{sampling}.fasta",
		low_aa=dirs_dict["VIRAL_DIR"]+ "/low_confidence_ORFs.{sampling}.fasta",
	output:
		high_dir=dirs_dict["VIRAL_DIR"]+ "/high_confidence_vContact.{sampling}",
		low_dir=dirs_dict["VIRAL_DIR"]+ "/low_confidence_vContact.{sampling}",
		high_genome_file=dirs_dict["VIRAL_DIR"]+ "/high_confidence_genome_file.{sampling}.csv",
		low_genome_file=dirs_dict["VIRAL_DIR"]+ "/low_confidence_genome_file.{sampling}.csv",
	params:
		clusterONE_dir=config["clusterONE_dir"],
		vcontact_dir=config["vcontact_dir"]
	message:
		"Clustering viral genomes with vContact2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 4
	shell:
		"""
		if [ ! -d {params.clusterONE_dir} ]
		then 
			mkdir -p {params.clusterONE_dir}
			curl -OL  http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar
			mv cluster_one-1.0.jar {params.clusterONE_dir}
		fi
		if [ ! -d {params.vcontact_dir} ]
		then
			git clone https://bitbucket.org/MAVERICLab/vcontact2/
			mv vcontact2 tools
			envir=$( which vcontact | rev | cut -d/ -f3 | rev)
			cp {params.vcontact_dir}/vcontact/data/ViralRefSeq-* .snakemake/conda/$envir/lib/python3.7/site-packages/vcontact/data/
		fi
		python ./{params.vcontact_dir}/vcontact/utilities/Gene2Genome.py -p {input.high_aa} -s Prodigal-FAA -o {output.high_genome_file}
		python ./{params.vcontact_dir}/vcontact/utilities/Gene2Genome.py -p {input.low_aa} -s Prodigal-FAA -o {output.low_genome_file}
		vcontact --raw-proteins {input.high_aa} --rel-mode 'Diamond' --proteins-fp {output.high_genome_file} \
		--db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode MCL \
		--output-dir {output.high_dir} --threads {threads}
		vcontact --raw-proteins {input.low_aa} --rel-mode 'Diamond' --proteins-fp {output.low_genome_file} \
		--db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode MCL \
		--output-dir {output.low_dir} --threads {threads}
		"""