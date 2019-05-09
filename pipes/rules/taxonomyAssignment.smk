rule getORFs:
	input:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.fasta",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.fasta"
	output:
		high_coords=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.coords",
		low_coords=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.coords",
		high_aa=dirs_dict["VIRAL_DIR"]+ "/high_confidence_ORFs.{sampling}.fasta",
		low_aa=dirs_dict["VIRAL_DIR"]+ "/low_confidence_ORFs.{sampling}.fasta",
		high_genome_file=dirs_dict["VIRAL_DIR"]+ "/high_confidence_genome_file.{sampling}.csv",
		low_genome_file=dirs_dict["VIRAL_DIR"]+ "/low_confidence_genome_file.{sampling}.csv",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		prodigal -i {input.high_contigs} -o {output.high_coords} -a {output.high_aa} -p meta
		prodigal -i {input.low_contigs} -o {output.low_coords} -a {output.low_aa} -p meta
		sed -i 's/;/|/g' {output.high_aa}
		sed -i 's/ # /|/g' {output.high_aa}
		sed -i 's/;/|/g' {output.low_aa}
		sed -i 's/ # /|/g' {output.low_aa}
		grep ">" {output.high_aa} | awk -F  "|" '{{ print substr($0,2,length($0))","substr($1,2,length($1))","$2";"$3";"$4";"$5";"$6}}' \
		> {output.high_genome_file}
		grep ">" {output.low_aa} | awk -F  "|" '{{ print substr($0,2,length($0))","substr($1,2,length($1))","$2";"$3";"$4";"$5";"$6}}' \
		> {output.low_genome_file}
		"""
rule clusterTaxonomy:
	input:
		high_aa=dirs_dict["VIRAL_DIR"]+ "/high_confidence_ORFs.{sampling}.fasta",
		low_aa=dirs_dict["VIRAL_DIR"]+ "/low_confidence_ORFs.{sampling}.fasta",
		high_genome_file=dirs_dict["VIRAL_DIR"]+ "/high_confidence_genome_file.{sampling}.csv",
		low_genome_file=dirs_dict["VIRAL_DIR"]+ "/low_confidence_genome_file.{sampling}.csv",
	output:
		high_dir=dirs_dict["VIRAL_DIR"]+ "/high_confidence_vContact.{sampling}",
		low_dir=dirs_dict["VIRAL_DIR"]+ "/low_confidence_vContact.{sampling}",
	params:
		clusterONE_dir=config["clusterONE_dir"]
	message:
		"Clustering viral genomes with vContact2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 1
	shell:
		"""
		if [ ! -d {params.clusterONE_dir} ]
		then 
			mkdir -p {params.clusterONE_dir}
			curl -OL  http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar
			mv cluster_one-1.0.jar {params.clusterONE_dir}
		fi
		vcontact --raw-proteins {input.high_aa} --rel-mode 'Diamond' --proteins-fp {input.high_genome_file} \
		--db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {params.clusterONE_dir} \
		--output-dir {output.high_dir}
		vcontact --raw-proteins {input.low_aa} --rel-mode 'Diamond' --proteins-fp {input.low_genome_file} \
		--db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {params.clusterONE_dir} \
		--output-dir {output.low_dir}
		"""