rule annotate_VIGA:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + REFERENCE_CONTIGS_BASE + ".tot.fasta",
		VIGA_dir=os.path.join(workflow.basedir, config['viga_dir']),
		piler_dir=os.path.join(workflow.basedir, (config['piler_dir'])),
		trf_dir=os.path.join(workflow.basedir, (config['trf_dir'])),
	output:
		modifiers=temp(dirs_dict["ANNOTATION"] + "/modifiers.txt"),
		temp_symlink=temp(dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot.fasta"),
		temp_viga_dir=temp(directory(dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + "_tempVIGA")),
		GenBank_file=dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + "_annotated.gbk",
		GenBank_table_temp1=temp(dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + "_annotated.tbl"),
		GenBank_table_temp2=temp(dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + "_annotated.tbl2"),
		GenBank_table=dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + ".tbl",
		GenBank_fasta=dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + "_annotated.fasta",
		csv=dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + "_annotated.csv",
		viga_log=dirs_dict["ANNOTATION"] + "/viga_log_" + REFERENCE_CONTIGS_BASE + ".tot.txt",
		viga_names=temp(dirs_dict["ANNOTATION"] + "/viga_names_" + REFERENCE_CONTIGS_BASE + ".tot.txt"),
		viga_topology_temp=temp(dirs_dict["ANNOTATION"] + "/viga_topology_temp" + REFERENCE_CONTIGS_BASE + "tot.txt"),
		viga_topology=(dirs_dict["ANNOTATION"] + "/viga_topology_" + REFERENCE_CONTIGS_BASE + "tot.txt"),
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REFERENCE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
		VIGA_dir=directory("../" + config['viga_dir']),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	message:
		"Annotating contigs with VIGA"
	threads: 64
	shell:
		"""
		PILER={input.piler_dir}
		PATH=$PILER:$PATH
		TRF={input.trf_dir}
		PATH=$TRF:$PATH
		ln -sfn {input.positive_contigs} {output.temp_symlink}
		mkdir -p {output.temp_viga_dir}
		cd {output.temp_viga_dir}
		touch {output.modifiers}
		{input.VIGA_dir}/VIGA.py --input {output.temp_symlink} --diamonddb {input.VIGA_dir}/databases/RefSeq_Viral_DIAMOND/refseq_viral_proteins.dmnd \
		--blastdb {input.VIGA_dir}/databases/RefSeq_Viral_BLAST/refseq_viral_proteins --hmmerdb {input.VIGA_dir}/databases/pvogs/pvogs.hmm \
		--rfamdb {input.VIGA_dir}/databases/rfam/Rfam.cm --modifiers {output.modifiers} --threads {threads} &> {output.viga_log}
		cat {output.viga_log} | grep "was renamed as" > {output.viga_names}
		cat {output.viga_log} | grep "according to LASTZ" > {output.viga_topology_temp}
		cat {output.viga_names} | while read line
		do
			stringarray=($line)
			new=${{stringarray[-1]}}
			old=${{stringarray[1]}}
			sed -i -e "s/${{new}}\t/${{old}}\t/g" -e "s/${{new}}_/${{old}}_/g" {output.csv}
			sed -i -e "s/${{new}}$/${{old}}/g" -e "s/${{new}} /${{old}} /g" -e "s/${{new}}_/${{old}}_/g" {output.GenBank_file}
			sed -i -e "s/${{new}}$/${{old}}/g" -e "s/${{new}} /${{old}} /g" -e "s/${{new}}_/${{old}}_/g" {output.GenBank_table_temp1}
			sed -i "s/>${{new}} $/>${{old}}/g" {output.GenBank_fasta}
			sed -i -e "s/${{new}} /${{old}} /g" {output.viga_topology_temp}
		done
		awk  '{{print $1 "\t" $6}}'  {output.viga_topology_temp} > {output.viga_topology}
		grep -v "gene$" {output.GenBank_table_temp1} > {output.GenBank_table_temp2}
		grep -n "CDS$" {output.GenBank_table_temp2} | cut -d : -f 1 | awk '{$1+=-1}1' | sed 's%$%d%' | sed -f - {output.GenBank_table_temp2} > {output.GenBank_table}
		"""

rule create_dbs_mmseqs2:
	input:
		MMseqs2_dir=(config['mmseqs_dir']),
		representatives=dirs_dict["vOUT_DIR"] + "/merged_scaffolds.tot_95-80.fna",
		reference=REFERENCE_CONTIGS
	output:
		index_representatives=dirs_dict["MMSEQS"] + "/representatives.index",
		index_reference=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE + ".index",
		idx_reference=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE + ".idx",
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REFERENCE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
	message:
		"Creating databases for reference and assembly mmseqs"
	threads: 4
	shell:
		"""
		{params.mmseqs}/mmseqs createdb {input.representatives} {params.representatives_name}
		{params.mmseqs}/mmseqs createdb {input.reference} {params.reference_name}
		{params.mmseqs}/mmseqs createindex {params.reference_name} tmp --search-type 2
		"""

rule search_contigs_mmseqs2:
	input:
		index_representatives=dirs_dict["MMSEQS"] + "/representatives.index",
		index_reference=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE + ".index",
		idx_reference=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE + ".idx",
	output:
		results_index=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE + "_search_results.index",
		results_table=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE + "_best_search_results.txt",
		temp_dir=temp(directory(dirs_dict["MMSEQS"] + "/tmp")),
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REFERENCE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
	message:
		"Comparing reference and assembly mmseqs"
	threads: 16
	shell:
		"""
		mkdir {output.temp_dir}
		{params.mmseqs}/mmseqs search {params.representatives_name} {params.reference_name} {params.results_name} {output.temp_dir} \
		--start-sens 1 --sens-steps 3 -s 7 --search-type 2 --threads {threads}
		{params.mmseqs}/mmseqs convertalis {params.representatives_name} {params.reference_name} {params.results_name} {output.results_table}
		"""
