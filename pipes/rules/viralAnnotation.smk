rule get_VIGA:
	output:
		VIGA_dir=directory(os.path.join(workflow.basedir, config['viga_dir'])),
		piler_dir=directory(os.path.join(workflow.basedir, config['piler_dir'])),
		trf_dir=directory(os.path.join(workflow.basedir, config['trf_dir'])),
	message:
		"Downloading MMseqs2"
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 4
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone --depth 1 https://github.com/EGTortuero/viga.git
		chmod 744 viga/create_dbs.sh viga/VIGA.py
		./viga/create_dbs.sh
		wget https://www.drive5.com/pilercr/pilercr1.06.tar.gz --no-check-certificate
		tar -xzvf pilercr1.06.tar.gz
		cd pilercr1.06
		make
		cd ..
		mkdir TRF
		cd TRF
		wget http://tandem.bu.edu/trf/downloads/trf409.linux64
		wget http://tandem.bu.edu/irf/downloads/irf307.linux.exe
		mv trf409.linux64 trf
		mv irf307.linux.exe irf
		chmod 744 trf irf
		"""

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
		GenBank_table=dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + "_annotated.tbl",
		GenBank_fasta=dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + "_annotated.fasta",
		csv=dirs_dict["ANNOTATION"] + "/" + REFERENCE_CONTIGS_BASE + ".tot" + "_annotated.csv",
		viga_log=dirs_dict["ANNOTATION"] + "/viga_log_" + REFERENCE_CONTIGS_BASE + ".tot.txt",
		viga_names=temp(dirs_dict["ANNOTATION"] + "/viga_names_" + REFERENCE_CONTIGS_BASE + ".tot.txt"),
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REFERENCE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
		VIGA_dir=directory("../" + config['viga_dir']),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	message:
		"Creating databases for reference and assembly mmseqs"
	threads: 16
	shell:
		"""
		WD=$(pwd)
		PILER=$WD/{input.piler_dir}
		PATH=$PILER:$PATH
		TRF=$WD/{input.trf_dir}
		PATH=$TRF:$PATH
		ln -sfn {input.positive_contigs} {output.temp_symlink}
		mkdir -p {output.temp_viga_dir}
		cd {output.temp_viga_dir}
		touch {output.modifiers}
		{input.VIGA_dir}/VIGA.py --input {output.temp_symlink} --diamonddb {input.VIGA_dir}/databases/RefSeq_Viral_DIAMOND/refseq_viral_proteins.dmnd \
		--blastdb {input.VIGA_dir}/databases/RefSeq_Viral_BLAST/refseq_viral_proteins --hmmerdb {input.VIGA_dir}/databases/pvogs/pvogs.hmm \
		--rfamdb {input.VIGA_dir}/databases/rfam/Rfam.cm --modifiers {output.modifiers} --threads {threads} &> {output.viga_log}
		cat {output.viga_log}| grep "was renamed as" > {output.viga_names}
		cat {output.viga_names} | while read line
		do
			stringarray=($line)
			new=${{stringarray[-1]}}
			old=${{stringarray[1]}}
			sed -i -e "s/${{new}}\t/${{old}}\t/g" -e "s/${{new}}_/${{old}}_/g" {output.csv}
			sed -i -e "s/${{new}}$/${{old}}/g" -e "s/${{new}} /${{old}} /g" -e "s/${{new}}_/${{old}}_/g" {output.GenBank_file}
			sed -i -e "s/${{new}}$/${{old}}/g" -e "s/${{new}} /${{old}} /g" -e "s/${{new}}_/${{old}}_/g" {output.GenBank_table}
			sed -i "s/>${{new}} $/>${{old}}/g" {output.GenBank_fasta}
		done
		"""

rule get_mmseqs:
	output:
		MMseqs2_dir=directory(config['mmseqs_dir']),
	message:
		"Downloading MMseqs2"
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 4
	shell:
		"""
		MM_dir={output.MMseqs2_dir}
		echo $MM_dir
		if [ ! -d $MM_dir ]
		then
			mkdir -p tools
			cd tools
			git clone https://github.com/soedinglab/MMseqs2.git
			cd MMseqs2
			mkdir build
			cd build
			cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
			make -j {threads}
			make install
		fi
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
rule get_VIBRANT:
	output:
		VIBRANT_dir=os.path.join(workflow.basedir, config['vibrant_dir']),
	message:
		"Downloading VIBRANT"
	threads: 4
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/AnantharamanLab/VIBRANT
		chmod -R 744 VIBRANT
		./VIBRANT/databases/VIBRANT_setup.sh
		git clone https://github.com/python/cpython
		cd cpython
		./configure
		make
		make test
		"""
rule annotate_VIBRANT:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/" + REFERENCE_CONTIGS_BASE + ".tot.fasta",
		VIBRANT_dir=os.path.join(workflow.basedir, config['vibrant_dir']),
	output:
		temp_vibrant=temp("tempVIBRANT.txt"),
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REFERENCE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
		VIGA_dir=directory("../" + config['viga_dir']),
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	message:
		"Annotating contigs with VIBRANT"
	threads: 16
	shell:
		"""
		{input.VIBRANT_dir}/VIBRANT_run.py -i {input.positive_contigs} -virome -t 16
		touch {output.temp_vibrant}
		"""
