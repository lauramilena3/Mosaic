rule get_VIGA:
	output:
		VIGA_dir=directory(config['viga_dir']),
		piler_dir=directory(config['piler_dir']),
		trf_dir=directory(config['trf_dir']),
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
		VIGA_dir=directory(config['viga_dir']),
		piler_dir=directory(config['piler_dir']),
		trf_dir=directory(config['trf_dir']),
	output:
		viga_results=dirs_dict["ANNOTATION"] + "/viga_results.txt",
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
		mkdir tempVIGA
		cd tempVIGA
		{params.VIGA_dir}/VIGA.py --input {input.positive_contigs} --diamonddb {params.VIGA_dir}/databases/RefSeq_Viral_DIAMOND/refseq_viral_proteins.dmnd \
		--blastdb {params.VIGA_dir}/databases/RefSeq_Viral_BLAST/refseq_viral_proteins --hmmerdb {params.VIGA_dir}/databases/pvogs/pvogs.hmm \
		--rfamdb {params.VIGA_dir}/databases/rfam/Rfam.cm --modifiers modifiers.txt --threads {threads}
		touch {output.viga_results}
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
