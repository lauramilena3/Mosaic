# cd TOOLS
# git clone https://github.com/EGTortuero/viga
# chmod 744 viga/create_dbs.sh viga/VIGA.py
# ./viga/create_dbs.sh
#
# cd pilercr1.06
# make
# PILER=/home/lmf/apps/Mosaic/pipes/tools/pilercr1.06
# PATH=$PILER:$PATH
# TRF=/home/lmf/apps/Mosaic/pipes/tools/TRF
# PATH=$TRF:$PATH
#
# cd VIGA
#
# /home/lmf/apps/Mosaic/pipes/tools/viga/VIGA.py --input /home/lmf/03_COLIPHAGES/FASTA/similar.fasta --diamonddb /home/lmf/apps/Mosaic/pipes/tools/viga/databases/RefSeq_Viral_DIAMOND/ --blastdb /home/lmf/apps/Mosaic/pipes/tools/viga/databases/RefSeq_Viral_BLAST/ --hmmerdb /home/lmf/apps/Mosaic/pipes/tools/viga/databases/pvogs/pvogs.hmm --rfamdb /home/lmf/apps/Mosaic/pipes/tools/viga/databases/rfam/Rfam.cm --modifiers modifiers.txt --threads 16

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
			git clone https://github.com/soedinglab/MMseqs2.git
			cd MMseqs2
			mkdir build
			cd build
			cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
			make -j {threads}
			make install
			export PATH=$(pwd)/bin/:$PATH
		fi
		"""

rule compare_contigs_mmseqs2:
	input:
		MMseqs2_dir=(config['mmseqs_dir']),
		representatives=dirs_dict["vOUT_DIR"] + "/merged_scaffolds.tot_95-80.fna",
		reference=REFERENCE_CONTIGS
	output:
		index_representatives=dirs_dict["MMSEQS"] + "/merged_scaffolds_95-80.index",
		index_reference=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE + ".index",
		idx_reference=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE + ".idx",
		temp_dir=temp(directory(dirs_dict["MMSEQS"] + "/tmp"))
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REFERENCE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" + "search_results"
	message:
		"Comparing reference and assembly mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 4
	shell:
		"""
		mmseqs createdb {input.representatives} {params.representatives_name}
		mmseqs createdb {input.reference} {params.reference_name}
		mmseqs createindex {params.reference_name} tmp --search-type 2
		mkdir {output.temp_dir}
		mmseqs search {params.representatives_name} {params.reference_name} {params.results_name} {output.temp_dir} \
		-a -s 7.0 --search-type 2 --threads {threads}
		"""
