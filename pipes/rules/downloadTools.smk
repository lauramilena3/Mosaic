rule downloadContaminants:
	output:
		contaminant_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.fasta",
		contaminant_dir=temp(directory(dirs_dict["CONTAMINANTS_DIR"] +"/temp_{contaminant}")),
	message:
		"Downloading contaminant genomes"
	params:
		contaminants_dir=dirs_dict["CONTAMINANTS_DIR"],
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml",
	threads:
		16
	shell:
		"""
		mkdir {output.contaminant_dir}
		cd {output.contaminant_dir}
		wget $(esearch -db "assembly" -query {wildcards.contaminant} | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{{print $0"/"$NF"_genomic.fna.gz"}}')
		gunzip -f *gz
		cat *fna >> {output.contaminant_fasta}
		"""

rule get_VIBRANT:
	output:
		VIBRANT_dir=directory(os.path.join(workflow.basedir, config['vibrant_dir'])),
	message:
		"Downloading VIBRANT"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/AnantharamanLab/VIBRANT
		chmod -R 744 VIBRANT
		cd VIBRANT/databases
		./VIBRANT_setup.py
		"""
rule get_mmseqs:
	output:
		MMseqs2_dir=directory(config['mmseqs_dir']),
	message:
		"Downloading MMseqs2"
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 1
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
rule get_VIGA:
	output:
		VIGA_dir=directory(os.path.join(workflow.basedir, config['viga_dir'])),
		piler_dir=directory(os.path.join(workflow.basedir, config['piler_dir'])),
		trf_dir=directory(os.path.join(workflow.basedir, config['trf_dir'])),
		viga_db_dir=directory("tools/databases"),
	message:
		"Downloading MMseqs2"
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 1
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
rule downloadViralTools:
	output:
		virSorter_dir=directory(config['virSorter_dir']),
		virFinder_dir=directory(config['virFinder_dir']),
	message:
		"Downloading required VirSorter and VirFinder"
	threads: 1
	shell:
		"""
		#VIRSORTER
		VS_dir="{config[virSorter_dir]}"
		echo $VS_dir
		if [ ! -d $VS_dir ]
		then
			mkdir -p tools
			cd tools
			git clone https://github.com/simroux/VirSorter.git
			cd VirSorter/Scripts
			make clean
			make
			cd ../../../
		fi
		#VIRFNDER
		VF_dir="{config[virFinder_dir]}"
		echo $VF_dir
   		if [ ! -d $VF_dir ]
		then
			if [ ! {config[operating_system]} == "linux" ]
			then
				curl -OL https://raw.github.com/jessieren/VirFinder/blob/master/mac/VirFinder_1.1.tar.gz?raw=true
			else
				curl -OL https://github.com/jessieren/VirFinder/blob/master/linux/VirFinder_1.1.tar.gz?raw=true
			fi
			mkdir -p {output.virFinder_dir}
			mv VirFinder*tar.gz* {output.virFinder_dir}/VirFinder_1.1.tar.gz
		fi
		"""

rule downloadViralDB:
	output:
		virSorter_db=directory(config['virSorter_db']),
	message:
		"Downloading VirSorter database"
	threads: 1
	params:
		virSorter_db="db/VirSorter"
	shell:
		"""
		VS_db="{config[virSorter_db]}"
		echo $VS_db
		if [ ! -d $VS_db ]
		then
			curl -OL https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
			mkdir -p {params.virSorter_db}
			tar -xvzf virsorter-data-v2.tar.gz -C {params.virSorter_db}
			rm virsorter-data-v2.tar.gz
		fi
		"""
