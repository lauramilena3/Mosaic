rule downloadViralFiles:
	output:
		virSorter_db=protected(directory(config['virSorter_db'])),
		virSorter_dir=directory(config['virSorter_dir']),
		virFinder_dir=directory(config['virFinder_dir'])
	message:
		"Downloading required VirSorter and VirFinder data"
	threads: 1
	params:
		virSorter_db="db/VirSorter"
	shell:
		"""
		if [ ! -d {config[virSorter_dir]} ]
		then
			echo "no existe"
			mkdir -p tools
			cd tools
			git clone https://github.com/simroux/VirSorter.git 
			cd VirSorter/Scripts 
			make clean
			make
			cd ../../../
		fi
		if [ ! -d {config[virSorter_db]} ]
		then
			curl -OL https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
			mkdir -p {params.virSorter_db}
			tar -xvzf virsorter-data-v2.tar.gz -C {params.virSorter_db}
			rm virsorter-data-v2.tar.gz
		fi
   		if [ ! -d {config[virFinder_dir]} ]
		then
			if [ ! {config[operating_system]} == "linux" ] 
			then
				curl -OL https://raw.github.com/jessieren/VirFinder/blob/master/mac/VirFinder_1.1.tar.gz?raw=true
				echo "mac"
			else
				curl -OL https://github.com/jessieren/VirFinder/blob/master/linux/VirFinder_1.1.tar.gz?raw=true
			fi
			mkdir -p {output.virFinder_dir}
			mv VirFinder*tar.gz* {output.virFinder_dir}/VirFinder_1.1.tar.gz
		fi
    	"""

rule virSorter:
	input:
		representatives=dirs_dict["vOUT_DIR"] + "/{sample}_merged_scaffolds_95-80.fna",
		virSorter_dir=config['virSorter_dir'],
		virSorter_db=config['virSorter_db']
	output:
		results=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter/VIRSorter_global-phage-signal.csv"
	params:
		out_folder=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter",
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		{config[virSorter_dir]}/wrapper_phage_contigs_sorter_iPlant.pl -f {input.representatives} \
			--db 2 \
			--wdir {params.out_folder} \
			--ncpu {threads} \
			--data-dir {input.virSorter_db} \
			--virome  
		"""

rule virFinder:
	input:
		scaffolds=dirs_dict["vOUT_DIR"] + "/{sample}_merged_scaffolds_95-80.fna",
		virFinder_dir=config['virFinder_dir']
	output:
		pvalues=dirs_dict["VIRAL_DIR"] + "/{sample}_virFinder_pvalues.txt"
	params:
		virFinder_script="scripts/virfinder_wrapper.R'"
	message: 
		"Scoring virus VirFinder"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		./{params.virFinder_script} {input.scaffolds} {output.pvalues} 
		"""

rule getViralTable:
	input:
		pvalues = dirs_dict["VIRAL_DIR"] + "/{sample}_virFinder_pvalues.txt",
		categories=dirs_dict["vOUT_DIR"] + "/{sample}_virSorter/VIRSorter_global-phage-signal.csv"
	output:
		pvalues = "{indir}/virfinder/{name}.pvalues.tsv"	
	params:
		virFinder_script="scripts/virfinder_wrapper.R'",
		virFinder_dir=config['virFinder_dir']
	message: 
		"Scoring virus VirFinder"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		"""

rule extractViralContigs:


                  

