rule downloadVirsorter:
	output:
		virSorter_db=config['virSorter_dir'],
		virSorter_dir=config['virSorter_db']
	message:
		"Downloading required VirSorter files"
	threads: 1
	params:
		virSorter_db="db/virSorter"
	shell:
		"""
		if [ ! -f {output.virSorter_dir} ]
		then
			mkdir -p tools
			git clone https://github.com/simroux/VirSorter.git 
			mv VirSorter tools
			cd tools/VirSorter/Scripts 
			make clean
			make
			cd ../../../
		fi
		if [ ! -f {output.virSorter_db} ]
		then
			curl -OL https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
			mkdir -p {params.virSorter_db}
			tar -xvzf virsorter-data-v2.tar.gz -C {params.virSorter_db}
			rm virsorter-data-v2.tar.gz
		fi
    	"""
rule virSorter:
	input:
		representatives=dirs_dict["vOUT_DIR"] + "/{sample}_merged_scaffolds_95-80.fna",
		virSorter_dir=config['virSorter_dir'],
		virSorter_db=config['virSorter_db']
	output:
		pvalues=dirs_dict["vOUT_DIR"] + "/{sample}_virSorter/VIRSorter_global-phage-signal.csv"
	params:
		out_folder=dirs_dict["vOUT_DIR"] + "/{sample}_virSorter",
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
			--data-dir {params.virSorter_db}   
			--virome  
		"""

#rule virFinder:
#	input:
#		scaffolds=directory("{ASSEMBLY_DIR}/{sample}/filtered_scaffolds.fasta")
#	output:
#		pvalues = "{indir}/virfinder/{name}.pvalues.tsv"	
#	params:
#		virFinder_script="scripts/virfinder_wrapper.R'"
#	messsage: 
#		"Scoring virus VirFinder"
#	conda:
#		"envs/main.yaml"
#	threads: 1
#	shell:
#		"""
#		Rscript {params.wrapper_script} {input.fasta} {output.pvalues} 
#		cat virFinder_list.txt |  sed '1d' | awk '{if($4 <= 0.05) print $1}' > virFinder_selection.txt
#		"""
                  
#rule getVirSorterDB:
#	input:
#		directory("{DB_DIR}/")
#	output:
#		db = "{DB_DIR}/virsorter-data-v2"	
#	messsage: 
#		"Downloading Virsorter Data"
#	shell:
#		"""
#		wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
#		tar -xvzf virsorter-data-v2.tar.gz
#       """
