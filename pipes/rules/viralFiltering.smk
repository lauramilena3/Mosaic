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
		virFinder_script="scripts/virfinder_wrapper.R"
	message: 
		"Scoring virus VirFinder"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		Rscript {params.virFinder_script} {input.scaffolds} {output.pvalues}
		"""

rule getViralTable:
	input:
		pvalues = dirs_dict["VIRAL_DIR"] + "/{sample}_virFinder_pvalues.txt",
		categories=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter/VIRSorter_global-phage-signal.csv"
	output:
		circular_H=dirs_dict["VIRAL_DIR"]+ "High_confidence_circular_list.txt",
		circular_L=dirs_dict["VIRAL_DIR"]+ "Low_confidence_circular_list.txt",
		non_circular_H=dirs_dict["VIRAL_DIR"]+ "High_confidence_non_circular_list.txt",
		non_circular_L=dirs_dict["VIRAL_DIR"]+ "Low_confidence_non_circular_list.txt"
	params:
		virFinder_script="scripts/virfinder_wrapper.R'",
		virFinder_dir=config['virFinder_dir']
	message: 
		"Scoring virus VirFinder"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	run:
		import pandas as pd
		results = pd.DataFrame(columns=('lenght', 'circular','VS_cat', 'VF_score', 'VF_pval', 'pCAT'))

		#VirSorter
		VS = {input.categories}  
		with open(VS) as fp:  
		    line = fp.readline()
		    cnt = 1
		    while line:
		        if line.startswith("#"):
		            if (line.strip().split()[1].isdigit()):
		                category=(line.strip().split()[1])
		        else:
		            circular="N"
		            if "cov_" in line.split(",")[0]:
		                contigNameA=line.split(",")[0].split("VIRSorter_")[1].split("cov_")[0]
		                contigNameB=line.split(",")[0].split("VIRSorter_")[1].split("cov_")[1].replace("_", ".")
		                contigName=contigNameA + "cov_" + contigNameB
		            else:
		                contigName=line.split(",")[0].split("VIRSorter_")[1]
		            if "-circular" in contigName:
		                contigName=contigName.split("-circular")[0]
		                circular="Y"
		            results.loc[contigName, 'VS_cat'] = int(category)
		            results.loc[contigName, 'circular'] = circular
		        line = fp.readline()
		        cnt += 1
		        
		#VirFinder
		VF = {input.pvalues} 
		with open(VF) as fp:  
		    line = fp.readline()
		    cnt = 1
		    while line:
		        if cnt != 1:
		            contigName=line.split("\t")[0].strip()
		            contigLenght=line.split("\t")[1]
		            contigScore=line.split("\t")[2]
		            contigPval=line.split("\t")[3].split("\n")[0]
		            results.loc[contigName, 'lenght'] = float(contigLenght)
		            results.loc[contigName, 'VF_score'] = float(contigScore)
		            results.loc[contigName, 'VF_pval'] = float(contigPval)
		                #check if circular also
		        line = fp.readline()
		        cnt += 1    

		#filtering DFs
		df_A_c=results[results['VS_cat']<3][results['circular']=="Y"]
		df_B_c=results[results['VF_score']>0.9][results['VF_pval']<0.05][results['circular']=="Y"]
		df_C_c=results[results['VF_score']>0.7][results['VF_pval']<0.05][results['VS_cat']>0][results['circular']=="Y"]

		df_A_nc=results[results['VS_cat']<3][results['circular']=="N"]
		df_B_nc=results[results['VF_score']>0.9][results['VF_pval']<0.05][results['circular']=="N"]
		df_C_nc=results[results['VF_score']>0.7][results['VF_pval']<0.05][results['VS_cat']>0][results['circular']=="N"]

		#joinin and dereplicating circular contigs
		lsA=(df_A_c.index.tolist())
		lsB=(df_B_c.index.tolist())
		lsC=(df_C_c.index.tolist())

		lsAnotB=set(lsA) - set(lsB)
		lsAB=lsA + list(lsAnotB)

		lsCnoAB=set(lsC) - set(lsAB)
		lsCf=list(lsCnoAB)
		print("circular")
		print("\n".join(lsCf))

		f=open({input.circular_H}, 'w')
		f.write("\n".join(lsAB))
		f.close()

		f=open({input.circular_L}, 'w')
		f.write("\n".join(lsCf))
		f.close()

		#joinin and dereplicating noncircular contigs
		lsA=(df_A_nc.index.tolist())
		lsB=(df_B_nc.index.tolist())
		lsC=(df_C_nc.index.tolist())

		lsAnotB=set(lsA) - set(lsB)
		lsAB=lsA + list(lsAnotB)

		lsCnoAB=set(lsC) - set(lsAB)
		lsCf=list(lsCnoAB)

		f=open({input.non_circular_H}, 'w')
		f.write("\n".join(lsAB))
		f.close()

		f=open({input.non_circular_L}, 'w')
		f.write("\n".join(lsCf))
		f.close()

rule extractViralContigs:


                  

