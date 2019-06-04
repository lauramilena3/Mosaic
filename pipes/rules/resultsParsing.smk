ruleorder: getAbundancesPE > getAbundancesSE

rule getAbundancesPE:
	input:
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
		unpaired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{{sampling}}.txt", sample=SAMPLES),
		paired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.{{sampling}}.txt", sample=SAMPLES),
	output:
		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	run:
		import pandas as pd
		import numpy as np

		lenght=7000
		percentage=0.7
		min_bases=5000
		for sampling in SAMPLING_TYPE:
			df_tpmean=pd.DataFrame()
			for sample in SAMPLES:
				#READ NUMBER
				paired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_paired_clean."+sampling+".txt")
				unpaired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_unpaired_clean."+sampling+".txt")
				paired=int(paired_size.readline())
				unpaired=int(unpaired_size.readline())
				reads=((paired*2)+unpaired)/1000000
				#NORMALIZE TP MEAN
				tpmean_file=dirs_dict["MAPPING_DIR"]+ "/" +sample+"_tpmean." + sampling + ".tsv"
				tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
				tpmean[sample] = tpmean[sample].apply(lambda x: x/reads)
				#REMOVE LOW COVERED CONTIGS
				breadth_file = dirs_dict["MAPPING_DIR"]+ "/" +sample+"_filtered_coverage." + sampling + ".txt"
				breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
				df=pd.merge(tpmean, breadth, on='contig', how='outer')
				#Divide dataframe in lenghts
				df['percentage']=df['breadth']/df['length']
				df=df.fillna(0)
				positive = df[(df['breadth']>7000) | (df['percentage']>percentage) ] 
				if df_tpmean.empty:
					positive.drop("breadth", axis=1, inplace=True)
					positive.drop("length", axis=1, inplace=True)
					positive.drop("percentage", axis=1, inplace=True)
					df_tpmean=positive
				else:
					positive.drop("length", axis=1, inplace=True)
					positive.drop("breadth", axis=1, inplace=True)
					positive.drop("percentage", axis=1, inplace=True)
					df_tpmean=pd.merge(positive, df_tpmean, on='contig', how='outer')
			filename="vOTU_abundance_table." + sampling + ".txt"
			df_tpmean=df_tpmean.fillna(0)
			df_tpmean.rename(columns={'contig':'#OTU ID'}, inplace=True)        
			df_tpmean.to_csv(dirs_dict["MAPPING_DIR"]+ "/" + filename, sep='\t', index=False, header=True)


rule getAbundancesSE:	
	input:
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
		unpaired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{{sampling}}.txt", sample=SAMPLES),
	output:
		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	run:
		import pandas as pd
		import numpy as np

		lenght=7000
		percentage=0.7
		min_bases=5000
		for sampling in SAMPLING_TYPE:
			df_tpmean=pd.DataFrame()
			for sample in SAMPLES:
				#READ NUMBER
				unpaired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_unpaired_clean."+sampling+".txt")
				unpaired=int(unpaired_size.readline())
				reads=(unpaired)/1000000
				#NORMALIZE TP MEAN
				tpmean_file=dirs_dict["MAPPING_DIR"]+ "/" +sample+"_tpmean." + sampling + ".tsv"
				tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
				tpmean[sample] = tpmean[sample].apply(lambda x: x/reads)
				#REMOVE LOW COVERED CONTIGS
				breadth_file = dirs_dict["MAPPING_DIR"]+ "/" +sample+"_filtered_coverage." + sampling + ".txt"
				breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
				df=pd.merge(tpmean, breadth, on='contig', how='outer')
				#Divide dataframe in lenghts
				df['percentage']=df['breadth']/df['length']
				df=df.fillna(0)
				positive = df[(df['breadth']>7000) | (df['percentage']>percentage) ] 
				if df_tpmean.empty:
					positive.drop("breadth", axis=1, inplace=True)
					positive.drop("length", axis=1, inplace=True)
					positive.drop("percentage", axis=1, inplace=True)
					df_tpmean=positive
				else:
					positive.drop("length", axis=1, inplace=True)
					positive.drop("breadth", axis=1, inplace=True)
					positive.drop("percentage", axis=1, inplace=True)
					df_tpmean=pd.merge(positive, df_tpmean, on='contig', how='outer')
			filename="vOTU_abundance_table." + sampling + ".txt"
			df_tpmean=df_tpmean.fillna(0)
			df_tpmean.rename(columns={'contig':'#OTU ID'}, inplace=True)        
			df_tpmean.to_csv(dirs_dict["MAPPING_DIR"]+ "/" + filename, sep='\t', index=False, header=True)

rule tabletoBIOM:
	input:
		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table.{sampling}.txt",
	output:
		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_json.{sampling}.biom",
	message:
		"Getting vOTU tables"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		biom convert -i {input.abundances} -o {output.abundances} --table-type="OTU table" --to-json		
		"""
rule getSummaryTable:
	input:
		hmm_results=dirs_dict["VIRAL_DIR"]+ "/hmm_parsed.{sampling}.out",
		table=dirs_dict["VIRAL_DIR"]+ "/viral_table.{sampling}.csv",
		high_dir=(dirs_dict["VIRAL_DIR"]+ "/high_confidence_vContact.{sampling}/genome_by_genome_overview.csv"),
		low_dir=(dirs_dict["VIRAL_DIR"]+ "/low_confidence_vContact.{sampling}/genome_by_genome_overview.csv"),
	output:
		summary=dirs_dict["MAPPING_DIR"]+ "/vOTU_summary.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	shell:
		"""
		touch {output.summary}
		"""