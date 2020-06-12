ruleorder: getAbundancesPE > getAbundancesSE

rule getAbundancesPE:
	input:
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
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
		df_tpmean=pd.DataFrame()
		sampling=wildcards.sampling
		for sample in SAMPLES:
			#READ NUMBER
			paired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_paired_clean."+sampling+".txt")
			unpaired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_unpaired_clean."+sampling+".txt")
			paired=int(paired_size.readline())
			unpaired=int(unpaired_size.readline())
			reads=((paired*2)+unpaired)/1000000
			#NORMALIZE TP MEAN
			tpmean_file=dirs_dict["MAPPING_DIR"]+ "/BamM_" +sample+"_tpmean." + sampling + ".tsv"
			tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
			tpmean[sample] = tpmean[sample].apply(lambda x: x/reads)
			#REMOVE LOW COVERED CONTIGS
			breadth_file = dirs_dict["MAPPING_DIR"]+ "/bedtools_" +sample+"_filtered_coverage." + sampling + ".txt"
			breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
			df=pd.merge(tpmean, breadth, on='contig', how='outer')
			#Divide dataframe in lenghts
			df['percentage' ]=df['breadth']/df['length']
			df=df.fillna(0)
			positive = df[(df['breadth']>7000) | (df['percentage']>percentage) ]
			if df_tpmean.empty:
				positive.drop("breadth", axis=1, inplace=True)
				positive.drop("length", axis=1, inplace=True)
				#positive.drop("percentage", axis=1, inplace=True)
				df_tpmean=positive
			else:
				positive.drop("length", axis=1, inplace=True)
				positive.drop("breadth", axis=1, inplace=True)
				#positive.drop("percentage", axis=1, inplace=True)
				df_tpmean=pd.merge(positive, df_tpmean, on='contig', how='outer')
		filename="vOTU_abundance_table." + sampling + ".txt"
		df_tpmean=df_tpmean.fillna(0)
		df_tpmean.rename(columns={'contig':'#OTU ID'}, inplace=True)
		df_tpmean.to_csv(dirs_dict["MAPPING_DIR"]+ "/" + filename, sep='\t', index=False, header=True)

# rule getAbundancesPE_user:
# 	input:
# 		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
# 		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
# 		unpaired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{{sampling}}.txt", sample=SAMPLES),
# 		paired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.{{sampling}}.txt", sample=SAMPLES),
# 	output:
# 		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table.{sampling}.txt",
# 	message:
# 		"Getting vOTU tables"
# 	threads: 1
# 	run:
# 		import pandas as pd
# 		import numpy as np
#
# 		lenght=7000
# 		percentage=0.7
# 		min_bases=5000
# 		df_tpmean=pd.DataFrame()
# 		sampling=wildcards.sampling
# 		for sample in SAMPLES:
# 			#READ NUMBER
# 			paired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_paired_clean."+sampling+".txt")
# 			unpaired_size=open(dirs_dict["CLEAN_DATA_DIR"]+ "/" +sample+"_unpaired_clean."+sampling+".txt")
# 			paired=int(paired_size.readline())
# 			unpaired=int(unpaired_size.readline())
# 			reads=((paired*2)+unpaired)/1000000
# 			#NORMALIZE TP MEAN
# 			tpmean_file=dirs_dict["MAPPING_DIR"]+ "/" +sample+"_tpmean." + sampling + ".tsv"
# 			tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
# 			tpmean[sample] = tpmean[sample].apply(lambda x: x/reads)
# 			#REMOVE LOW COVERED CONTIGS
# 			breadth_file = dirs_dict["MAPPING_DIR"]+ "/" +sample+"_filtered_coverage." + sampling + ".txt"
# 			breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
# 			df=pd.merge(tpmean, breadth, on='contig', how='outer')
# 			#Divide dataframe in lenghts
# 			df['percentage']=df['breadth']/df['length']
# 			df=df.fillna(0)
# 			positive = df[(df['breadth']>7000) | (df['percentage']>percentage) ]
# 			if df_tpmean.empty:
# 				positive.drop("breadth", axis=1, inplace=True)
# 				positive.drop("length", axis=1, inplace=True)
# 				#positive.drop("percentage", axis=1, inplace=True)
# 				df_tpmean=positive
# 			else:
# 				positive.drop("length", axis=1, inplace=True)
# 				positive.drop("breadth", axis=1, inplace=True)
# 				#positive.drop("percentage", axis=1, inplace=True)
# 				df_tpmean=pd.merge(positive, df_tpmean, on='contig', how='outer')
# 		filename="vOTU_abundance_table." + sampling + ".txt"
# 		df_tpmean=df_tpmean.fillna(0)
# 		df_tpmean.rename(columns={'contig':'#OTU ID'}, inplace=True)
# 		df_tpmean.to_csv(dirs_dict["MAPPING_DIR"]+ "/" + filename, sep='\t', index=False, header=True)

rule getAbundancesSE:
	input:
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
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
		SAMPLING_TYPE={sampl}
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

rule getAbundancesDB:
	input:
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
		paired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.{{sampling}}.txt", sample=SAMPLES),
	output:
		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_DB.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	run:
		import pandas as pd
		import numpy as np
		import os
		os.chdir(RESULTS_DIR)
		sampling=wildcards.sampling
		df_tpmean=pd.DataFrame()
		for sample in SAMPLES:
		    #READ NUMBER
		    paired_size=open("02_CLEAN_DATA"+ "/" +sample+"_paired_clean."+sampling+".txt")
		    unpaired_size=open("02_CLEAN_DATA"+ "/" +sample+"_unpaired_clean."+sampling+".txt")
		    paired=int(paired_size.readline())
		    unpaired=int(unpaired_size.readline())
		    reads=((paired*2)+unpaired)/1000000
		    #NORMALIZE TP MEAN
		    tpmean_file="06_MAPPING"+ "/BamM_" +sample+"_tpmean." + sampling + ".tsv"
		    tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
		    tpmean[sample] = tpmean[sample].apply(lambda x: x/paired)
		    tpmean["contig"] = tpmean["contig"].str.strip()

		    breadth_file = "06_MAPPING"+ "/bedtools_" +sample+"_filtered_coverage." + sampling + ".txt"
		    with open(breadth_file) as f:
		        first_line = f.readline().strip()
		        if first_line=="":
		            splited=["", ""]
		        else:
		            splited=first_line.split(" ", 1)
		    breadth = pd.DataFrame([splited], columns=['breadth', 'contig'])
		    #print(tpmean['contig'].tolist())
		    #print(breadth['contig'].tolist())

		    df=pd.merge(tpmean,breadth,left_on='contig',right_on='contig')
		    #pd.merge(tpmean, breadth, on='contig', how='outer')
		    print(df)
		    print(sample)

		    df["breadth"] = pd.to_numeric(df["breadth"])
		    df['percentage' ]=df['breadth']/df['length']
		    df=df.fillna(0)
		    df.drop("breadth", axis=1, inplace=True)
		    df.drop("length", axis=1, inplace=True)
		    df.columns = ['contig', sample + "_depth", sample + "_breadth" ]
		    #REMOVE LOW COVERED CONTIGS
		    df=df[df[sample + "_breadth"]>0]
		    if df_tpmean.empty:
		        df_tpmean=df
		    else:
		        df_tpmean=pd.merge(df, df_tpmean, on='contig', how='outer')
		df_tpmean=df_tpmean.fillna(0)
		df_tpmean.rename(columns={'contig':'OTU'}, inplace=True)


		a_series = (df_tpmean != 0).any(axis=1)
		df_tpmean = df_tpmean.loc[a_series]
		df_tpmean.set_index('OTU', inplace=True)

		df_tpmean_70=df_tpmean.loc[(df_tpmean >= 0.6).any(axis=1)]

		filename="06_MAPPING/vOTU_abundance_table_DB." + sampling + ".txt"
		filename_70="06_MAPPING/vOTU_abundance_table_DB_70." + sampling + ".txt"

		df_tpmean.to_csv(filename, sep='\t', header=True)
		df_tpmean_70.to_csv(filename_70, sep='\t', header=True)

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
		genome_file=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/genome_by_genome_overview.csv",
	output:
		summary=dirs_dict["MAPPING_DIR"]+ "/vOTU_summary.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	shell:
		"""
		touch {output.summary}
		"""
