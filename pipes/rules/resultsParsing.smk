rule getAbundancesPE:
	input:
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_{{confidence}}_confidence_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_{{confidence}}_confidence_tpmean.{{sampling}}.tsv", sample=SAMPLES),
		unpaired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{{sampling}}.txt", sample=SAMPLES),
		paired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.{{sampling}}.txt", sample=SAMPLES),
	output:
		abundances=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence_vOTU_abundance_table.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	run:
		import pandas as pd
		import numpy as np

		lenght=7000
		percentage=0.7
		min_bases=5000
		SAMPLING=["tot", "sub"]
		CONFIDENCES=["high", "low"]
		for sampling in SAMPLING:
		    for confidence in CONFIDENCES:
		        df_tpmean=pd.DataFrame()
		        for sample in SAMPLES:
		            #READ NUMBER
		            paired_size=open(sample+"_paired_clean."+sampling+".txt")
		            unpaired_size=open(sample+"_unpaired_clean."+sampling+".txt")
		            paired=int(paired_size.readline())
		            unpaired=int(unpaired_size.readline())
		            reads=((paired*2)+unpaired)/1000000
		            #NORMALIZE TP MEAN
		            tpmean_file=sample+"_"+ confidence + "_confidence_tpmean." + sampling + ".tsv"
		            tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
		            tpmean[sample] = tpmean[sample].apply(lambda x: x/reads)
		            #REMOVE LOW COVERED CONTIGS
		            breadth_file = sample+"_"+ confidence + "_confidence_filtered_coverage." + sampling + ".txt"
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
		        filename=confidence+ "_confidence_vOTU_abundance_table." + sampling + ".txt"
		        df_tpmean=df_tpmean.fillna(0)
		        df_tpmean.to_csv(filename, sep='\t', index=True, header=False)
rule getAbundancesSE:	
	input:
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_{{confidence}}_confidence_filtered_coverage.{{sampling}}.txt", sample=SAMPLES),
		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/{sample}_{{confidence}}_confidence_tpmean.{{sampling}}.tsv", sample=SAMPLES),
		unpaired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.sub.txt", sample=SAMPLES),
	output:
		abundances=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence_vOTU_abundance_table.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	shell:
		"""
		touch {output.abundances}
		"""
rule tabletoBIOM:
	input:
		abundances=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence_vOTU_abundance_table.{sampling}.txt",
	output:
		abundances=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence_vOTU_abundance_table.{sampling}.biom",
	message:
		"Getting vOTU tables"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		touch {output.abundances}
		"""
rule getSummaryTable:
	input:
		hmm_results=dirs_dict["VIRAL_DIR"]+ "/hmm_parsed.{sampling}.out",
		pvalues = dirs_dict["VIRAL_DIR"] + "/virFinder_pvalues.{sampling}.txt",
		categories=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/VIRSorter_global-phage-signal.csv",
		high_dir=(dirs_dict["VIRAL_DIR"]+ "/high_confidence_vContact.{sampling}"),
		low_dir=(dirs_dict["VIRAL_DIR"]+ "/low_confidence_vContact.{sampling}"),
		representative_lenghts=dirs_dict["vOUT_DIR"] + "/representative_lengths.{sampling}.txt"
	output:
		summary=dirs_dict["MAPPING_DIR"]+ "/vOTU_summary.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	shell:
		"""
		touch {output.summary}
		"""