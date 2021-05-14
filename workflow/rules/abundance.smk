ruleorder: mapReadsToContigsPE > mapReadsToContigsSE

rule createContigBowtieDb:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
		contigs_info=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta.fai",
		contigs_lenght=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_lenght.{sampling}.txt",
	params:
		prefix=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bowtie2-build -f {input.representatives} {params.prefix}
		#Get genome file
		samtools faidx {input.representatives}
		awk -F' ' '{{print $1"	"$2}}' {output.contigs_info} > {output.contigs_lenght}
		"""
# rule createContigBBDb:
# 	input:
# 		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
# 	output:
# 		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
# 		contigs_info=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta.fai",
# 		contigs_lenght=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_lenght.{sampling}.txt",
# 	params:
# 		prefix=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
# 	message:
# 		"Creating contig DB with Bowtie2"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		bowtie2-build -f {input.representatives} {params.prefix}
# 		#Get genome file
# 		samtools faidx {input.representatives}
# 		bbmap.sh ref={input.representatives}
# 		awk -F' ' '{{print $1"	"$2}}' {output.contigs_info} > {output.contigs_lenght}
#
# 		"""
rule mapReadsToContigsPE:
	input:
#		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		sam=dirs_dict["MAPPING_DIR"]+ "/bbmap_{sample}.{sampling}.sam",
		bam=dirs_dict["MAPPING_DIR"]+ "/bbmap_{sample}.{sampling}.bam",
 		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
 		bam_indexed=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam.bai",
		covstats=dirs_dict["MAPPING_DIR"]+ "/bbmap_covstats_{sample}.{sampling}.txt",
		covhist=dirs_dict["MAPPING_DIR"]+ "/bbmap_covhist_{sample}.{sampling}.txt",
		basecov=dirs_dict["MAPPING_DIR"]+ "/bbmap_basecov_{sample}.{sampling}.txt",
		bincov=dirs_dict["MAPPING_DIR"]+ "/bbmap_bincov_{sample}.{sampling}.txt",
		scafstats=dirs_dict["MAPPING_DIR"]+ "/bbmap_scafstats_{sample}.{sampling}.txt",
	params:
		contigs=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
 	resources:
 		mem_mb=12000
	shell:
		"""
		#bowtie2 --non-deterministic -x {params.contigs} -1 {input.forward_paired} \
		#-2 {input.reverse_paired} -U {input.unpaired} -S {output.sam} -p {threads}
		bbmap.sh -Xmx{resources.mem_mb}m ref={input.filtered_representatives} nodisk in1={input.forward_paired} in2={input.reverse_paired}  \
		outm={output.sam} threads={threads} covstats={output.covstats} covhist={output.covhist} basecov={output.basecov} \
		bincov={output.bincov} scafstats={output.scafstats}
		#Sam to Bam
		samtools view -b -S {output.sam} > {output.bam}
		samtools sort {output.bam} -o {output.bam_sorted}
		samtools index {output.bam_sorted}
		"""
rule mapReadsToContigsSE:
	input:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		sam=dirs_dict["MAPPING_DIR"]+ "/bbmap_{sample}.{sampling}.sam",
		bam=dirs_dict["MAPPING_DIR"]+ "/bbmap_{sample}.{sampling}.bam",
		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
		bam_indexed=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam.bai",
	params:
		contigs=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		#Mapping reads to contigs
		#bowtie2 --non-deterministic -x {params.contigs} -U {input.unpaired} -S {output.sam} -p {threads}
		#Sam to Bam
		samtools view -b -S {output.sam} > {output.bam}
		samtools sort {output.bam} -o {output.bam_sorted}
		samtools index {output.bam_sorted}
		"""

# rule filterBAM:
# 	input:
# 		bam=dirs_dict["MAPPING_DIR"]+ "/bowtie_{sample}.{sampling}.bam",
# 	output:
# 		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
# 		bam_filtered=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted_filtered.{sampling}.bam",
# 	params:
# 		out_dir=dirs_dict["MAPPING_DIR"],
# 		temp_bam_filtered=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}_filtered.bam",
# 		p_ident=config['p_ident'],
# 	message:
# 		"Filtering reads in Bam file with BamM"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env2.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		samtools sort {input.bam} -o {output.bam_sorted}
# 		bamm filter --bamfile {output.bam_sorted} --percentage_id {params.p_ident} -o {params.out_dir}
# 		mv {params.temp_bam_filtered} {output.bam_filtered}
# 		"""

rule tpmeanPerConfidence:
	input:
		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
		bam_indexed=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam.bai",
	output:
		tpmean=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{sampling}.tsv",
		counts=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_counts.{sampling}.tsv",
	message:
		"Calculating tpmean depth coverage"
	conda:
		dirs_dict["ENVS_DIR"] + "/env2.yaml"
	threads: 1
	shell:
		"""
		bamm parse -c {output.tpmean} -m tpmean -b {input.bam_sorted}
		bamm parse -c {output.counts} -m counts -b {input.bam_sorted}
		"""
rule getBreadthCoverage:
	input:
		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
	output:
		bam_cov=dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_genomecov.{sampling}.txt",
		cov_final=dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_coverage.{sampling}.txt",
		tmp_sort=temp(directory(dirs_dict["MAPPING_DIR"]+ "/temp_{sample}_{sampling}")),
	message:
		"Calculating breadth coverage contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bedtools genomecov -dz -ibam {input.bam_sorted} > {output.bam_cov}
		mkdir {output.tmp_sort}
		cut -f 1 {output.bam_cov} | sort -T {output.tmp_sort} | uniq -c | sort -nr -T {output.tmp_sort} | sed -e 's/^[[:space:]]*//' > {output.cov_final}
		"""
ruleorder: getAbundancesPE > getAbundancesSE

rule getAbundancesPE:
	input:
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_coverage.{{sampling}}.txt", sample=SAMPLES),
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
			reads=((paired))/1000000
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
		cov_final=expand(dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_coverage.{{sampling}}.txt", sample=SAMPLES),
		tpmean=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_tpmean.{{sampling}}.tsv", sample=SAMPLES),
		counts=expand(dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_counts.{{sampling}}.tsv", sample=SAMPLES),
		paired_size=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.{{sampling}}.txt", sample=SAMPLES),
	output:
		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_DB.{sampling}.txt",
		parsed_abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_DB_70.{sampling}.txt",
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
		df_counts=pd.DataFrame()

		for sample in SAMPLES:
		    #READ NUMBER
		    paired_size=open("02_CLEAN_DATA"+ "/" +sample+"_paired_clean."+sampling+".txt")
		    unpaired_size=open("02_CLEAN_DATA"+ "/" +sample+"_unpaired_clean."+sampling+".txt")
		    paired=int(paired_size.readline())
		    unpaired=int(unpaired_size.readline())
		    #reads=((paired*2)+unpaired)/1000000
		    #NORMALIZE TP MEAN
		    tpmean_file="06_MAPPING"+ "/BamM_" +sample+"_tpmean." + sampling + ".tsv"
		    tpmean = pd.read_csv(tpmean_file, sep="\t", header=0, names=("contig", "length", sample))
		    tpmean[sample] = tpmean[sample].apply(lambda x: x*1000000/paired)
		    tpmean["contig"] = tpmean["contig"].str.strip()
			#READ COUNTS
		    print(paired)
		    counts_file="06_MAPPING"+ "/BamM_" +sample+"_counts." + sampling + ".tsv"
		    counts = pd.read_csv(counts_file, sep="\t", header=0, names=("contig", "length", sample))
		    counts["contig"] = counts["contig"].str.strip()

		    breadth_file = "06_MAPPING"+ "/bedtools_" +sample+"_coverage." + sampling + ".txt"
		    #breadth = pd.read_csv(breadth_file, sep=" ", header=0, names=("breadth", "contig"))
		    contig=[]
		    brdth=[]
		    with open(breadth_file) as fp:
		        for line in fp:
		            first_line = line.strip()
		            if first_line=="":
		                contig.append("")
		                brdth.append("")
		            else:
		                brdth.append(first_line.split(" ", 1)[0])
		                contig.append(first_line.split(" ", 1)[1])
		        breadth = pd.DataFrame({'contig': contig,'breadth': brdth})
			# TPMEAN
		    df=pd.merge(tpmean,breadth,left_on='contig',right_on='contig')

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

			#COUNTS
		    df2=pd.merge(counts,breadth,left_on='contig',right_on='contig')

		    df2["breadth"] = pd.to_numeric(df2["breadth"])
		    df2['percentage']=df2['breadth']/df2['length']
		    df2=df2.fillna(0)
		    df2.drop("breadth", axis=1, inplace=True)
		    df2.drop("length", axis=1, inplace=True)
		    df2.columns = ['contig', sample + "_depth", sample + "_breadth" ]
		    #REMOVE NON COVERED CONTIGS
		    df2=df2[df2[sample + "_breadth"]>0]
		    if df_counts.empty:
		        df_counts=df2
		    else:
		        df_counts=pd.merge(df2, df_counts, on='contig', how='outer')
		#TP MEAN
		df_tpmean=df_tpmean.fillna(0)
		df_tpmean.rename(columns={'contig':'OTU'}, inplace=True)

		a_series = (df_tpmean != 0).any(axis=1)
		df_tpmean = df_tpmean.loc[a_series]
		df_tpmean.set_index('OTU', inplace=True)

		df_tpmean_70=df_tpmean.loc[(df_tpmean >= 0.7).any(axis=1)]
		cols_d = [c for c in df_tpmean_70.columns if c.lower()[-5:] == 'depth']
		cols_b = [c for c in df_tpmean_70.columns if c.lower()[-5:] != 'depth']

		df_tpmean_70_d=df_tpmean_70[cols_d]
		df_tpmean_70_b=df_tpmean_70[cols_b]

		filter_df=(df_tpmean_70_b >= 0.7)
		filter_df.columns=df_tpmean_70_d.columns
		filtered_df_tpmean_70_d=df_tpmean_70_d[filter_df].fillna(0)

		#COUNTS
		df_counts=df_counts.fillna(0)
		df_counts.rename(columns={'contig':'OTU'}, inplace=True)

		a_series = (df_counts != 0).any(axis=1)
		df_counts = df_counts.loc[a_series]
		df_counts.set_index('OTU', inplace=True)

		df_counts_70=df_counts.loc[(df_counts >= 0.7).any(axis=1)]
		cols_d = [c for c in df_counts_70.columns if c.lower()[-5:] == 'depth']
		cols_b = [c for c in df_counts_70.columns if c.lower()[-5:] != 'depth']

		df_counts_70_d=df_counts_70[cols_d]
		df_counts_70_b=df_counts_70[cols_b]

		filter_df2=(df_counts_70_b >= 0.7)
		filter_df2.columns=df_counts_70_d.columns
		filtered_df_counts_70_d=df_counts_70_d[filter_df2].fillna(0)

		#SAVE FILES

		filename="06_MAPPING/vOTU_abundance_table_DB." + sampling + ".txt"
		filename_70="06_MAPPING/vOTU_abundance_table_DB_70." + sampling + ".txt"
		filename_70_counts="06_MAPPING/vOTU_abundance_counts_DB_70." + sampling + ".txt"

		df_tpmean.to_csv(filename, sep='\t', header=True)
		filtered_df_tpmean_70_d.to_csv(filename_70, sep='\t', header=True)
		filtered_df_counts_70_d.to_csv(filename_70_counts, sep='\t', header=True)

rule tabletoBIOM:
	input:
		abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_DB.{sampling}.txt",
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
#		hmm_results=dirs_dict["VIRAL_DIR"]+ "/hmm_parsed.{sampling}.out",
#		table=dirs_dict["VIRAL_DIR"]+ "/viral_table.{sampling}.csv",
#		genome_file=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vContact.{sampling}/genome_by_genome_overview.csv",
	output:
		summary=dirs_dict["MAPPING_DIR"]+ "/vOTU_summary.{sampling}.txt",
	message:
		"Getting vOTU tables"
	threads: 1
	shell:
		"""
		touch {output.summary}
		"""
