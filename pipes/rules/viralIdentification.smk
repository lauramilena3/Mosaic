rule virSorter:
	input:
		representatives=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		virSorter_dir=config['virSorter_dir'],
		virSorter_db=config['virSorter_db']
	output:
		results=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/VIRSorter_global-phage-signal.csv"
	params:
		out_folder=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 32
	shell:
		"""
		WD=$(pwd)
		{config[virSorter_dir]}/wrapper_phage_contigs_sorter_iPlant.pl -f {input.representatives} \
			--db 2 \
			--wdir {params.out_folder} \
			--ncpu {threads} \
			--data-dir $WD/{input.virSorter_db} \
			--virome
		"""

rule virFinder:
	input:
		representatives=dirs_dict["vOUT_DIR"] + "/representative_contigs.{sampling}.fasta",
		virFinder_dir=config['virFinder_dir'],
	output:
		pvalues=dirs_dict["VIRAL_DIR"] + "/virFinder_pvalues.{sampling}.txt"
	params:
		virFinder_script="scripts/virfinder_wrapper.R"
	message:
		"Scoring virus VirFinder"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		Rscript {params.virFinder_script} {input.representatives} {output.pvalues}
		"""

rule annotate_VIBRANT:
	input:
		representatives=dirs_dict["vOUT_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		VIBRANT_dir=os.path.join(workflow.basedir, config['vibrant_dir']),
	output:
		vibrant=directory(dirs_dict["VIRAL_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}"),
	params:
		viral_dir=directory(dirs_dict["VIRAL_DIR"]),
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	message:
		"Annotating contigs with VIBRANT"
	threads: 32
	shell:
		"""
		cd {params.viral_dir}
		{input.VIBRANT_dir}/VIBRANT_run.py -i {input.representatives} -t 1
		#{input.VIBRANT_dir}/VIBRANT_run.py -i {input.representatives} -virome -t 1
		"""
#		vibrant_figures=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_figures_" +REFERENCE_CONTIGS_BASE + ".tot"),
#		vibrant_tables_parsed=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_parsed_" +REFERENCE_CONTIGS_BASE + ".tot"),
#		vibrant_tables_unformated=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +REFERENCE_CONTIGS_BASE + ".tot"),
#		vibrant_phages=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +REFERENCE_CONTIGS_BASE + ".tot"),
#		vibrant_results=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +REFERENCE_CONTIGS_BASE + ".tot"),


rule parseViralTable:
	input:
		pvalues = dirs_dict["VIRAL_DIR"] + "/virFinder_pvalues.{sampling}.txt",
		categories=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/VIRSorter_global-phage-signal.csv",
		vibrant=(dirs_dict["VIRAL_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS),
	output:
		circular_H=dirs_dict["VIRAL_DIR"]+ "/high_confidence_circular_list.{sampling}.txt",
		circular_L=dirs_dict["VIRAL_DIR"]+ "/low_confidence_circular_list.{sampling}.txt",
		non_circular_H=dirs_dict["VIRAL_DIR"]+ "/high_confidence_non_circular_list.{sampling}.txt",
		non_circular_L=dirs_dict["VIRAL_DIR"]+ "/low_confidence_non_circular_list.{sampling}.txt",
		circular_unk=dirs_dict["VIRAL_DIR"]+ "/unknown_circular_list.{sampling}.txt",
		table=dirs_dict["VIRAL_DIR"]+ "/viral_table.{sampling}.csv"
	params:
		vibrant_results=dirs_dict["VIRAL_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS,
	message:
		"Parsing VirSorter and VirFinder results"
	threads: 1
	run:
		import pandas as pd
		results = pd.DataFrame(columns=('lenght', 'circular','type', 'VS_cat', 'VF_score', 'VF_pval'))


		#VirSorter
		VS = input.categories
		with open(VS) as fp:
			line = fp.readline()
			cnt = 1
			while line:
				if line.startswith("#"):
					if (line.strip().split()[1].isdigit()):
						contig_type=(line.strip().split("-")[1])
				else:
					circular="N"
					contigName=line.split(",")[0].split("VIRSorter_")[1].replace(".", "_")
					category=line.split(",")[4]
					if "-circular" in contigName:
						contigName=contigName.split("-circular")[0]
						circular="Y"
					if "suggestCircular=yes" in contigName:
						circular="Y"
					results.loc[contigName, 'VS_cat'] = int(category)
					results.loc[contigName, 'circular'] = circular
					results.loc[contigName, 'type'] = contig_type
				line = fp.readline()
				cnt += 1

		#VirFinder
		VF = input.pvalues
		with open(VF) as fp:
			line = fp.readline()
			cnt = 1
			while line:
				if cnt != 1:
					contigName=line.split("\t")[0].strip().replace(".", "_")
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

		df_A_nc=results[results['VS_cat']<3][results['circular']!="Y"]
		df_B_nc=results[results['VF_score']>0.9][results['VF_pval']<0.05][results['circular']!="N"]
		df_C_nc=results[results['VF_score']>0.7][results['VF_pval']<0.05][results['VS_cat']>0][results['circular']=="N"]

		df_circular=results[results['circular']=="Y"]


		#joinin and dereplicating circular contigs
		lsA=(df_A_c.index.tolist())
		lsB=(df_B_c.index.tolist())
		lsC=(df_C_c.index.tolist())

		lsAnotB=set(lsA) - set(lsB)
		lsAB=lsA + list(lsAnotB)

		lsCnoAB=set(lsC) - set(lsAB)
		lsCf=list(lsCnoAB)

		f=open(output.circular_H, 'w')
		f.write("\n".join(lsAB))
		f.close()

		f=open(output.circular_L, 'w')
		f.write("\n".join(lsCf))
		f.close()
		#dereplicating circular contigs

		lsCirc=df_circular.index.tolist()
		lsCircNothers=set(lsCirc) - set(lsA) - set(lsB) - set(lsC)

		f=open(output.circular_unk, 'w')
		f.write("\n".join(lsCircNothers))
		f.close()

		#joinin and dereplicating noncircular contigs
		lsA=(df_A_nc.index.tolist())
		lsB=(df_B_nc.index.tolist())
		lsC=(df_C_nc.index.tolist())

		lsAnotB=set(lsA) - set(lsB)
		lsAB=lsA + list(lsAnotB)

		lsCnoAB=set(lsC) - set(lsAB)
		lsCf=list(lsCnoAB)

		f=open(output.non_circular_H, 'w')
		f.write("\n".join(lsAB))
		f.close()

		f=open(output.non_circular_L, 'w')
		f.write("\n".join(lsCf))
		f.close()

		results.to_csv(output.table)

rule hmmCircularContigs:
	input:
		circular_unk=dirs_dict["VIRAL_DIR"]+ "/unknown_circular_list.{sampling}.txt",
		representatives=dirs_dict["vOUT_DIR"] + "/merged_scaffolds.{sampling}_95-80.fna",
	output:
		edited_fasta=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_95-80.{sampling}.fna",
		circular_unk_fasta=dirs_dict["VIRAL_DIR"]+ "/unknown_circular.{sampling}.fasta",
		hmm_results=dirs_dict["VIRAL_DIR"]+ "/hmm_parsed.{sampling}.out",
		hmm_out=dirs_dict["VIRAL_DIR"]+ "/hmmsearch.{sampling}.out",
		hmm_list=dirs_dict["VIRAL_DIR"]+ "/positive_rep_list.{sampling}.txt",
	params:
		hmm="db/hmm/ssDNA.hmm",
		min_score=50,
		min_eval=0.001
	message:
		"Selecting Viral Circular Contigs with hmmsearch"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		sed 's/\./_/g' {input.representatives} > {output.edited_fasta}
		seqtk subseq {output.edited_fasta} {input.circular_unk} > {output.circular_unk_fasta}
		if [ -s {output.circular_unk_fasta} ]
		then
			hmmsearch --tblout {output.hmm_out} -E {params.min_eval} {params.hmm} {output.circular_unk_fasta}
			cat {output.hmm_out} | grep -v '^#' | awk '{{ if ( $6 > {params.min_score} ) {{print $1,$3,$5,$6}} }}' > {output.hmm_results} || true
			cut -d' ' -f1 {output.hmm_results} | sort | uniq > {output.hmm_list}
		else
			touch {output.hmm_out}
			touch {output.hmm_results}
			touch {output.hmm_list}
		fi
		"""
rule extractViralContigs:
	input:
		circular_H=dirs_dict["VIRAL_DIR"]+ "/high_confidence_circular_list.{sampling}.txt",
		circular_L=dirs_dict["VIRAL_DIR"]+ "/low_confidence_circular_list.{sampling}.txt",
		non_circular_H=dirs_dict["VIRAL_DIR"]+ "/high_confidence_non_circular_list.{sampling}.txt",
		non_circular_L=dirs_dict["VIRAL_DIR"]+ "/low_confidence_non_circular_list.{sampling}.txt",
		hmm_list=dirs_dict["VIRAL_DIR"]+ "/positive_rep_list.{sampling}.txt",
		edited_fasta=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_95-80.{sampling}.fna"
	output:
		high_contigs_dup=dirs_dict["VIRAL_DIR"]+ "/high_confidence_dup.{sampling}.fasta",
		low_contigs_dup=dirs_dict["VIRAL_DIR"]+ "/low_confidence_dup.{sampling}.fasta",
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.fasta",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.fasta",
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/positive_contigs.{sampling}.fasta",
	message:
		"Selecting Viral Contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 1
	shell:
		"""
		#High Confidence
		#circular
		seqtk subseq {input.edited_fasta} {input.circular_H} > {output.high_contigs_dup}
		#non-circular
		seqtk subseq {input.edited_fasta} {input.non_circular_H} >> {output.high_contigs_dup}
		seqtk subseq {input.edited_fasta} {input.hmm_list} >> {output.high_contigs_dup}
		#Low Confidence
		#circular
		seqtk subseq {input.edited_fasta} {input.circular_L} > {output.low_contigs_dup}
		#non-circular
		seqtk subseq {input.edited_fasta} {input.non_circular_L} >> {output.low_contigs_dup}
		#filter duplicated sequences
		awk '/^>/{{f=!d[$1];d[$1]=1}}f' {output.high_contigs_dup} > {output.high_contigs}
		awk '/^>/{{f=!d[$1];d[$1]=1}}f' {output.low_contigs_dup} > {output.low_contigs}
		sed -i 's/_/-/g' {output.high_contigs}
		sed -i 's/_/-/g' {output.low_contigs}
		cat {output.high_contigs} {output.low_contigs} > {output.positive_contigs}
		"""
