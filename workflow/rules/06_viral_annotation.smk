rule annotate_VIGA:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		VIGA_dir=os.path.join(workflow.basedir, config['viga_dir']),
		piler_dir=os.path.join(workflow.basedir, (config['piler_dir'])),
		trf_dir=os.path.join(workflow.basedir, (config['trf_dir'])),
	output:
		modifiers=temp(dirs_dict["ANNOTATION"] + "/modifiers.txt"),
		temp_symlink=temp(dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta"),
		temp_viga_dir=temp(directory(dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_tempVIGA")),
		GenBank_file=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot" + "_annotated.gbk",
		GenBank_table_temp1=temp(dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot" + "_annotated.tbl"),
		GenBank_table_temp2=temp(dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot" + "_annotated.tbl2"),
		GenBank_table=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot" + ".tbl",
		GenBank_fasta=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot" + "_annotated.fasta",
		csv=dirs_dict["ANNOTATION"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot" + "_annotated.csv",
		viga_names=temp(dirs_dict["ANNOTATION"] + "/viga_names_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.txt"),
		viga_topology_temp=temp(dirs_dict["ANNOTATION"] + "/viga_topology_temp" + REPRESENTATIVE_CONTIGS_BASE + "tot.txt"),
		viga_topology=(dirs_dict["ANNOTATION"] + "/viga_topology_" + REPRESENTATIVE_CONTIGS_BASE + "tot.txt"),
	params:
		viga_log=dirs_dict["ANNOTATION"] + "/viga_log_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.txt",
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REPRESENTATIVE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
		VIGA_dir=directory("../" + config['viga_dir']),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	message:
		"Annotating contigs with VIGA"
	threads: 8
	shell:
		"""
		PILER={input.piler_dir}
		PATH=$PILER:$PATH
		TRF={input.trf_dir}
		PATH=$TRF:$PATH
		ln -sfn {input.representatives} {output.temp_symlink}
		mkdir -p {output.temp_viga_dir}
		cd {output.temp_viga_dir}
		touch {output.modifiers}
		echo "viga.1"
		{input.VIGA_dir}/VIGA.py --input {output.temp_symlink} --diamonddb {input.VIGA_dir}/databases/RefSeq_Viral_DIAMOND/refseq_viral_proteins.dmnd \
		--blastdb {input.VIGA_dir}/databases/RefSeq_Viral_BLAST/refseq_viral_proteins --hmmerdb {input.VIGA_dir}/databases/pvogs/pvogs.hmm \
		--rfamdb {input.VIGA_dir}/databases/rfam/Rfam.cm --modifiers {output.modifiers} --threads {threads} &> {params.viga_log}
		echo "viga.2"
		cat {params.viga_log} | grep "was renamed as" > {output.viga_names}
		echo "viga.3"
		cat {params.viga_log} | grep "according to LASTZ" > {output.viga_topology_temp}
		echo "viga.4"
		cat {output.viga_names} | while read line
		do
			stringarray=($line)
			new=${{stringarray[-1]}}
			old=${{stringarray[1]}}
			sed -i -e "s/${{new}}\t/${{old}}\t/g" -e "s/${{new}}_/${{old}}_/g" {output.csv}
			sed -i -e "s/${{new}}$/${{old}}/g" -e "s/${{new}} /${{old}} /g" -e "s/${{new}}_/${{old}}_/g" {output.GenBank_file}
			sed -i -e "s/${{new}}$/${{old}}/g" -e "s/${{new}} /${{old}} /g" -e "s/${{new}}_/${{old}}_/g" {output.GenBank_table_temp1}
			sed -i "s/>${{new}} $/>${{old}}/g" {output.GenBank_fasta}
			sed -i -e "s/${{new}} /${{old}} /g" {output.viga_topology_temp}
		done
		echo "viga.5"
		awk  '{{print $1 "\t" $6}}'  {output.viga_topology_temp} > {output.viga_topology}
		echo "viga.6"
		grep -v "gene$" {output.GenBank_table_temp1} > {output.GenBank_table_temp2}
		echo "viga.7"
		grep -n "CDS$" {output.GenBank_table_temp2} | cut -d : -f 1 | awk '{{$1+=-1}}1' | sed 's%$%d%' | sed -f - {output.GenBank_table_temp2} > {output.GenBank_table}
		echo "viga.8"
		sed -i "s/tRNA-?(Asp|Gly)(atcc)/tRNA-Xxx/g" {output.GenBank_table}
		echo "viga.9"
		"""
rule annotate_BLAST:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		blast=(os.path.join(workflow.basedir,"db/ncbi/NCBI_viral_proteins.faa")),
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + "_blast_viralRefSeq.{sampling}.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	message:
		"Annotating contigs with BLAST"
	threads: 8
	shell:
		"""
		blastp -num_threads {threads} -db {input.blast} -query {input.representatives} \
		-outfmt "6 qseqid sseqid stitle qstart qend qlen slen qcovs evalue length" > {output.blast_output}
		"""

rule blasToIMGVR:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		img_vr_db="/home/lmf/" + (config['img_vr_db']) + "IMGVR.fasta",
	output:
		blast_output=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + "_blast_output_IMG_VR.tot.csv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	message:
		"Blast contigs agaist IMG/VR database"
	threads: 8
	shell:
		"""
		blastn -num_threads {threads} -db {input.img_vr_db} -query {input.representatives} \
		-outfmt "6 qseqid sseqid qstart qend qlen slen qcovs evalue length" > {output.blast_output}
		"""

rule create_dbs_mmseqs2:
	input:
		MMseqs2_dir=(config['mmseqs_dir']),
		representatives=dirs_dict["vOUT_DIR"] + "/merged_scaffolds.tot_95-80.fna",
		reference=REPRESENTATIVE_CONTIGS
	output:
		index_representatives=dirs_dict["MMSEQS"] + "/representatives.index",
		index_reference=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".index",
		idx_reference=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".idx",
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REPRESENTATIVE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
	message:
		"Creating databases for reference and assembly mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 4
	shell:
		"""
		mmseqs createdb {input.representatives} {params.representatives_name}
		mmseqs createdb {input.reference} {params.reference_name}
		mmseqs createindex {params.reference_name} tmp --search-type 3
		"""
rule search_contigs_mmseqs2:
	input:
		index_representatives=dirs_dict["MMSEQS"] + "/representatives.index",
		index_reference=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".index",
		idx_reference=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".idx",
	output:
		results_index=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_search_results.index",
		results_table=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE + "_best_search_results.txt",
		temp_dir=temp(directory(dirs_dict["MMSEQS"] + "/tmp")),
	params:
		representatives_name=dirs_dict["MMSEQS"] + "/" + "representatives",
		reference_name=dirs_dict["MMSEQS"] + "/" + REPRESENTATIVE_CONTIGS_BASE,
		results_name=dirs_dict["MMSEQS"] + "/" +  REPRESENTATIVE_CONTIGS_BASE + "_search_results",
		mmseqs= "./" + config['mmseqs_dir'] + "/build/bin",
	message:
		"Comparing reference and assembly mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	threads: 16
	shell:
		"""
		mkdir {output.temp_dir}
		mmseqs search {params.representatives_name} {params.reference_name} {params.results_name} {output.temp_dir} \
		--start-sens 1 --sens-steps 3 -s 7 --search-type 3 --threads {threads}
		mmseqs convertalis {params.representatives_name} {params.reference_name} {params.results_name} {output.results_table}
		"""

rule create_WIsH_models:
	input:
		wish_dir=os.path.join(workflow.basedir, (config['wish_dir'])),
		FNA=directory("db/PATRIC/FNA"),
	output:
		model_dir=directory("db/PATRIC/FNA/wish_modelDir"),
	params:
		model_dir_ln="db/PATRIC/FNA/wish_modelDir/wish_modelDir_ln",
	message:
		"Create WIsH bacterial DB"
	threads: 1
	shell:
		"""
		{input.wish_dir}/WIsH -c build -g {input.FNA} -m {output.model_dir}
		mkdir {params.model_dir_ln}
		cd {params.model_dir_ln}
		ln -s ../*.mm .
		"""

rule hostID_WIsH:
	input:
		wish_dir=os.path.join(workflow.basedir, (config['wish_dir'])),
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		model_dir=("db/PATRIC/FNA/wish_modelDir"),
	output:
		results_dir=directory(dirs_dict["VIRAL_DIR"] + "/wish/wish_" + REPRESENTATIVE_CONTIGS_BASE + "_resultsDir"),
		phages_dir=directory(dirs_dict["VIRAL_DIR"] + "/wish/wish_" + REPRESENTATIVE_CONTIGS_BASE + "_phagesDir"),
	params:
		model_dir_ln="db/PATRIC/FNA/wish_modelDir/wish_modelDir_ln",
		phages_dir=dirs_dict["VIRAL_DIR"] + "/wish_modelDir_ln",
	message:
		"Host finding with WIsH"
	threads: 1
	shell:
		"""
		mkdir {output.phages_dir}
		cd {output.phages_dir}
		awk -F '>' '/^>/ {{F=sprintf("%s.fa", $2); print > F;next;}} {{print F; close(F)}}' < {input.representatives}
		cd {workflow.basedir}
		mkdir {output.results_dir}
		{input.wish_dir}/WIsH -c predict -g {output.phages_dir} -m {params.model_dir_ln} -r {output.results_dir} -b
		"""

rule mapReadstoContigsPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{contigs}.fasta",
	output:
		sam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired.{sampling}_to_{contigs}.sam",
		bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired.{sampling}_to_{contigs}.bam",
		sorted_bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{sampling}_to_{contigs}.bam",
		sorted_bam_paired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{sampling}_to_{contigs}.bam.bai",
		sam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{sampling}_to_{contigs}.sam",
		bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{sampling}_to_{contigs}.bam",
		sorted_bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{sampling}_to_{contigs}.bam",
		sorted_bam_unpaired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{sampling}_to_{contigs}.bam.bai",
	params:
		db_name=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_{sampling}_to_{contigs}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		bowtie2-build -f {input.scaffolds} {params.db_name} --threads {threads}
		#paired
		bowtie2 -x {params.db_name} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam_paired} --threads {threads}
		samtools view -b -S {output.sam_paired} > {output.bam_paired}
		samtools sort {output.bam_paired} -o {output.sorted_bam_paired}
		samtools index {output.sorted_bam_paired}
		#unpaired
		bowtie2 -x {params.db_name} -U {input.unpaired} -S {output.sam_unpaired} --threads {threads}
		samtools view -b -S {output.sam_unpaired} > {output.bam_unpaired}
		samtools sort {output.bam_unpaired} -o {output.sorted_bam_unpaired}
		samtools index {output.sorted_bam_unpaired}
		"""
rule annotate_VIBRANT:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		VIBRANT_dir=os.path.join(workflow.basedir, config['vibrant_dir']),
	output:
		vibrant=directory(dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + ".{sampling}"),
		vibrant_circular=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_circular.{sampling}.txt",
		vibrant_positive=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_list.{sampling}.txt",
		vibrant_quality=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + "_positive_quality.{sampling}.txt",
	params:
		viral_dir=directory(dirs_dict["vOUT_DIR"]),
		minlen=5000,
		name_circular=dirs_dict["vOUT_DIR"] + "/VIBRANT_" + REPRESENTATIVE_CONTIGS_BASE  + ".{sampling}/VIBRANT_results*/VIBRANT_complete_circular*.{sampling}.tsv"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	message:
		"Annotating viral contigs with VIBRANT"
	threads: 8
	shell:
		"""
		cd {params.viral_dir}
		{input.VIBRANT_dir}/VIBRANT_run.py -i {input.representatives} -t {threads} -virome
		cut -f1 {params.name_circular} > {output.vibrant_circular}
		touch {output.vibrant_circular}
		cp {output.vibrant}/VIBRANT_phages_*/*phages_combined.txt {output.vibrant_positive}
		cp {output.vibrant}/VIBRANT_results*/VIBRANT_genome_quality*.tsv {output.vibrant_quality}
		#vibrant_figures=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_figures_" +VIRAL_CONTIGS_BASE + ".tot"),
		#vibrant_tables_parsed=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_parsed_" +VIRAL_CONTIGS_BASE + ".tot"),
		#vibrant_tables_unformated=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +VIRAL_CONTIGS_BASE + ".tot"),
		#vibrant_phages=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +VIRAL_CONTIGS_BASE + ".tot"),
		#vibrant_results=(directory(dirs_dict["ANNOTATION"] + "/VIBRANT_HMM_tables_unformatted_" +VIRAL_CONTIGS_BASE + ".tot"),
		"""

rule detectNucleotideModifications:
	input:
		fastq_file=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore.fastq",
		fast5_dir=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore_single_fast5",
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
	output:
		plus_wig=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + ".fraction_modified_reads.plus.wig"),
		minus_wig=(dirs_dict["ANNOTATION"] + "/"+ REPRESENTATIVE_CONTIGS_BASE + ".fraction_modified_reads.minus.wig"),
	params:
		representative_basename=REPRESENTATIVE_CONTIGS_BASE,
	message:
		"Detecting nucleotide modifications with tombo"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 16
	shell:
		"""
		tombo preprocess annotate_raw_with_fastqs --fast5-basedir {input.fast5_dir} --fastq-filenames {input.fastq_file} --overwrite --processes {theads}
		tombo resquiggle {input.fast5_dir} {input.representatives} --processes {threads}
		tombo detect_modifications de_novo --fast5-basedirs {input.fast5_dir} --statistics-file-basename {params.representative_basename}.de_novo --processes {threads}
		tombo text_output browser_files --fast5-basedirs {input.fast5_dir} --statistics-filename {params.representative_basename}.de_novo.tombo.stats \
		--genome-fasta {input.representatives} --browser-file-basename {params.representative_basename} --file-types fraction
		"""
rule estimateGenomeCompletness:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		quality_summary=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/quality_summary.tsv",
		completeness=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/completeness.tsv",
		contamination=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/contamination.tsv",
		repeats=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/repeats.tsv",
	params:
		checkv_outdir=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}",
		checkv_db=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}",
	message:
		"Estimating genome completeness with CheckV "
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	threads: 4
	shell:
		"""
		checkv contamination {input.representatives} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv completeness {input.representatives} {params.checkv_outdir} -t {threads} -d {config[checkv_db]}
		checkv repeats {input.representatives} {params.checkv_outdir}
		checkv quality_summary {input.representatives} {params.checkv_outdir}
		"""

rule parseSummary:
	input:
		quality_summary=dirs_dict["vOUT_DIR"] + "/checkV_{sampling}/quality_summary.tsv",
		viral_boundary=dirs_dict["VIRAL_DIR"] + "/virSorter_{sampling}/final-viral-boundary.tsv",
		tsv=(dirs_dict["vOUT_DIR"] + "/taxonomy_report_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.tsv"),
		taxonomy_results=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + "_vcontact2_taxonomy.{sampling}.csv",
		parsed_abundances=dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_DB_70.{sampling}.txt",
	output:
		summary=dirs_dict["ANNOTATION"] + "/summary_information.{sampling}.csv",
	message:
		"Assigning viral taxonomy with vContact2 results"
	threads: 1
	run:
		import pandas as pd
		df1=pd.read_csv(input.quality_summary, sep="\t")
		df1=df1[["contig_id","contig_length","gene_count","viral_genes","host_genes","checkv_quality","provirus","termini"]]
		df2=pd.read_csv(input.viral_boundary, sep="\t")
		df2=df2[["seqname","group"]]
		df3=pd.read_csv(input.tsv, sep="\t",header=None, names=["name", "id", "rank", "taxonomy_mmseqs"])
		df3=df3[["name","rank", "taxonomy_mmseqs"]]
		df4=pd.read_csv(input.taxonomy_results, sep="\t",header=None, names=["name", "taxonomy_vcontact2"])
		df5=pd.read_csv(input.parsed_abundances,sep="\t")
		df1.merge(df2, left_on='contig_id', right_on='seqname', how="outer").merge(df3, left_on='contig_id', right_on='name', how="outer").merge(df4, left_on='contig_id', right_on='name', how="outer").merge(df5, left_on='contig_id', right_on='OTU', how="outer").to_csv(output.summary)
