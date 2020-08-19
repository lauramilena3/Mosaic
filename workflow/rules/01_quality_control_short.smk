ruleorder: trim_adapters_quality_illumina_PE > trim_adapters_quality_illumina_SE
#ruleorder: listContaminants_PE > listContaminants_SE
#ruleorder: removeContaminants_PE > removeContaminants_SE
ruleorder: subsampleReadsIllumina_PE > subsampleReadsIllumina_SE
ruleorder: normalizeReads_PE > normalizeReads_SE
ruleorder: postQualityCheckIlluminaPE > postQualityCheckIlluminaSE

rule download_SRA:
	input:
		sratoolkit="tools/sratoolkit.2.10.0-ubuntu64"
	output:
		forward=(dirs_dict["RAW_DATA_DIR"] + "/{SRA}_pass_1.fastq"),
		reverse=(dirs_dict["RAW_DATA_DIR"] + "/{SRA}_pass_2.fastq"),
	params:
		SRA_dir=dirs_dict["RAW_DATA_DIR"],
	message:
		"Downloading SRA run"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
#	threads: 1
	shell:
		"""
		{input.sratoolkit}/bin/fastq-dump --outdir {params.SRA_dir} --skip-technical --readids --read-filter pass \\
		--dumpbase --split-files --clip -N 0 -M 0 {wildcards.SRA}
		"""

rule qualityCheckIllumina:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"]+"/{sample}_{reads}.fastq"
	output:
		html=temp(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{reads}_fastqc.html"),
		zipped=temp(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{reads}_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
#	threads: 1
	shell:
		"""
		fastqc {input}
		"""

rule multiQC:
	input:
		html=expand(dirs_dict["RAW_DATA_DIR"]+"/{sample}_{reads}_fastqc.html", sample=SAMPLES, reads=READ_TYPES),
		zipped=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{reads}_fastqc.zip", sample=SAMPLES, reads=READ_TYPES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
	params:
		fastqc_dir=dirs_dict["RAW_DATA_DIR"],
		html_name="preQC_illumina_report.html",
		multiqc_dir=dirs_dict["QC_DIR"],
	message:
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	shell:
		"""
		multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule trim_adapters_quality_illumina_PE:
	input:
		forward=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq",
		reverse=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq",
		qc_report=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html"
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
		forward_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq"),
		reverse_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq"),
		trimmomatic_values=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_trimmomatic_values.txt"),
	params:
		adapters=dirs_dict["ADAPTERS_DIR"] + "/" + config['adapters_file']
	message:
		"Trimming Illumina Adapters with Trimmomatic"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 8
	shell:
		"""
		echo leading {config[trimmomatic_leading]} trailing {config[trimmomatic_trailing]} winsize {config[trimmomatic_window_size]} winqual {config[trimmomatic_window_quality]} minlnth {config[trimmomatic_minlen]} > {output.trimmomatic_values}
		trimmomatic PE -threads {threads} -phred33 {input.forward} {input.reverse} \
		{output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
		ILLUMINACLIP:{params.adapters}:2:30:10:2:true LEADING:{config[trimmomatic_leading]} TRAILING:{config[trimmomatic_trailing]} \
		SLIDINGWINDOW:{config[trimmomatic_window_size]}:{config[trimmomatic_window_quality]} MINLEN:{config[trimmomatic_minlen]}
		"""

rule trim_adapters_quality_illumina_SE:
	input:
		forward=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq",
		qc_report=dirs_dict["QC_DIR"]+ "/pre_processing_multiqc_report.html"
	output:
		forward_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq"),
	params:
		adapters=dirs_dict["ADAPTERS_DIR"] + "/" + config['adapters_file']
	message:
		"Trimming Illumina Adapters with Trimmomatic"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 2
	shell:
		"""
		trimmomatic SE -threads {threads} -phred33 {input.forward} {output.forward_unpaired} \
		ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:{config[trimmomatic_leading]} TRAILING:{config[trimmomatic_trailing]} \
		SLIDINGWINDOW:{config[trimmomatic_window_size]}:{config[trimmomatic_window_quality]} MINLEN:{config[trimmomatic_minlen]}
		"""
# rule formatContaminants:
# 	input:
# 		contaminant_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.fasta",
# 	output:
# 		contaminant_bitmask=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.bitmask",
# 		contaminant_srprism=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.srprism.idx",
# 		contaminant_blastdb=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.fasta.nhr",
# 	message:
# 		"Creating contaminant databases"
# 	params:
# 		contaminants_dir=dirs_dict["CONTAMINANTS_DIR"],
# 		contaminant_srprism=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.srprism",
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	resources:
# 		mem_mb=32000
# 	shell:
# 		"""
# 		bmtool -d {input.contaminant_fasta} -o {output.contaminant_bitmask}  -w 18 -z
# 		srprism mkindex -i {input.contaminant_fasta} -o {params.contaminant_srprism} -M {resources.mem_mb}
# 		makeblastdb -in {input.contaminant_fasta} -dbtype nucl
# 		"""
#
# rule listContaminants_PE:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
# 		forward_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq"),
# 		reverse_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq"),
# 		contaminant_bitmask=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.bitmask",
# 		contaminant_srprism=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.srprism.idx",
# 		contaminant_blastdb=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.fasta.nhr",
# 	output:
# 		bmtagger_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-{contaminant}-BMTagger_paired.txt",
# 		bmtagger_unpaired_forward=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-{contaminant}-BMTagger_unpaired_forward.txt",
# 		bmtagger_unpaired_reverse=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-{contaminant}-BMTagger_unpaired_reverse.txt",
# 		temp_dir=temp(directory(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-{contaminant}_temp"))
# 	params:
# 		contaminant_srprism=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.srprism",
# 	message:
# 		"Removing contaminants with BMTagger"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		#PE
# 		#paired
# 		mkdir {output.temp_dir}
# 		bmtagger.sh -b {input.contaminant_bitmask} -x {params.contaminant_srprism} -T {output.temp_dir} -q 1 \
# 		-1 {input.forward_paired} -2 {input.reverse_paired} -o {output.bmtagger_paired}
# 		#unpaired
# 		bmtagger.sh -b {input.contaminant_bitmask} -x {params.contaminant_srprism} -T {output.temp_dir} -q 1 \
# 		-1 {input.forward_unpaired} -o {output.bmtagger_unpaired_forward}
# 		bmtagger.sh -b {input.contaminant_bitmask} -x {params.contaminant_srprism} -T {output.temp_dir} -q 1 \
# 		-1 {input.reverse_unpaired} -o {output.bmtagger_unpaired_reverse}
# 		"""
#
# rule removeContaminants_PE:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
# 		forward_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq"),
# 		reverse_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq"),
# 		bmtagger_paired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{{sample}}-{contaminant}-BMTagger_paired.txt", contaminant=CONTAMINANTS),
# 		bmtagger_unpaired_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{{sample}}-{contaminant}-BMTagger_unpaired_forward.txt", contaminant=CONTAMINANTS),
# 		bmtagger_unpaired_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{{sample}}-{contaminant}-BMTagger_unpaired_reverse.txt", contaminant=CONTAMINANTS),
# 	output:
# 		bmtagger_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-all-BMTagger_paired.txt"),
# 		bmtagger_unpaired_forward=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-all-BMTagger_unpaired_forward.txt"),
# 		bmtagger_unpaired_reverse=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-all-BMTagger_unpaired_reverse.txt"),
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.survived"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.survived"),
# 		forward_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq.survived",
# 		reverse_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq.survived",
# 	params:
# 		clean_data_dir=dirs_dict["CLEAN_DATA_DIR"],
# 	message:
# 		"Removing contaminants with BMTagger"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	threads: 1
# 	resources:
# 		mem_mb=48000
# 	shell:
# 		"""
# 		if [ -z {CONTAMINANTS} ]
# 		then
# 			ln -s {input.forward_paired} {output.forward_paired}
# 			ln -s {input.reverse_paired} {output.reverse_paired}
# 			ln -s {input.forward_unpaired} {output.forward_unpaired}
# 			ln -s {input.reverse_unpaired} {output.reverse_unpaired}
# 			touch {output.bmtagger_paired}
# 			touch {output.bmtagger_unpaired_forward}
# 			touch {output.bmtagger_unpaired_reverse}
# 		else
# 		#PE
# 		#paired
# 		cat {params.clean_data_dir}/{wildcards.sample}*BMTagger_paired.txt | sort | uniq > {output.bmtagger_paired}
# 		iu-remove-ids-from-fastq -i {input.forward_paired} -l {output.bmtagger_paired} -d " "
# 		iu-remove-ids-from-fastq -i {input.reverse_paired} -l {output.bmtagger_paired} -d " "
# 		#forward
# 		cat {params.clean_data_dir}/{wildcards.sample}*BMTagger_unpaired_forward.txt | sort | uniq > {output.bmtagger_unpaired_forward}
# 		iu-remove-ids-from-fastq -i {input.forward_unpaired} -l {output.bmtagger_unpaired_forward} -d " "
# 		#reverse
# 		cat {params.clean_data_dir}/{wildcards.sample}*BMTagger_unpaired_reverse.txt | sort | uniq > {output.bmtagger_unpaired_reverse}
# 		iu-remove-ids-from-fastq -i {input.reverse_unpaired} -l {output.bmtagger_unpaired_reverse} -d " "
# 		fi
# 		"""

rule remove_phiX174_PE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
		forward_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq",
		reverse_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq",
		phiX_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/GCF_000819615.1.fasta",
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq"),
		forward_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired_clean.tot.fastq"),
		reverse_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired_clean.tot.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq",
		paired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.tot.txt",
		unpaired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.txt",
	message:
		"Removing phiX174 with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 4
	resources:
		mem_mb=4000
	shell:
		"""
		#PE
		#PAIRED
		bbduk.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} ref={input.phiX_fasta} k=31 hdist=1 threads={threads}
		grep -c "^@" {output.forward_paired} > {output.paired_size}
		#UNPAIRED
		bbduk.sh -Xmx{resources.mem_mb}m in={input.forward_unpaired} out={output.forward_unpaired} ref={input.phiX_fasta} k=31 hdist=1 threads={threads}
		bbduk.sh -Xmx{resources.mem_mb}m in={input.reverse_unpaired} out={output.reverse_unpaired} ref={input.phiX_fasta} k=31 hdist=1 threads={threads}
		cat {output.forward_unpaired} {output.reverse_unpaired} > {output.unpaired}
		grep -c "^@" {output.unpaired} > {output.unpaired_size} ||  echo "0" > {output.unpaired_size}
		"""

rule postQualityCheckIlluminaPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq",
	output:
		html_forward=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_forward_paired_clean.tot_fastqc.html"),
		zipped_forward=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_forward_paired_clean.tot_fastqc.zip"),
		html_reverse=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_clean.tot_fastqc.html"),
		zipped_reverse=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_clean.tot_fastqc.zip"),
		html_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_unpaired_clean.tot_fastqc.html"),
		zipped_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_unpaired_clean.tot_fastqc.zip"),
	params:
		postQC_dir=dirs_dict["CLEAN_DATA_DIR"] +"/postQC",
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
#	threads: 1
	shell:
		"""
		fastqc {input.forward_paired} -o {params.postQC_dir}
		fastqc {input.reverse_paired} -o {params.postQC_dir}
		fastqc {input.unpaired} -o {params.postQC_dir}
		"""

rule postQualityCheckIlluminaSE:
	input:
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/post_{sample}_unpaired_clean.tot.fastq",
	output:
		html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/post_{sample}_unpaired_clean.tot_fastqc.html"),
		zipped=temp(dirs_dict["CLEAN_DATA_DIR"] + "/post_{sample}_unpaired_clean.tot_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
#	threads: 1
	shell:
		"""
		fastqc {input}
		"""

rule postMultiQC:
	input:
		html_forward=expand(dirs_dict["CLEAN_DATA_DIR"]  + "/postQC" + "/{sample}_forward_paired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_forward_paired_clean.tot_fastqc.zip", sample=SAMPLES),
		html_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_clean.tot_fastqc.zip", sample=SAMPLES),
		html_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_unpaired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"]  + "/postQC" + "/{sample}_unpaired_clean.tot_fastqc.zip", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/postQC_illumina_report.html"
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"] +  "/postQC",
		html_name="postQC_illumina_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	shell:
		"""
		multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
		"""
rule subsampleReadsIllumina_PE:
	input:
		unpaired_sizes=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.txt", sample=SAMPLES),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq"
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.sub.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.sub.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.sub.fastq",
		paired_size=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.sub.txt"),
		unpaired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.sub.txt"
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	params:
		max_subsample=int(int(config['max_subsample'])/2),
		files_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/*_unpaired_clean.tot.txt",
		files_paired=dirs_dict["CLEAN_DATA_DIR"] + "/*_paired_clean.tot.txt"
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		#paired
		paired=$( cat {params.files_paired} | sort -n | head -1 )
		p=$([ $paired -le {params.max_subsample} ] && echo "$paired" || echo {params.max_subsample})
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads=$p
		#unpaired
		unpaired_temp=$( cat {params.files_unpaired} | sort -n | head -1 )
		un=$([ $unpaired_temp -le {params.max_subsample} ] && echo "$unpaired_temp" || echo {params.max_subsample})
		reads_left=$(({params.max_subsample} - ($paired*2)))
		unpaired=$([ $un -le $reads_left ] && echo "$un" || echo $reads_left )
		reformat.sh in={input.unpaired} out={output.unpaired} reads=$unpaired
		grep -c "^@" {output.unpaired} > {output.unpaired_size}
		grep -c "^@" {output.forward_paired} > {output.paired_size}
		"""

rule subsampleReadsIllumina_SE:
	input:
		unpaired_sizes=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.txt", sample=SAMPLES),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq"
	output:
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.sub.fastq",
		unpaired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.sub.txt"
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	params:
		max_subsample=config['max_subsample'],
		files_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/*_unpaired_clean.tot.txt"
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		#unpairedte
		unpaired_temp=$( cat {params.files_unpaired} | sort -n | head -1 )
		echo $unpaired_temp
		un=$([ $unpaired_temp -le {params.max_subsample} ] && echo "$unpaired_temp" || echo {params.max_subsample})
		echo $un
		reformat.sh in={input.unpaired} out={output.unpaired} reads=$un
		grep -c "^@" {output.unpaired} > {output.unpaired_size}
		"""

rule normalizeReads_PE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq",
	message:
		"Normalizing reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	params:
		min_depth=config['min_norm'],
		max_depth=config['max_norm']
	threads: 4
	resources:
		mem_mb=6000
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
		target={params.max_depth} mindepth={params.min_depth} t={threads}
		#unpaired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth}
		"""

rule normalizeReads_SE:
	input:
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq"
	output:
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq"
	message:
		"Normalizing and error correcting reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	params:
		min_depth=config['min_norm'],
		max_depth=config['max_norm']
	threads: 4
	resources:
		mem_mb=8000
	shell:
		"""
		#SE
		#unpaired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth} t={threads}
		"""

rule kmer_rarefraction:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
	output:
		histogram=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{sampling}.csv"),
	message:
		"Counting unique reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 1
	shell:
		"""
		bbcountunique.sh in1={input.forward_paired} in2={input.reverse_paired} out={output.histogram}
		"""

rule plot_kmer:
	input:
		histograms=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{{sampling}}.csv", sample=SAMPLES),
	output:
		plot=(dirs_dict["CLEAN_DATA_DIR"] + "/kmer_rarefraction_plot.{sampling}.png"),
	message:
		"Plot unique reads with BBtools"
#	conda:
#		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 1
	run:
		import pandas as pd
		import seaborn as sns; sns.set()
		import matplotlib.pyplot as plt

		plt.figure(figsize=(12,12))
		sns.set(font_scale=2)
		sns.set_style("whitegrid")

		read_max=0

		for h in input.histograms:
		    df=pd.read_csv(h, sep="\t")
		    df.columns=["count", "percent", "c", "d", "e", "f", "g", "h", "i", "j"]
		    df=df[["count", "percent"]]
		    ax = sns.lineplot(x="count", y="percent", data=df,err_style='band', label=h.split("/")[-1].split("kmer")[0])
		    read_max=max(read_max,df["count"].max())

		ax.set(ylim=(0, 100))
		ax.set(xlim=(0, read_max*1.2))

		ax.set_xlabel("Read count",fontsize=20)
		ax.set_ylabel("New k-mers (%)",fontsize=20)
		ax.figure.savefig(output.plot)
