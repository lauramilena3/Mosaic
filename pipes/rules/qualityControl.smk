rule qualityCheckIllumina:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"]+"/{sample}_{type}.fastq"
	output:
		html=dirs_dict["RAW_DATA_DIR"] + "/{sample}_{type}_fastqc.html",
		zipped=(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{type}_fastqc.zip")
	message: 
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
#	threads: 1
	shell:
		"""
		fastqc {input}
		"""
rule qualityCheckNanopore:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"]+"/{sample}_nanopore.fastq"
	output:
		nanoqc_dir=dirs_dict["RAW_DATA_DIR"] + "/{sample}_nanoplot",
		nanoqc=dirs_dict["QC_DIR"] + "/{sample}_nanopore.html"
	message: 
		"Performing nanoQC statistics"
	params:
		nanopore="nanopore"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
#	threads: 1
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.raw_fastq}
		mv {output.nanoqc_dir}/summary.html {output.nanoqc}
		"""
rule multiQC:
	input:
		forward=expand(dirs_dict["RAW_DATA_DIR"]+"/{sample}_{type}_fastqc.html", sample=SAMPLES, type=READ_TYPES)
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/pre_processing_multiqc_report.html"
	params:
		fastqc_dir=dirs_dict["RAW_DATA_DIR"],
		html_name="pre_processing_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message: 
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	shell:
		"""
		multiqc {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
		"""
rule trim_adapters_quality_illumina:
	input:
		forward=dirs_dict["RAW_DATA_DIR"] + "/{sample}_R1.fastq",
		reverse=dirs_dict["RAW_DATA_DIR"] + "/{sample}_R2.fastq",
		qc_report=dirs_dict["QC_DIR"]+ "/pre_processing_multiqc_report.html"
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
		forward_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq"),
		reverse_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq")
	params:
		adapters=dirs_dict["ADAPTERS_DIR"] + "/" + config['adapters_file']
	message: 
		"Trimming Illumina Adapters with Trimmomatic"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 2
	shell:
		"""
		trimmomatic PE -threads {threads} -phred33 {input.forward} {input.reverse} \
		{output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
		ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:{config[trimmomatic_leading]} TRAILING:{config[trimmomatic_trailing]} \
		SLIDINGWINDOW:{config[trimmomatic_window_size]}:{config[trimmomatic_window_quality]} MINLEN:{config[trimmomatic_minlen]}
		"""

rule remove_adapters_quality_nanopore:
	input:
		raw_data=dirs_dict["RAW_DATA_DIR"] + "/{sample}_nanopore.fastq",
		qc_report=dirs_dict["QC_DIR"]+ "/pre_processing_multiqc_report.html"
	output:
		(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore.fastq")
	message: 
		"Trimming Nanopore Adapters with Porechop"
	params:
		porechopped=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_porechopped.fastq"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env3.yaml"
	threads: 2
	shell:
		"""
		porechop -i {input.raw_data} -o {params.porechopped} --threads {threads}
		NanoFilt -q 10 -l 1000 --headcrop 50 {params.porechopped} > {output}
		"""

rule getContaminants:
	output:
		contaminants_file=dirs_dict["CONTAMINANTS_DIR"] +"/contaminants.fasta"
	message: 
		"Downloading contaminant genomes"
	params:
		contaminants_dir=dirs_dict["CONTAMINANTS_DIR"]
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	shell:
		"""
		line=$(echo {CONTAMINANTS})
        for contaminant in $line
        do
        	echo $contaminant
			wget $(esearch -db "assembly" -query $contaminant | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{{print $0"/"$NF"_genomic.fna.gz"}}')
			gunzip -f *$contaminant*gz
			cat *$contaminant*fna >> {output.contaminants_file}
			rm *$contaminant*fna
        done
		"""

rule removeContaminants:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
		forward_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq"),
		reverse_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq"),
		contaminants_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/contaminants.fasta"
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.fastq",
		singletons=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_singletons.fastq",
		temp_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_temp_unpaired.fastq"
	message: 
		"Removing contaminants with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		#PE
		#paired
		bbduk.sh -Xmx{resources.mem_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
		outs={output.singletons} ref={input.contaminants_fasta} k=31 hdist=1 threads={threads}
		#unpaired
		cat {input.forward_unpaired} {input.reverse_unpaired} > {output.temp_unpaired}
		bbduk.sh -Xmx{resources.mem_mb}m in={output.temp_unpaired} out={output.unpaired} ref={input.contaminants_fasta} k=31 hdist=1 threads={threads} 
		#singletons
		cat {output.singletons} >> {output.unpaired}
		"""
rule normalizeReads:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.fastq"
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.fastq"
	message: 
		"Normalizing reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	params: 
		min_depth=config['min_norm'],
		max_depth=config['max_norm']
	threads: 1
	resources:
		mem_mb=4000
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
		target={params.max_depth} mindepth={params.min_depth}
		#unpaired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth}		
		"""

