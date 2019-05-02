ruleorder: hybridAsemblySpades > shortReadAsemblySpadesPE > shortReadAsemblySpadesSE

rule hybridAsemblySpades:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{type}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{type}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{type}.fastq",
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{type}.fastq"
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{type}.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}/filtered_list.txt")
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}")
	message: 
		"Assembling hybrid reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		--meta -t {threads} --only-assembler --nanopore {input.nanopore}
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		"""
rule shortReadAsemblySpadesPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{type}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{type}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{type}.fastq",
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{type}.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}/filtered_list.txt")
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}")
	message: 
		"Assembling PE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		--meta -t {threads} --only-assembler
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		"""

rule shortReadAsemblySpadesSE:
	input:
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{type}.fastq"
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{type}.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}/filtered_list.txt")
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{type}")
	message: 
		"Assembling SE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		spades.py  --ss-1 {input.unpaired} -o {params.assembly_dir} \
		--meta -t {threads} --only-assembler
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		"""

rule downloadCanu:
	output:
		canu_dir=config['canu_dir'],
	message:
		"Downloading required VirSorter files"
	threads: 1
	shell:
		"""
		if [ ! -d {output.canu_dir} ]
		then
			if [ {config[operating_system]} == "macOs" ]
			then
				curl -OL https://github.com/marbl/canu/releases/download/v1.8/canu-1.8.Darwin-amd64.tar.xz
			else
				curl -OL https://github.com/marbl/canu/releases/download/v1.8/canu-1.8.Linux-amd64.tar.xz
			fi
		fi
		tar -xJf canu-1.8.*.tar.xz tools
		"""

rule asemblyCanu:
	input:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{type}.fastq",
		canu_dir=config['canu_dir']
	output:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu/{sample}.contigs.{type}.fasta"
	message:
		"Assembling nanopore reads with canu"
	params: 
		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu_{type}"
	threads: 4
	shell:
		"""
		./{config[canu_dir]}/canu genomeSize=5m minReadLength=1000 -p \
		contigFilter="{config[min_cov]} {config[min_len]} 1.0 1.0 2" \
		corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 \
		redMemory=32 oeaMemory=32 batMemory=200 -nanopore-raw {input.nanopore} \
		-d {params.assembly_dir} -p {wildcards.sample} useGrid=false executiveThreads={threads}
		"""

rule errorCorrectCanuPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{type}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{type}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{type}.fastq",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu_{type}/{sample}.contigs.fasta"
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu_filtered_scaffolds.{type}.fasta"),
		sam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired.{type}.sam",
		bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired.{type}.bam",
		sorted_bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{type}.bam",
		sorted_bam_paired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{type}.bam.bai",
		sam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{type}.sam",
		bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{type}.bam",
		sorted_bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{type}.bam",
		sorted_bam_unpaired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{type}.bam.bai"
	params:
		pilon_dir=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_{type}",
		scaffolds_pilon=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_{type}/pilon.fasta"),
		db_name=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_{type}"
	message:
		"Correcting nanopore assembly with pilon"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bowtie2-build -f {input.scaffolds} {params.db_name}
		#paired
		bowtie2 -x {params.db_name} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam_paired}
		samtools view -b -S {output.sam_paired} > {output.bam_paired}
		samtools sort {output.bam_paired} -o {output.sorted_bam_paired}
		samtools index {output.sorted_bam_paired} 
		#unpaired
		bowtie2 -x {params.db_name} -U {input.unpaired} -S {output.sam_unpaired}
		samtools view -b -S {output.sam_unpaired} > {output.bam_unpaired}
		samtools sort {output.bam_unpaired} -o {output.sorted_bam_unpaired}
		samtools index {output.sorted_bam_unpaired} 
		#PILON
		pilon --genome {input.scaffolds} --frags {output.sorted_bam_paired} --unpaired {output.sorted_bam_unpaired} \
		--outdir {params.pilon_dir}
		cp {params.scaffolds_pilon} {output.scaffolds}
		"""
rule errorCorrectCanuSE:
	input:
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{type}.fastq",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu_{type}/{sample}.contigs.fasta"
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu_filtered_scaffolds.{type}.fasta"),
		sam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{type}.sam",
		bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{type}.bam",
		sorted_bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{type}.bam",
		sorted_bam_unpaired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{type}.bam.bai"
	params:
		pilon_dir=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_{type}",
		scaffolds_pilon=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_{type}/pilon.fasta"),
		db_name=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_{type}"
	message:
		"Correcting nanopore assembly with pilon"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bowtie2-build -f {input.scaffolds} {params.db_name}
		#unpaired
		bowtie2 -x {params.db_name} -U {input.unpaired} -S {output.sam_unpaired}
		samtools view -b -S {output.sam_unpaired} > {output.bam_unpaired}
		samtools sort {output.bam_unpaired} -o {output.sorted_bam_unpaired}
		samtools index {output.sorted_bam_unpaired} 
		#PILON
		pilon --genome {input.scaffolds} --unpaired {output.sorted_bam_unpaired} --outdir {params.pilon_dir}
		cp {params.scaffolds_pilon} {output.scaffolds}
		"""

rule assemblyStats:
	input:
		scaffolds_canu=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu_filtered_scaffolds.{type}.fasta"),
		scaffolds_spades=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{type}.fasta")
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_quast_{type}"),
		quast_txt=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_quast_report.{type}.txt"
	message:
		"Creating assembly stats with quast"
	threads: 1
	shell:
		"""
		if [ -f {config[quast_dir]} ]
		then
			curl -OL https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
    		tar -xzf quast-5.0.2.tar.gz tools
    	fi
		./{config[quast_dir]}/quast.py {input.scaffolds_canu} {input.scaffolds_spades} -o {output.quast_report_dir}
		cp {output.quast_report_dir}/report.txt {output.quast_txt}
		"""
