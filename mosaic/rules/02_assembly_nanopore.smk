#ruleorder: asemblyCanuPOOLED > asemblyCanu
ruleorder: hybridAsemblySpades > shortReadAsemblySpadesPE
#ruleorder: errorCorrectPE > errorCorrectSE
# ruleorder: assemblyStatsHYBRID > assemblyStatsILLUMINA
#ruleorder: mergeAssembliesHYBRID > mergeAssembliesSHORT


if POOLED==True:
	#ruleorder: symlinkPooled>subsampleReadsNanopore
	#ruleorder: symlinkPooled>remove_adapters_quality_nanopore
	ruleorder: hybridAsemblySpadesPooled>hybridAsemblySpades>shortReadAsemblySpadesPE
	rule hybridAsemblySpadesPooled:
		benchmark:
			dirs_dict["BENCHMARKS"] + "/hybridAsemblySpadesPooled/sample={sample}__sampling={sampling}.tsv"
		input:
			forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq.gz"),
			reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq.gz"),
			unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq.gz",
			nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/"+ NANOPORE_SAMPLES +"_nanopore_clean.{sampling}.fastq.gz"
		output:
			scaffolds=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.unrenamed.fasta"),
			filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/filtered_list.txt")
		params:
			raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/scaffolds.fasta",
			assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}"),
			metagenomic_flag=METAGENOME_FLAG,
			mem_gb=HYBRID_MEM_GB,
		message:
			"Assembling hybrid reads with metaSpades"
		conda:
			dirs_dict["ENVS_DIR"] + "/env3.yaml"
		threads: 16
		resources:
			mem_mb=HYBRID_MEM_MB
		shell:
			"""
			spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
			{params.metagenomic_flag}  -t {threads} --nanopore {input.nanopore} --memory {params.mem_gb}
			grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
			| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
			seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
			"""

rule hybridAsemblySpades:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq.gz",
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq.gz"
	output:
		scaffolds=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.unrenamed.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/filtered_list.txt")
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}"),
		metagenomic_flag=METAGENOME_FLAG,
		mem_gb=HYBRID_MEM_GB,
	message:
		"Assembling hybrid reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/hybridAsemblySpades/sample={sample}__sampling={sampling}.tsv"
	threads: 16
	resources:
		mem_mb=HYBRID_MEM_MB
	shell:
		"""
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		{params.metagenomic_flag}  -t {threads} --nanopore {input.nanopore} --memory {params.mem_gb}
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds}
		"""

rule asemblyCanu:
	input:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq.gz",
		#canu_dir=config['canu_dir'],
	output:
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/canu_{sample}_{sampling}"),
		scaffolds_final=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_canu.{sampling}.unrenamed.fasta"),
	message:
		"Assembling Nanopore reads with Canu"
	params:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/canu_{sample}_{sampling}/{sample}.contigs.fasta",
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_longread.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/asemblyCanu/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		# ./{config[canu_dir]}/canu genomeSize=100k minReadLength=1000 -p \
		# contigFilter="{config[min_cov]} {config[min_len]} 1.0 1.0 2" \
		# corOutCoverage=all corMhapSensitivity=high correctedErrorRate=0.105 corMinCoverage=0 \
		# corMaxEvidenceCoverageLocal=10 corMaxEvidenceCoverageGlobal=10 \
		# redMemory=32 oeaMemory=32 batMemory=200 -nanopore {input.nanopore} \
		# -d {output.assembly_dir} -p {wildcards.sample} useGrid=false maxThreads={threads}

		canu genomeSize=5m minReadLength=1000 -fast\
		contigFilter="{config[min_cov]} {config[min_len]} 1.0 1.0 2" \
		corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 \
		redMemory=32 oeaMemory=32 batMemory=200 -nanopore {input.nanopore} \
		-d {output.assembly_dir} -p {wildcards.sample} useGrid=false maxThreads={threads}
		cp {params.scaffolds} {output.scaffolds_final}
		sed -i s"/ /_/"g {output.scaffolds_final}
		"""

rule asemblyFlye:
	input:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq.gz",
	output:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample}_{sampling}/assembly.fasta",
		scaffolds_final=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_flye.{sampling}.unrenamed.fasta")
	message:
		"Assembling Nanopore reads with Flye"
	params:
		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample}_{sampling}",
		genome_size=config["genome_size"],
		metagenomic_flag=METAGENOME_FLAG,
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_longread.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/asemblyFlye/sample={sample}__sampling={sampling}.tsv"
	threads: 32
	shell:
		"""
		rm -rf {params.assembly_dir}
		flye --nano-hq {input.nanopore} --out-dir {params.assembly_dir} --genome-size {params.genome_size} --threads {threads} {params.metagenomic_flag}
		cp {output.scaffolds} {output.scaffolds_final}
		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds_final}
		"""

rule errorCorrectMedaka:
	input:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
		#corrected4=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_4_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq.gz",
	output:
		corrected_medaka=temp(dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.unrenamed.fasta"),
		#temp=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/medaka_temp_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}")),
		fai=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.fasta.fai"),
		mmi=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.fasta.map-ont.mmi"),
		medaka_dir1=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}")),
		medaka_dir2=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}")),
	params:
		medaka_model=config["medaka_model"], # "r941_min_high_g360"
		medaka_fasta1=dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}/consensus.fasta",
		medaka_fasta2=dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}/consensus.fasta",
	message:
		"Correcting nanopore assembly with long reads using two rounds of Medaka"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/errorCorrectMedaka/sample={sample}__sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		medaka_consensus -i {input.nanopore} -d {input.scaffolds} -o {output.medaka_dir1} -m {params.medaka_model} -f
		medaka_consensus -i {input.nanopore} -d {params.medaka_fasta1} -o {output.medaka_dir2} -m {params.medaka_model} -f
		cp {params.medaka_fasta2} {output.corrected_medaka}
		sed "s/>/>medaka_/g" -i {output.corrected_medaka}
		"""
		
rule errorCorrectRacon_2rounds:
	input:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq.gz",
		corrected_medaka=dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
	output:
		overlap1=dirs_dict["ASSEMBLY_DIR"] + "/minimap2_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}.paf",
		corrected1=temp(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}.unrenamed.fasta"),
		overlap2=dirs_dict["ASSEMBLY_DIR"] + "/minimap2_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.paf",
		corrected2=temp(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.unrenamed.fasta"),
	message:
		"Correcting nanopore assembly with long reads using four rounds of Racon "
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_longread.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/errorCorrectRacon_2rounds/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		#Racon round 1
		minimap2 -t {threads} {input.corrected_medaka} {input.nanopore} > {output.overlap1}
		racon -m 8 -x -6 -g -8 -w 500 -t {threads}  {input.nanopore} {output.overlap1} {input.corrected_medaka} > {output.corrected1}
		sed "s/>/>racon_1_/g" -i {output.corrected1}

		#Racon round 2
		minimap2 -t {threads} {output.corrected1} {input.nanopore} > {output.overlap2}
		racon -m 8 -x -6 -g -8 -w 500 -t {threads}  {input.nanopore} {output.overlap2} {output.corrected1} > {output.corrected2}
		sed "s/>/>racon_2_/g" -i {output.corrected2}
		"""


rule errorCorrectPilonPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq.gz",
		corrected2_racon=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
	output:
		bowtie_index=temp(expand(dirs_dict["ASSEMBLY_DIR"] + "/{{sample}}_bowtieDB_{round}_{{sampling}}.{part}.bt2", round=[1, 2, 3, 4], part=["1", "2", "3", "4", "rev.1", "rev.2"])),
		#round1
		sam_paired1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_1.{sampling}.sam",
		bam_paired1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_1.{sampling}.bam",
		sorted_bam_paired1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted_1.{sampling}.bam",
		sorted_bam_paired_ix1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted_1.{sampling}.bam.bai",
		sam_unpaired1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_1.{sampling}.sam",
		bam_unpaired1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_1.{sampling}.bam",
		sorted_bam_unpaired1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted_1.{sampling}.bam",
		sorted_bam_unpaired_ix1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted_1.{sampling}.bam.bai",
		scaffolds_pilon1=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_1_{sampling}/pilon.fasta"),
		scaffolds_pilon1_final=temp(dirs_dict["ASSEMBLY_DIR"] + "/pilon_1_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.unrenamed.fasta"),
		#round2
		sam_paired2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_2.{sampling}.sam",
		bam_paired2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_2.{sampling}.bam",
		sorted_bam_paired2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted_2.{sampling}.bam",
		sorted_bam_paired_ix2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted_2.{sampling}.bam.bai",
		sam_unpaired2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_2.{sampling}.sam",
		bam_unpaired2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_2.{sampling}.bam",
		sorted_bam_unpaired2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted_2.{sampling}.bam",
		sorted_bam_unpaired_ix2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted_2.{sampling}.bam.bai",
		scaffolds_pilon2=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_2_{sampling}/pilon.fasta"),
		scaffolds_pilon2_final=temp(dirs_dict["ASSEMBLY_DIR"] + "/pilon_2_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.unrenamed.fasta"),
		#round3
		sam_paired3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_3.{sampling}.sam",
		bam_paired3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_3.{sampling}.bam",
		sorted_bam_paired3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted_3.{sampling}.bam",
		sorted_bam_paired_ix3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted_3.{sampling}.bam.bai",
		sam_unpaired3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_3.{sampling}.sam",
		bam_unpaired3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_3.{sampling}.bam",
		sorted_bam_unpaired3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted_3.{sampling}.bam",
		sorted_bam_unpaired_ix3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted_3.{sampling}.bam.bai",
		scaffolds_pilon3=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_3_{sampling}/pilon.fasta"),
		scaffolds_pilon3_final=temp(dirs_dict["ASSEMBLY_DIR"] + "/pilon_3_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.unrenamed.fasta"),
		#round4
		sam_paired4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_4.{sampling}.sam",
		bam_paired4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_4.{sampling}.bam",
		sorted_bam_paired4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted_4.{sampling}.bam",
		sorted_bam_paired_ix4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted_4.{sampling}.bam.bai",
		sam_unpaired4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_4.{sampling}.sam",
		bam_unpaired4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_4.{sampling}.bam",
		sorted_bam_unpaired4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted_4.{sampling}.bam",
		sorted_bam_unpaired_ix4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted_4.{sampling}.bam.bai",
		scaffolds_pilon4=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_4_{sampling}/pilon.fasta"),
		scaffolds_pilon4_final=temp(dirs_dict["ASSEMBLY_DIR"] + "/pilon_4_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.unrenamed.fasta"),
		#final
		scaffolds=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}.unrenamed.fasta"),
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		#round1
		pilon_dir1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_1_{sampling}",
		db_name1=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_1_{sampling}",
		#round2
		pilon_dir2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_2_{sampling}",
		db_name2=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_2_{sampling}",
		#round3
		pilon_dir3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_3_{sampling}",
		db_name3=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_3_{sampling}",
		#round4
		pilon_dir4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_4_{sampling}",
		db_name4=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_4_{sampling}",
	message:
		"Correcting nanopore assembly using four rounds of Pilon"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/errorCorrectPilonPE/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	resources:
		mem_mb=16384
	shell:
		"""
		## ROUND 1
		# Index 1
		bowtie2-build -f {input.corrected2_racon} {params.db_name1} --threads {threads}
		# paired 1
		bowtie2 -X {params.max_insert_bp} -x {params.db_name1} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam_paired1} --threads {threads}
		samtools view -b -S {output.sam_paired1} > {output.bam_paired1}
		samtools sort {output.bam_paired1} -o {output.sorted_bam_paired1} -@ {threads}
		samtools index {output.sorted_bam_paired1}
		# unpaired 1
		bowtie2 -x {params.db_name1} -U {input.unpaired} -S {output.sam_unpaired1} --threads {threads}
		samtools view -b -S {output.sam_unpaired1} > {output.bam_unpaired1}
		samtools sort {output.bam_unpaired1} -o {output.sorted_bam_unpaired1}
		samtools index {output.sorted_bam_unpaired1}
		# PILON 1
		pilon -Xmx{resources.mem_mb}m --genome {input.corrected2_racon} --frags {output.sorted_bam_paired1} \
		--unpaired {output.sorted_bam_unpaired1} --outdir {params.pilon_dir1}

		## ROUND 2
		# Index 2
		bowtie2-build -f {output.scaffolds_pilon1} {params.db_name2} --threads {threads}
		# paired 2
		bowtie2 -X {params.max_insert_bp} -x {params.db_name2} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam_paired2} --threads {threads}
		samtools view -b -S {output.sam_paired2} > {output.bam_paired2}
		samtools sort {output.bam_paired2} -o {output.sorted_bam_paired2} -@ {threads}
		samtools index {output.sorted_bam_paired2}
		# unpaired 2
		bowtie2 -x {params.db_name2} -U {input.unpaired} -S {output.sam_unpaired2} --threads {threads}
		samtools view -b -S {output.sam_unpaired2} > {output.bam_unpaired2}
		samtools sort {output.bam_unpaired2} -o {output.sorted_bam_unpaired2}
		samtools index {output.sorted_bam_unpaired2}
		# PILON 2
		pilon -Xmx{resources.mem_mb}m --genome {output.scaffolds_pilon1} --frags {output.sorted_bam_paired2} \
		--unpaired {output.sorted_bam_unpaired2} --outdir {params.pilon_dir2}

		## ROUND 3
		# Index 3
		bowtie2-build -f {output.scaffolds_pilon2} {params.db_name3} --threads {threads}
		# paired 3
		bowtie2 -X {params.max_insert_bp} -x {params.db_name3} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam_paired3} --threads {threads}
		samtools view -b -S {output.sam_paired3} > {output.bam_paired3}
		samtools sort {output.bam_paired3} -o {output.sorted_bam_paired3} -@ {threads}
		samtools index {output.sorted_bam_paired3}
		# unpaired 3
		bowtie2 -x {params.db_name3} -U {input.unpaired} -S {output.sam_unpaired3} --threads {threads}
		samtools view -b -S {output.sam_unpaired3} > {output.bam_unpaired3}
		samtools sort {output.bam_unpaired3} -o {output.sorted_bam_unpaired3}
		samtools index {output.sorted_bam_unpaired3}
		# PILON 3
		pilon -Xmx{resources.mem_mb}m --genome {output.scaffolds_pilon2} --frags {output.sorted_bam_paired3} \
		--unpaired {output.sorted_bam_unpaired3} --outdir {params.pilon_dir3}

		## ROUND 4
		# Index 4
		bowtie2-build -f {output.scaffolds_pilon3} {params.db_name4} --threads {threads}
		# paired 4
		bowtie2 -X {params.max_insert_bp} -x {params.db_name4} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam_paired4} --threads {threads}
		samtools view -b -S {output.sam_paired4} > {output.bam_paired4}
		samtools sort {output.bam_paired4} -o {output.sorted_bam_paired4} -@ {threads}
		samtools index {output.sorted_bam_paired4}
		# unpaired 4
		bowtie2 -x {params.db_name4} -U {input.unpaired} -S {output.sam_unpaired4} --threads {threads}
		samtools view -b -S {output.sam_unpaired4} > {output.bam_unpaired4}
		samtools sort {output.bam_unpaired4} -o {output.sorted_bam_unpaired4}
		samtools index {output.sorted_bam_unpaired4}
		# PILON 4
		pilon -Xmx{resources.mem_mb}m --genome {output.scaffolds_pilon3} --frags {output.sorted_bam_paired4} \
		--unpaired {output.sorted_bam_unpaired4} --outdir {params.pilon_dir4}

		# Format scaffold names
		sed "s/>/>{wildcards.sample}_/g" {output.scaffolds_pilon4} > {output.scaffolds}
		cp {output.scaffolds_pilon1} {output.scaffolds_pilon1_final}
		cp {output.scaffolds_pilon2} {output.scaffolds_pilon2_final}
		cp {output.scaffolds_pilon3} {output.scaffolds_pilon3_final}
		cp {output.scaffolds_pilon4} {output.scaffolds_pilon4_final}
		"""

rule mergeAssembliesHYBRID:
	input:
		corrected_scaffolds=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample_nanopore}_"+ LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{{sampling}}.fasta", sample_nanopore=NANOPORE_SAMPLES),
		hybrid_contigs=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES),
	output:
		corrected_scaffolds=temp(dirs_dict["ASSEMBLY_DIR"] + "/merged_"+ LONG_ASSEMBLER + "_corrected_scaffolds.{{sampling}}.fasta"),
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
		merged_assembly_len=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_lengths.{sampling}.txt",
	message:
		"Merging assembled contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mergeAssembliesHYBRID/sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		cat {input.corrected_scaffolds} > {output.corrected_scaffolds}
		sed -i "s/=/_/g" {output.corrected_scaffolds}
		cat {input.hybrid_contigs} {output.corrected_scaffolds} > {output.merged_assembly}
		cat {output.merged_assembly} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
		$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.merged_assembly_len}
		"""
