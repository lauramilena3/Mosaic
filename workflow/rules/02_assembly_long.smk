#ruleorder: asemblyCanuPOOLED > asemblyCanu
ruleorder: hybridAsemblySpades > shortReadAsemblySpadesPE
#ruleorder: errorCorrectPE > errorCorrectSE
ruleorder: assemblyStatsHYBRID > assemblyStatsILLUMINA
ruleorder: mergeAssembliesHYBRID > mergeAssembliesSHORT


if POOLED==True:
	ruleorder: symlinkPooled>subsampleReadsNanopore

	rule symlinkPooled:
		input:
			pooled=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.{{sampling}}.fastq", sample_nanopore=NANOPORE_SAMPLES),
		output:
			expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{{sampling}}.fastq", sample=SAMPLES),
		message:
			"Creating symbolic links from pooled sample"
		threads: 11
		shell:
			"""
			for destination in {output}
			do
				ln -s {input.pooled} "$destination"
			done
			#add a merged illumina
			"""
rule hybridAsemblySpades:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq",
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq"
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/filtered_list.txt")
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}/scaffolds.fasta",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_{sampling}")
	message:
		"Assembling hybrid reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 16
	shell:
		"""
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		--meta -t {threads} --only-assembler --nanopore {input.nanopore} --memory 350
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		"""

# rule asemblyCanuPOOLED:
# 	input:
# 		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq",
# 		canu_dir=config['canu_dir']
# 	output:
# 		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_canu_{sampling}/" +config['nanopore_pooled_name'] + ".contigs.fasta",
# 		scaffolds_all=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_canu.{{sampling}}.fasta", sample=SAMPLES)
# 	message:
# 		"Assembling Nanopore reads with Canu"
# 	params:
# 		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/"+ config['nanopore_pooled_name']+ "_canu_{sampling}",
# 		assembly=dirs_dict["ASSEMBLY_DIR"],
# 		sample_list=" ".join(SAMPLES),
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 4
# 	shell:
# 		"""
# 		./{config[canu_dir]}/canu genomeSize=45m minReadLength=1000 -p \
# 		contigFilter="{config[min_cov]} {config[min_len]} 1.0 1.0 2" \
# 		corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 \
# 		redMemory=32 oeaMemory=32 batMemory=200 -nanopore-raw {input.nanopore} \
# 		-d {params.assembly_dir} -p {config[nanopore_pooled_name]} useGrid=false executiveThreads={threads}
# 		for sample in {params.sample_list}
# 		do
# 			cat {output.scaffolds} | sed s"/ /_/"g  > {params.assembly}/${{sample}}_contigs_canu.{wildcards.sampling}.fasta
# 		done
# 		"""

rule asemblyCanu:
	input:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq",
		canu_dir=config['canu_dir'],
	output:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/canu_{sample}_{sampling}/{sample}.contigs.fasta",
		scaffolds_final=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_canu.{sampling}.fasta"
	message:
		"Assembling Nanopore reads with Canu"
	params:
		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/canu_{sample}_{sampling}"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		./{config[canu_dir]}/canu genomeSize=100k minReadLength=1000 -p \
		contigFilter="{config[min_cov]} {config[min_len]} 1.0 1.0 2" \
		corOutCoverage=all corMhapSensitivity=high correctedErrorRate=0.105 corMinCoverage=0 \
		corMaxEvidenceCoverageLocal=10 corMaxEvidenceCoverageGlobal=10 \
		redMemory=32 oeaMemory=32 batMemory=200 -nanopore {input.nanopore} \
		-d {params.assembly_dir} -p {wildcards.sample} useGrid=false maxThreads={threads}
		cp {output.scaffolds} {output.scaffolds_final}
		sed -i s"/ /_/"g {output.scaffolds_final}
		"""

rule asemblyFlye:
	input:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq",
	output:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample}_{sampling}/assembly.fasta",
		scaffolds_final=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_flye.{sampling}.fasta"
	message:
		"Assembling Nanopore reads with Flye"
	params:
		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample}_{sampling}",
		genome_size="20m"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		flye --nano-raw {input.nanopore} --out-dir {params.assembly_dir} --genome-size {params.genome_size} --meta --threads {threads}
		cp {output.scaffolds} {output.scaffolds_final}
		"""
rule errorCorrectRacon_1st:
	input:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
	output:
		overlap=dirs_dict["ASSEMBLY_DIR"] + "/minimap2_1_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.paf",
		corrected=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
	message:
		"Correcting nanopore assembly with Racon, long reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		minimap2 -t {threads} {input.scaffolds} {input.nanopore} > {output.overlap}
		racon -t {threads}  {input.nanopore} {output.overlap} {input.scaffolds} > {output.corrected}
		"""
rule errorCorrectRacon_2nd:
	input:
		corrected1=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		illumina=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_racon_illumina_clean.{sampling}.fastq"),
		overlap=dirs_dict["ASSEMBLY_DIR"] + "/minimap2_2_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.sam",
		corrected=dirs_dict["ASSEMBLY_DIR"] + "/racon_polished_{sample}_contigs_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
	params:
		racon_merge="./scripts/merge_illumina_racon.py",
	message:
		"Correcting nanopore assembly with Racon, short reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		{params.racon_merge} {input.forward_paired} {input.reverse_paired} {input.unpaired} > {output.illumina}
		minimap2 -t {threads} -ax sr {input.corrected1} {output.illumina} > {output.overlap}
		racon -t {threads} {output.illumina} {output.overlap} {input.corrected1} > {output.corrected}
		"""

rule errorCorrectPilonPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_1_"+ LONG_ASSEMBLER + ".{sampling}.fasta",
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds.{sampling}.fasta"),
		sam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired.{sampling}.sam",
		bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired.{sampling}.bam",
		sorted_bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{sampling}.bam",
		sorted_bam_paired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{sampling}.bam.bai",
		sam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{sampling}.sam",
		bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{sampling}.bam",
		sorted_bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{sampling}.bam",
		sorted_bam_unpaired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{sampling}.bam.bai"
	params:
		pilon_dir=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_{sampling}",
		scaffolds_pilon=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_{sampling}/pilon.fasta"),
		db_name=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_{sampling}"
	message:
		"Correcting nanopore assembly with Pilon"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	resources:
		mem_mb=16384
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
		#PILON
		pilon -Xmx{resources.mem_mb}m --genome {input.scaffolds} --frags {output.sorted_bam_paired} \
		--unpaired {output.sorted_bam_unpaired} --outdir {params.pilon_dir}
		cp {params.scaffolds_pilon} {output.scaffolds}
		"""


# rule asemblyFlye2nd:
# 	input:
# 		corrected_scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds.{sampling}.fasta"),
# 		hybrid_contigs=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta"),
# 	output:
# 		scaffolds_flye2=dirs_dict["ASSEMBLY_DIR"] + "/flye_combined_assembly_{sample}.{sampling}.fasta",
# 		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/flye_combined_assembly_{sample}_{sampling}/assembly.fasta",
# 	message:
# 		"Assembling metaSPAdes and canu subassemblies with Flye"
# 	params:
# 		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/flye_combined_assembly_{sample}_{sampling}",
# 		genome_size="20m"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 4
# 	shell:
# 		"""
# 		flye --subassemblies {input.corrected_scaffolds} {input.hybrid_contigs} --out-dir {params.assembly_dir} --genome-size {params.genome_size} --meta --threads {threads}
# 		cp {output.scaffolds} {output.scaffolds_flye2}
# 		sed -i "s/>/>{wildcards.sample}_/g" {output.scaffolds_flye2}
# 		"""

# rule errorCorrectSE:
# 	input:
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
# 		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_canu.{sampling}.fasta"
# 	output:
# 		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds.{sampling}.fasta"),
# 		sam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{sampling}.sam",
# 		bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired.{sampling}.bam",
# 		sorted_bam_unpaired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{sampling}.bam",
# 		sorted_bam_unpaired_ix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_unpaired_sorted.{sampling}.bam.bai",
# 	params:
# 		pilon_dir=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_{sampling}",
# 		scaffolds_pilon=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_pilon_{sampling}/pilon.fasta"),
# 		db_name=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_bowtieDB_{sampling}",
# 	message:
# 		"Correcting nanopore assembly with Pilon"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		bowtie2-build -f {input.scaffolds} {params.db_name}
# 		#unpaired
# 		bowtie2 -x {params.db_name} -U {input.unpaired} -S {output.sam_unpaired}
# 		samtools view -b -S {output.sam_unpaired} > {output.bam_unpaired}
# 		samtools sort {output.bam_unpaired} -o {output.sorted_bam_unpaired}
# 		samtools index {output.sorted_bam_unpaired}
# 		#PILON
# 		pilon --genome {input.scaffolds} --unpaired {output.sorted_bam_unpaired} --outdir {params.pilon_dir}
# 		cp {params.scaffolds_pilon} {output.scaffolds}
# 		"""

rule assemblyStatsHYBRID:
	input:
		quast_dir=(config["quast_dir"]),
		scaffolds_spades=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES),
		scaffolds_long=expand(dirs_dict["ASSEMBLY_DIR"] + "/flye_combined_assembly_{sample}.{{sampling}}.fasta", sample=SAMPLES),
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast_{sampling}"),
		quast_txt=dirs_dict["ASSEMBLY_DIR"] + "/assembly_quast_report.{sampling}.txt",
	message:
		"Creating assembly stats with quast"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		{input.quast_dir}/quast.py {input.scaffolds_long} {input.scaffolds_spades} -o {output.quast_report_dir}
		cp {output.quast_report_dir}/report.txt {output.quast_txt}
		"""

rule scoreALE:
	input:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{contigs}.fasta",
		sorted_bam_paired=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{sampling}_to_{contigs}.bam",
		ALE_dir=directory(config['ALE_dir']),
	output:
		ale=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_paired_sorted.{sampling}_to_{contigs}.ale",
	message:
		"Creating assembly stats with ALE"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		{input.ALE_dir}/src/ALE {input.sorted_bam_paired} {input.scaffolds} {output.ale}
		"""

rule mergeAssembliesHYBRID:
	input:
		corrected_scaffolds=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER + "_corrected_scaffolds.{{sampling}}.fasta", sample=SAMPLES),
		hybrid_contigs=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta", sample=SAMPLES),
	output:
		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
		merged_assembly_len=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_lengths.{sampling}.txt",
	message:
		"Merging assembled contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		cat {input.hybrid_contigs} {input.corrected_scaffolds} > {output.merged_assembly}
		cat {output.merged_assembly} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
		$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.merged_assembly_len}
		"""

# rule mergeAssembliesPOOLED:
# 	input:
# 		scaffolds_spades=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{{sampling}}.fasta",sample=SAMPLES),
# 		scaffolds_long=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample_nanopore}_"+ LONG_ASSEMBLER + "_filtered_scaffolds.{{sampling}}.fasta", sample_nanopore=NANOPORE_SAMPLES),
# 	output:
# 		merged_assembly=(dirs_dict["VIRAL_DIR"] + "/merged_scaffolds.{sampling}.fasta"),
# 		merged_assembly_len=dirs_dict["VIRAL_DIR"] + "/merged_scaffolds_lengths.{sampling}.txt",
# 	message:
# 		"Merging assembled contigs"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	threads: 1
# 	shell:
# 		"""
# 		cat {input.scaffolds_long} {input.scaffolds_spades} > {output.merged_assembly}
# 		cat {output.merged_assembly} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} \
# 		$0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output.merged_assembly_len}
# 		"""
