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
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
		unpaired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.txt"
	output:
		sam=dirs_dict["MAPPING_DIR"]+ "/bowtie_{sample}.{sampling}.sam",
		bam=dirs_dict["MAPPING_DIR"]+ "/bowtie_{sample}.{sampling}.bam",
	params:
		contigs=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		#bowtie2 --non-deterministic -x {params.contigs} -1 {input.forward_paired} \
		#-2 {input.reverse_paired} -U {input.unpaired} -S {output.sam} -p {threads}
		bbmap.sh ref={input.representatives} nodisk in1={input.forward_paired} in2={input.reverse_paired}  \
		outm={output.sam} threads={threads}
		#Sam to Bam
		samtools view -b -S {output.sam} > {output.bam}
		"""
rule mapReadsToContigsSE:
	input:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.1.bt2",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		sam=dirs_dict["MAPPING_DIR"]+ "/bowtie_{sample}.{sampling}.sam",
		bam=dirs_dict["MAPPING_DIR"]+ "/bowtie_{sample}.{sampling}.bam",
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
	message:
		"Calculating tpmean depth coverage"
	conda:
		dirs_dict["ENVS_DIR"] + "/env2.yaml"
	threads: 1
	shell:
		"""
		bamm parse -c {output.tpmean} -m tpmean -b {input.bam_sorted}
		"""
rule getBreadthCoverage:
	input:
		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/BamM_{sample}_sorted.{sampling}.bam",
	output:
		bam_cov=dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_genomecov.{sampling}.txt",
		cov_final=dirs_dict["MAPPING_DIR"]+ "/bedtools_{sample}_coverage.{sampling}.txt",
	message:
		"Calculating breadth coverage contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bedtools genomecov -dz -ibam {input.bam_sorted} > {output.bam_cov}
		cut -f 1 {output.bam_cov} | sort| uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > {output.cov_final}
		"""
