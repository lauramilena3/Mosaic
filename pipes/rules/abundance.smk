ruleorder: mapReadsToContigsPE > mapReadsToContigsSE
ruleorder: getAbundancesPE > getAbundancesSE

rule createContigBowtieDb:
	input:
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.fasta",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.fasta"
	output:
		high_contigs=dirs_dict["MAPPING_DIR"]+ "/high_confidence.{sampling}.1.bt2",
		low_contigs=dirs_dict["MAPPING_DIR"]+ "/low_confidence.{sampling}.1.bt2",
		high_contigs_info=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.fasta.fai",
		low_contigs_info=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.fasta.fai",
		high_contigs_lenght=dirs_dict["VIRAL_DIR"]+ "/high_confidence_lenghts.{sampling}.txt",
		low_contigs_lenght=dirs_dict["VIRAL_DIR"]+ "/low_confidence_lenghts.{sampling}.txt",
	params:
		high_contigs=dirs_dict["MAPPING_DIR"]+ "/high_confidence.{sampling}",
		low_contigs=dirs_dict["MAPPING_DIR"]+ "/low_confidence.{sampling}"
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bowtie2-build -f {input.high_contigs} {params.high_contigs}
		bowtie2-build -f {input.low_contigs} {params.low_contigs}
		#Get genome file
		samtools faidx {input.high_contigs}
		samtools faidx {input.low_contigs}
		awk -F' ' '{{print $1"	"$2}}' {output.high_contigs_info} > {output.high_contigs_lenght}
		awk -F' ' '{{print $1"	"$2}}' {output.low_contigs_info} > {output.low_contigs_lenght}
		"""

rule mapReadsToContigsPE:
	input:
		high_bt2=dirs_dict["MAPPING_DIR"]+ "/high_confidence.{sampling}.1.bt2",
		low_bt2=dirs_dict["MAPPING_DIR"]+ "/low_confidence.{sampling}.1.bt2",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		high_sam=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence.{sampling}.sam",
		low_sam=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence.{sampling}.sam",
	params:
		high_contigs=dirs_dict["MAPPING_DIR"]+ "/high_confidence.{sampling}",
		low_contigs=dirs_dict["MAPPING_DIR"]+ "/low_confidence.{sampling}"
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		bowtie2 --non-deterministic -x {params.high_contigs} -1 {input.forward_paired} -2 {input.reverse_paired} \
		-U {input.unpaired} -S {output.high_sam} -p {threads}
		bowtie2 --non-deterministic -x {params.low_contigs} -1 {input.forward_paired} -2 {input.reverse_paired} \
		-U {input.unpaired} -S {output.low_sam} -p {threads}
		"""
rule mapReadsToContigsSE:
	input:
		high_bt2=dirs_dict["MAPPING_DIR"]+ "/high_confidence.{sampling}.1.bt2",
		low_bt2=dirs_dict["MAPPING_DIR"]+ "/low_confidence.{sampling}.1.bt2",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		high_sam=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence.{sampling}.sam",
		low_sam=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence.{sampling}.sam",
		high_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence.{sampling}.bam",
		low_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence.{sampling}.bam",
	params:
		high_contigs=dirs_dict["MAPPING_DIR"]+ "/high_confidence.{sampling}",
		low_contigs=dirs_dict["MAPPING_DIR"]+ "/low_confidence.{sampling}"
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		#Mapping reads to contigs
		bowtie2 --non-deterministic -x {params.high_contigs} -U {input.unpaired} -S {output.high_sam} -p {threads}
		bowtie2 --non-deterministic -x {params.low_contigs} -U {input.unpaired} -S {output.low_sam} -p {threads}
		#Sam to Bam
		samtools view -b -S {output.high_sam} > {output.high_bam}
		samtools view -b -S {output.low_sam} > {output.low_bam}
		"""
		
rule filterBAM:
	input:
		high_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence.{sampling}.bam",
		low_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence.{sampling}.bam",
	output:
		high_bam_sorted=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence_sorted.{sampling}.bam",
		low_bam_sorted=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence_sorted.{sampling}.bam",
		high_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence_sorted.{sampling}_filtered.bam",
		low_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence_sorted.{sampling}_filtered.bam",
		high_bam_cov=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence_sorted.{sampling}.tsv",
		low_bam_cov=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence_sorted.{sampling}.tsv",
	params:
		out_dir=dirs_dict["MAPPING_DIR"]
	message:
		"Filtering reads in Bam file with BamM"
	conda:
		dirs_dict["ENVS_DIR"] + "/env2.yaml"
	threads: 1
	shell:
		"""
		samtools sort {input.high_bam} -o {output.high_bam_sorted}
		samtools sort {input.low_bam} -o {output.low_bam_sorted}
		bamm filter --bamfile {output.high_bam_sorted} --percentage_id 0.95 --percentage_aln 0.9 -o {params.out_dir}
		bamm filter --bamfile {output.low_bam_sorted} --percentage_id 0.95 --percentage_aln 0.9 -o {params.out_dir}
		bamm parse -c {output.high_bam_cov} -m tpmean -b {output.high_bam}
		bamm parse -c {output.low_bam_cov} -m tpmean -b {output.low_bam}
		"""
rule filterContigs:
	input:
		high_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence_sorted.{sampling}_filtered.bam",
		low_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence_sorted.{sampling}_filtered.bam",
		high_contigs=dirs_dict["VIRAL_DIR"]+ "/high_confidence.{sampling}.fasta",
		low_contigs=dirs_dict["VIRAL_DIR"]+ "/low_confidence.{sampling}.fasta",
		high_contigs_lenght=dirs_dict["VIRAL_DIR"]+ "/high_confidence_lenghts.{sampling}.txt",
		low_contigs_lenght=dirs_dict["VIRAL_DIR"]+ "/low_confidence_lenghts.{sampling}.txt",
	output:
		high_bam_cov=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence_filtered_genomecov.{sampling}.txt",
		low_bam_cov=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence_filtered_genomecov.{sampling}.txt",
		high_bam_final=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence_filtered_coverage.{sampling}.txt",
		low_bam_final=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence_filtered_coverage.{sampling}.txt",
	message:
		"Calculating breadth coverage contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bedtools genomecov -dz -ibam {input.high_bam} > {output.high_bam_cov}
		bedtools genomecov -dz -ibam {input.low_bam} > {output.low_bam_cov}
		cut -f 1 {output.high_bam_cov} | sort| uniq -c | sort -nr > {output.high_bam_final}
		cut -f 1 {output.low_bam_cov} | sort| uniq -c | sort -nr > {output.low_bam_final}
		"""

rule getAbundancesPE:
	input:
		high_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence_filtered_coverage.{sampling}.bam",
		low_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence_filtered_coverage.{sampling}.bam",
		high_contigs=dirs_dict["MAPPING_DIR"]+ "/high_confidence.{sampling}.fasta",
		low_contigs=dirs_dict["MAPPING_DIR"]+ "/low_confidence.{sampling}.fasta"
	output:
		high_bam_sorted=dirs_dict["MAPPING_DIR"]+ "/{sample}_final_results.{sampling}.csv",
	message:
		"Getting vOTU tables"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bamm 
		"""
rule getAbundancesSE:
	input:
		high_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_high_confidence_filtered_coverage.{sampling}.bam",
		low_bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_low_confidence_filtered_coverage.{sampling}.bam",
		high_contigs=dirs_dict["MAPPING_DIR"]+ "/high_confidence.{sampling}.fasta",
		low_contigs=dirs_dict["MAPPING_DIR"]+ "/low_confidence.{sampling}.fasta"
	output:
		high_bam_sorted=dirs_dict["MAPPING_DIR"]+ "/{sample}_final_results.{sampling}.csv",
	message:
		"Getting vOTU tables"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bamm 
		"""
