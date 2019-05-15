ruleorder: mapReadsToContigsPE > mapReadsToContigsSE
ruleorder: getAbundancesPE > getAbundancesSE

rule createContigBowtieDb:
	input:
		contigs=dirs_dict["VIRAL_DIR"]+ "/{confidence}_confidence.{sampling}.fasta",
	output:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence.{sampling}.1.bt2",
		contigs_info=dirs_dict["VIRAL_DIR"]+ "/{confidence}_confidence.{sampling}.fasta.fai",
		contigs_lenght=dirs_dict["VIRAL_DIR"]+ "/{confidence}_confidence_lenghts.{sampling}.txt",
	params:
		prefix=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence.{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bowtie2-build -f {input.contigs} {params.prefix}
		#Get genome file
		samtools faidx {input.contigs}
		awk -F' ' '{{print $1"	"$2}}' {output.contigs_info} > {output.contigs_lenght}
		"""

rule mapReadsToContigsPE:
	input:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence.{sampling}.1.bt2",
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		sam=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence.{sampling}.sam",
	params:
		contigs=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence.{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		bowtie2 --non-deterministic -x {params.contigs} -1 {input.forward_paired} \
		-2 {input.reverse_paired} -U {input.unpaired} -S {output.sam} -p {threads}
		"""
rule mapReadsToContigsSE:
	input:
		contigs_bt2=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence.{sampling}.1.bt2",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq",
	output:
		sam=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence.{sampling}.sam",
		bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence.{sampling}.bam",
	params:
		contigs=dirs_dict["MAPPING_DIR"]+ "/{confidence}_confidence.{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		#Mapping reads to contigs
		bowtie2 --non-deterministic -x {params.contigs} -U {input.unpaired} -S {output.sam} -p {threads}
		#Sam to Bam
		samtools view -b -S {output.sam} > {output.bam}
		"""
		
rule filterBAM:
	input:
		bam=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence.{sampling}.bam",
	output:
		bam_sorted=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence_sorted.{sampling}.bam",
		bam_filtered=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence_sorted.{sampling}_filtered.bam",
		tpmean=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence_tpmean.{sampling}.tsv",
	params:
		out_dir=dirs_dict["MAPPING_DIR"]
	message:
		"Filtering reads in Bam file with BamM"
	conda:
		dirs_dict["ENVS_DIR"] + "/env2.yaml"
	threads: 1
	shell:
		"""
		samtools sort {input.bam} -o {output.bam_sorted}
		bamm filter --bamfile {output.bam_sorted} --percentage_id 0.95 --percentage_aln 0.9 -o {params.out_dir}
		bamm parse -c {output.tpmean} -m tpmean -b {output.bam_filtered}
		"""
rule filterContigs:
	input:
		bam_filtered=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence_sorted.{sampling}_filtered.bam",
		contigs=dirs_dict["VIRAL_DIR"]+ "/{confidence}_confidence.{sampling}.fasta",
		contigs_lenght=dirs_dict["VIRAL_DIR"]+ "/{confidence}_confidence_lenghts.{sampling}.txt",
	output:
		bam_cov=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence_filtered_genomecov.{sampling}.txt",
		cov_final=dirs_dict["MAPPING_DIR"]+ "/{sample}_{confidence}_confidence_filtered_coverage.{sampling}.txt",
	message:
		"Calculating breadth coverage contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 1
	shell:
		"""
		bedtools genomecov -d -ibam {input.bam_filtered} > {output.bam_cov}
		grep -v "0$" {output.bam_cov} | cut -f 1 | sort| uniq -c | sort -nr > {output.cov_final}
		sed -ie 's/^[[:space:]]*//' {output.cov_final}
		"""
