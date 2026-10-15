# Isolate mapping always uses every QC-passed read, irrespective of sampling flags.
def mapping_read_path(wildcards, side, subset=False):
	if subset and not ISOLATES:
		return dirs_dict["ASSEMBLY_TEST"] + f"/2M_{wildcards.sample}_{side}_paired_clean.{getattr(wildcards, 'sampling', 'tot')}.fastq.gz"
	return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_{side}_paired_clean.{('tot' if ISOLATES else getattr(wildcards, 'sampling', 'tot'))}.fastq.gz"

def mapping_orphan_reads(wildcards):
	return [dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_unpaired_clean.tot.fastq.gz"] if ISOLATES else []

def mapping_orphan_option(wildcards, input):
	import shlex
	return " -U " + shlex.quote(input.unpaired[0]) if input.unpaired else ""

#ruleorder: mapReadsToContigsPE > mapReadsToContigsSE

# Track all six project index files in both builders and consumers. Snakemake
# removes temporary indexes only after every dependent mapping job finishes.

rule subsampleReadsIllumina_PE_mapping:
	input:
		paired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}_read_count.txt",),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		forward_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["ASSEMBLY_TEST"] + "/2M_{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	params:
		n_reads=2000000
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/subsampleReadsIllumina_PE_mapping/sample={sample}__sampling={sampling}.tsv"
	threads: 1
	# resources:
	# 	mem_mb=4000
	shell:
		"""
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads={params.n_reads} sampleseed=1
		"""

rule subsampleReadsIllumina_PE_mapping_7M:
	input:
		paired_sizes=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}_read_count.txt",),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		forward_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/7M_{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=temp(dirs_dict["ASSEMBLY_TEST"] + "/7M_{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	params:
		n_reads=7000000
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/subsampleReadsIllumina_PE_mapping_7M/sample={sample}__sampling={sampling}.tsv"
	threads: 1
	# resources:
	# 	mem_mb=4000
	shell:
		"""
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads={params.n_reads} sampleseed=1
		"""

rule subsampleReadsIllumina_PE_vOTU_mapping:
	input:
		viral_subsampling=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sub_sampling_reads.txt"),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
	output:
		forward_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.sub.fastq.gz"),
		reverse_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.sub.fastq.gz"),
		viral_subsampling=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sub_sampling_reads_final.txt"),
	message:
		"Subsampling Illumina reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/subsampleReadsIllumina_PE_vOTU_mapping/sample={sample}.tsv"
	threads: 1
	shell:
		"""
		reads=$(cat {input.viral_subsampling})
		reformat.sh in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} reads=$reads sampleseed=1
		echo $(( $(zgrep -Ec "$" {output.forward_paired}) / 4 )) > {output.viral_subsampling}
		"""

rule buildBowtieDB_assembly:
	input:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
	output:
		contigs_bt2_1=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.1.bt2"),
		contigs_bt2_2=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.2.bt2"),
		contigs_bt2_3=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.3.bt2"),
		contigs_bt2_4=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.4.bt2"),
		contigs_bt2_rev_1=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.rev.1.bt2"),
		contigs_bt2_rev_2=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.rev.2.bt2"),
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_assembly/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		if [ -s {input.scaffolds:q} ]; then
			bowtie2-build {input.scaffolds:q} {params.prefix:q} --threads {threads}
		else
			touch {output:q}
		fi
		"""

rule buildBowtieDB_genes:
	input:
		NR_fna_150=dirs_dict["ANNOTATION"]+ "/predicted_genes_NR_95_85_150bp_tot.fna",
	output:
		NR_bt2_150=temp(expand(dirs_dict["ANNOTATION"] + "/predicted_genes_NR_95_85_150bp_tot.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"])),
	params:
		prefix=dirs_dict["ANNOTATION"] + "/predicted_genes_NR_95_85_150bp_tot",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_genes/tot.tsv"
	threads: 8
	shell:
		"""
		bowtie2-build {input.NR_fna_150} {params.prefix} --threads {threads}
		"""

def input_mapping(wildcards):
	input_list = []
	if METAGENOME:
		input_list.append(dirs_dict["ASSEMBLY_TEST"] + f"/2M_{wildcards.sample}_forward_paired_clean.{wildcards.sampling}.fastq.gz")
		input_list.append(dirs_dict["ASSEMBLY_TEST"] + f"/2M_{wildcards.sample}_reverse_paired_clean.{wildcards.sampling}.fastq.gz")
	else:
		input_list.append(dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_forward_paired_clean.tot.fastq.gz")
		input_list.append(dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_reverse_paired_clean.tot.fastq.gz")
	return input_list
	
rule stat_mapReadsToAssembly:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
		contigs_bt2_1=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.1.bt2"),
		contigs_bt2_2=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.2.bt2"),
		contigs_bt2_3=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.3.bt2"),
		contigs_bt2_4=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.4.bt2"),
		contigs_bt2_rev_1=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.rev.1.bt2"),
		contigs_bt2_rev_2=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.rev.2.bt2"),
		forward_paired=lambda wildcards: input_mapping(wildcards)[0],
		reverse_paired=lambda wildcards: input_mapping(wildcards)[1],
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}_sorted.bam.bai"),
		filtered_bam=(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}_filtered.bam") if ISOLATES else temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_{sample}_assembled_contigs.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_assembled_contigs.{sampling}.txt",
		basecov=(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs.{sampling}_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		covstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_assembled_contigs.{sampling}_covstats.txt",
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unpaired=mapping_orphan_option,
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}",
		alignment_evidence="--include-secondary --exclude-supplementary" if ISOLATES else "",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/stat_mapReadsToAssembly/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		if [ -s {input.fasta:q} ]; then
			bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --all --very-sensitive
		else
			printf '@HD\tVN:1.6\tSO:coordinate\n' > {output.sam:q}
		fi
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		if [ -s {input.fasta:q} ]; then
			coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 {params.alignment_evidence} -t {threads}
		else
			cp {output.sorted_bam:q} {output.filtered_bam:q}
		fi
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		if [ -s {input.fasta:q} ]; then
			for basecov in {output.basecov:q}; do
				bedtools genomecov -dz -ibam {output.filtered_bam:q} | gzip -c > "$basecov"
			done
			coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm -o {output.covstats}
		else
			for basecov in {output.basecov:q}; do
				printf '' | gzip -c > "$basecov"
			done
			printf 'Contig\tMean\tLength\tCovered Bases\tRead Count\tVariance\tTrimmed Mean\tRPKM\n' > {output.covstats:q}
		fi
		"""

rule buildBowtieDB_viral:
	input:
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		contigs_bt2=temp(expand(dirs_dict["VIRAL_DIR"] + "/{{sample}}_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"])),
	params:
		prefix=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_viral/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2-build {input.positive_contigs} {params.prefix} --threads {threads}
		"""

rule stat_mapReadsToViral:
	input:
		contigs_bt2=expand(dirs_dict["VIRAL_DIR"] + "/{{sample}}_" + VIRAL_CONTIGS_BASE + ".{{sampling}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"]),
		forward_paired=lambda wc: mapping_read_path(wc, "forward", subset=True),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse", subset=True),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_{sample}_viral_contigs.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_viral_contigs.{sampling}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_viral_contigs.{sampling}_covstats.txt",
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unpaired=mapping_orphan_option,
		prefix=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/stat_mapReadsToViral/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --all --very-sensitive
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		"""

rule buildBowtieDB_derreplicated:
	input:
		derreplicated_positive_contigs=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}.fasta",
	output:
		contigs_bt2=temp(expand(dirs_dict["vOUT_DIR"] + "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{{sampling}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"])),
	params:
		prefix=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_derreplicated/sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2-build {input.derreplicated_positive_contigs} {params.prefix} --threads {threads}
		"""

rule stat_mapReadsToDerreplicated:
	input:
		contigs_bt2=expand(dirs_dict["vOUT_DIR"] + "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{{sampling}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"]),
		forward_paired=lambda wc: mapping_read_path(wc, "forward", subset=True),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse", subset=True),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_{sample}_derreplicated_contigs.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_derreplicated_contigs.{sampling}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_derreplicated_contigs.{sampling}_covstats.txt",
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unpaired=mapping_orphan_option,
		prefix=dirs_dict["vOUT_DIR"]+ "/combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq.{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/stat_mapReadsToDerreplicated/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --all --very-sensitive
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		"""

rule buildBowtieDB_unfiltered:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
	output:
		contigs_bt2=temp(expand(dirs_dict["MAPPING_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".{{sampling}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"])),
	params:
		prefix=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_unfiltered/sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2-build {input.representatives} {params.prefix} --threads {threads}
		"""

rule stat_mapReadsToUnfiltered:
	input:
		contigs_bt2=expand(dirs_dict["MAPPING_DIR"] + "/" + REPRESENTATIVE_CONTIGS_BASE + ".{{sampling}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"]),
		forward_paired=lambda wc: mapping_read_path(wc, "forward", subset=True),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse", subset=True),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}.bam"),
		sorted_bam=(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_{sample}_unfiltered_contigs.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_flagstats_filtered_{sample}_unfiltered_contigs.{sampling}.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/STATS_FILES/bowtie2_{sample}_unfiltered_contigs.{sampling}_covstats.txt",
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unpaired=mapping_orphan_option,
		prefix=dirs_dict["MAPPING_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/stat_mapReadsToUnfiltered/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --all --very-sensitive
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		"""

rule buildBowtieDB_filtered:
	input:
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
	output:
		contigs_bt2=temp(expand(dirs_dict["MAPPING_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{{sampling}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"])),
	params:
		prefix=dirs_dict["MAPPING_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_filtered/sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2-build {input.filtered_representatives} {params.prefix} --threads {threads}
		"""

rule mapReadsToContigsPE:
	input:
		contigs_bt2=expand(dirs_dict["MAPPING_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"]),
		forward_paired=lambda wc: mapping_read_path(wc, "forward"),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse"),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}.sam"),
		bam=(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_{sample}.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_filtered_{sample}.{sampling}.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/bowtie2_flagstats_unique_{sample}.{sampling}.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique_sorted.bam"),
		unique_sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique_sorted.bam.bai"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique_covstats.txt",
		basecov=(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		unique_basecov=(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_{sampling}_unique_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unique_flags="-F 2308" if ISOLATES else "-f 0x2",
		unpaired=mapping_orphan_option,
		prefix=dirs_dict["MAPPING_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapReadsToContigsPE/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --all --very-sensitive
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ {threads} -h {params.unique_flags} {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ {threads} -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ {threads} {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		for basecov in {output.basecov:q}; do
			bedtools genomecov -dz -ibam {output.filtered_bam} | gzip -c > "$basecov"
		done
		for basecov in {output.unique_basecov:q}; do
			bedtools genomecov -dz -ibam {output.unique_sorted_bam} | gzip -c > "$basecov"
		done
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""
rule call_SNPs_sub:
	input:
		filtered_representatives=dirs_dict["vOUT_DIR"]+ "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
		sorted_bam=(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_sub_sorted.bam"),
	output:
		snp_temp=temp(dirs_dict["MAPPING_DIR"]+ "/{sample}_sub_SNP_calls.bcf"),
		snp=(dirs_dict["MAPPING_DIR"]+ "/{sample}_sub_SNP_calls.tsv"),
	message:
		"Calling SNPs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/call_SNPs_sub/sample={sample}.tsv"
	threads: 8
	shell:
		"""
		bcftools mpileup -Ou -I -d 1000 -f {input.filtered_representatives} {input.sorted_bam} | bcftools call -mv -Ob --ploidy 1 -o {output.snp_temp}
		bcftools view -i '%QUAL>=20' {output.snp_temp} | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" -o {output.snp}
		"""

rule buildBowtieDB_contaminants:
	input:
		contaminants=dirs_dict["CONTAMINANTS_DIR_POST"]+ "/{contaminant}.fasta",
	output:
		contigs_bt2=temp(expand(dirs_dict["CONTAMINANTS_DIR_POST"] + "/{{contaminant}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"])),
	params:
		prefix=dirs_dict["CONTAMINANTS_DIR_POST"]+ "/{contaminant}",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_contaminants/contaminant={contaminant}.tsv"
	threads: 32
	shell:
		"""
		bowtie2-build {input.contaminants} {params.prefix} --threads {threads}
		"""


rule mapReads_contaminants:
	input:
		contigs_bt2=expand(dirs_dict["CONTAMINANTS_DIR_POST"] + "/{{contaminant}}.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"]),
		forward_paired=lambda wc: mapping_read_path(wc, "forward", subset=True),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse", subset=True),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_flagstats_{sample}_{contaminant}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_flagstats_filtered_{sample}.{contaminant}.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_flagstats_unique_{sample}.{contaminant}.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique_sorted.bam"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique_covstats.txt",
		basecov=(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		unique_basecov=(dirs_dict["MAPPING_DIR"]+ "/CONTAMINANTS/bowtie2_{sample}_{contaminant}_unique_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unique_flags="-F 2308" if ISOLATES else "-f 0x2",
		unpaired=mapping_orphan_option,
		prefix=dirs_dict["CONTAMINANTS_DIR_POST"]+ "/{contaminant}",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapReads_contaminants/contaminant={contaminant}__sample={sample}.tsv"
	threads: 16
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --fast
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ {threads} -h {params.unique_flags} {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ {threads} -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ {threads} {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		for basecov in {output.basecov:q}; do
			bedtools genomecov -dz -ibam {output.filtered_bam} | gzip -c > "$basecov"
		done
		for basecov in {output.unique_basecov:q}; do
			bedtools genomecov -dz -ibam {output.unique_sorted_bam} | gzip -c > "$basecov"
		done
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""

# RefSeq mapping is separate from both the vOTU catalogue and user references.
# Keep the index in the results directory, not beside the shared source database.
REFSEQ_MAPPING_DIR=dirs_dict["MAPPING_DIR"] + "/REFSEQ_VIRAL"
REFSEQ_INDEX_PREFIX=REFSEQ_MAPPING_DIR + "/INDEX/RefSeqViral"
REFSEQ_INDEX_FILES=expand(REFSEQ_INDEX_PREFIX + ".{part}.bt2l", part=["1", "2", "3", "4", "rev.1", "rev.2"])

rule buildBowtieDB_RefSeq:
	input:
		fasta=lambda wildcards: config["RefSeqViral_db"],
	output:
		index=temp(REFSEQ_INDEX_FILES),
	params:
		prefix=REFSEQ_INDEX_PREFIX,
	message:
		"Creating the RefSeq Viral Bowtie2 index"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_RefSeq/tot.tsv"
	threads: 16
	shell:
		"""
		bowtie2-build --large-index --threads {threads} {input.fasta:q} {params.prefix:q}
		"""

rule mapReads_RefSeq:
	input:
		index=REFSEQ_INDEX_FILES,
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot.sam"),
		bam=REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot.bam",
		sorted_bam=temp(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_sorted.bam"),
		sorted_bam_idx=temp(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_sorted.bam.bai"),
		filtered_bam=temp(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_filtered.bam"),
		flagstats=REFSEQ_MAPPING_DIR + "/bowtie2_flagstats_RefSeqViral_{sample}.tot.txt",
		flagstats_filtered=REFSEQ_MAPPING_DIR + "/bowtie2_flagstats_filtered_RefSeqViral_{sample}.tot.txt",
		flagstats_unique=REFSEQ_MAPPING_DIR + "/bowtie2_flagstats_unique_RefSeqViral_{sample}.tot.txt",
		unique_bam=temp(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique.bam"),
		unique_sorted_bam=temp(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_sorted.bam"),
		unique_sorted_bam_idx=temp(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_sorted.bam.bai"),
		covstats=REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_covstats.txt",
		covstats_unique=REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_covstats.txt",
		basecov=(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		unique_basecov=(REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot_unique_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unique_flags="-F 2308" if ISOLATES else "-f 0x2",
		unpaired=mapping_orphan_option,
		prefix=REFSEQ_INDEX_PREFIX,
	message:
		"Mapping cleaned paired-end reads to RefSeq Viral"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapReads_RefSeq/sample={sample}.tsv"
	log:
		bowtie2=REFSEQ_MAPPING_DIR + "/bowtie2_RefSeqViral_{sample}_tot.log",
	threads: 16
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix:q} -1 {input.forward_paired:q} -2 {input.reverse_paired:q} {params.unpaired} -S {output.sam:q} --threads {threads} --no-unal --fast 2> {log.bowtie2:q}
		samtools view -@ {threads} -bS {output.sam:q} -o {output.bam:q}
		samtools sort -@ {threads} {output.bam:q} -o {output.sorted_bam:q}
		samtools index {output.sorted_bam:q}
		samtools flagstat {output.sorted_bam:q} > {output.flagstats:q}
		coverm filter -b {output.sorted_bam:q} -o {output.filtered_bam:q} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam:q} > {output.flagstats_filtered:q}
		# No-XS primary mappings include orphans/discordant reads in isolate mode; other modes retain the proper-pair heuristic.
		samtools view -@ {threads} -h {params.unique_flags} {output.filtered_bam:q} | awk '/^@/ || !/XS:i:/' | samtools view -@ {threads} -b -o {output.unique_bam:q} -
		samtools sort -@ {threads} {output.unique_bam:q} -o {output.unique_sorted_bam:q}
		samtools index {output.unique_sorted_bam:q}
		samtools flagstat {output.unique_bam:q} > {output.flagstats_unique:q}
		for basecov in {output.basecov:q}; do
			bedtools genomecov -dz -ibam {output.filtered_bam:q} | gzip -c > "$basecov"
		done
		for basecov in {output.unique_basecov:q}; do
			bedtools genomecov -dz -ibam {output.unique_sorted_bam:q} | gzip -c > "$basecov"
		done
		coverm contig -b {output.filtered_bam:q} -m mean length covered_bases count variance trimmed_mean rpkm -t {threads} -o {output.covstats:q}
		coverm contig -b {output.unique_sorted_bam:q} -m mean length covered_bases count variance trimmed_mean rpkm -t {threads} -o {output.covstats_unique:q}
		"""

ALL_ASSEMBLED_INDEX_PREFIX=ALL_ASSEMBLED_MAPPING_DIR + "/INDEX/{assembled_reference}"
ALL_ASSEMBLED_INDEX_FILES=expand(ALL_ASSEMBLED_MAPPING_DIR + "/INDEX/{{assembled_reference}}.{part}.bt2l", part=["1", "2", "3", "4", "rev.1", "rev.2"])

def all_assembled_reference_fasta(wildcards):
	if wildcards.assembled_reference == "all_assembled_top":
		return ALL_ASSEMBLED_TOP_PREFIX + ".fasta"
	return ALL_ASSEMBLED_DIR + "/all_assembled_contigs_derreplicated_rep_seq.tot.fasta"

rule buildBowtieDB_all_assembled:
	input:
		fasta=all_assembled_reference_fasta,
	output:
		index=temp(ALL_ASSEMBLED_INDEX_FILES),
	params:
		prefix=ALL_ASSEMBLED_INDEX_PREFIX,
	wildcard_constraints:
		assembled_reference="all_assembled|all_assembled_top",
	message:
		"Creating the {wildcards.assembled_reference} Bowtie2 index"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_all_assembled/reference={assembled_reference}.tsv"
	threads: 16
	shell:
		"""
		if [ -s {input.fasta:q} ]; then
			bowtie2-build --large-index --threads {threads} {input.fasta:q} {params.prefix:q}
		else
			# Empty-catalogue markers, not valid indices; the mapping rule skips Bowtie2.
			touch {output.index:q}
		fi
		"""

rule mapReads_all_assembled:
	input:
		index=ALL_ASSEMBLED_INDEX_FILES,
		fasta=all_assembled_reference_fasta,
		forward_paired=lambda wc: mapping_read_path(wc, "forward", subset=True),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse", subset=True),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot.sam"),
		bam=ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot.bam",
		sorted_bam=temp(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_sorted.bam"),
		sorted_bam_idx=temp(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_sorted.bam.bai"),
		filtered_bam=temp(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_filtered.bam"),
		flagstats=ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_flagstats_{assembled_reference}_{sample}.tot.txt",
		flagstats_filtered=ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_flagstats_filtered_{assembled_reference}_{sample}.tot.txt",
		mapped_pairs=ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_mapped_pairs_filtered_{assembled_reference}_{sample}.tot.txt",
		flagstats_unique=ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_flagstats_unique_{assembled_reference}_{sample}.tot.txt",
		unique_bam=temp(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_unique.bam"),
		unique_sorted_bam=temp(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_unique_sorted.bam"),
		unique_sorted_bam_idx=temp(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_unique_sorted.bam.bai"),
		covstats=ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_covstats.txt",
		covstats_unique=ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_unique_covstats.txt",
		basecov=(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		unique_basecov=(ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot_unique_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unique_flags="-F 2308" if ISOLATES else "-f 0x2",
		unpaired=mapping_orphan_option,
		prefix=ALL_ASSEMBLED_INDEX_PREFIX,
	wildcard_constraints:
		assembled_reference="all_assembled|all_assembled_top",
		sample="|".join(re.escape(sample) for sample in SAMPLES) or "(?!)",
	message:
		"Mapping cleaned reads to {wildcards.assembled_reference} (all paired and orphan reads in isolate mode)"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapReads_all_assembled/reference={assembled_reference}__sample={sample}.tsv"
	log:
		bowtie2=ALL_ASSEMBLED_MAPPING_DIR + "/bowtie2_{assembled_reference}_{sample}_tot.log",
	threads: 16
	shell:
		"""
		if [ -s {input.fasta:q} ]; then
			bowtie2 -X {params.max_insert_bp} -x {params.prefix:q} -1 {input.forward_paired:q} -2 {input.reverse_paired:q} {params.unpaired} \
				-S {output.sam:q} --threads {threads} --no-unal --very-sensitive 2> {log.bowtie2:q}
		else
			printf '@HD\tVN:1.6\tSO:coordinate\n' > {output.sam:q}
			printf 'Empty full catalogue; mapping skipped.\n' > {log.bowtie2:q}
		fi
		samtools view -@ {threads} -bS {output.sam:q} -o {output.bam:q}
		samtools sort -@ {threads} {output.bam:q} -o {output.sorted_bam:q}
		samtools index {output.sorted_bam:q}
		samtools flagstat {output.sorted_bam:q} > {output.flagstats:q}
		if [ -s {input.fasta:q} ]; then
			coverm filter -b {output.sorted_bam:q} -o {output.filtered_bam:q} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		else
			cp {output.sorted_bam:q} {output.filtered_bam:q}
		fi
		samtools flagstat {output.filtered_bam:q} > {output.flagstats_filtered:q}
		samtools view -@ {threads} -c -f 0x42 -F 0x104 {output.filtered_bam:q} > {output.mapped_pairs:q}
		samtools view -@ {threads} -h {params.unique_flags} {output.filtered_bam:q} | awk '/^@/ || !/XS:i:/' | samtools view -@ {threads} -b -o {output.unique_bam:q} -
		samtools sort -@ {threads} {output.unique_bam:q} -o {output.unique_sorted_bam:q}
		samtools index {output.unique_sorted_bam:q}
		samtools flagstat {output.unique_bam:q} > {output.flagstats_unique:q}
		if [ -s {input.fasta:q} ]; then
			coverm contig -b {output.filtered_bam:q} -m mean length covered_bases count variance trimmed_mean rpkm -t {threads} -o {output.covstats:q}
			coverm contig -b {output.unique_sorted_bam:q} -m mean length covered_bases count variance trimmed_mean rpkm -t {threads} -o {output.covstats_unique:q}
			for basecov in {output.basecov:q}; do
				bedtools genomecov -dz -ibam {output.filtered_bam:q} | gzip -c > "$basecov"
			done
			for basecov in {output.unique_basecov:q}; do
				bedtools genomecov -dz -ibam {output.unique_sorted_bam:q} | gzip -c > "$basecov"
			done
		else
			printf 'Contig\tMean\tLength\tCovered Bases\tRead Count\tVariance\tTrimmed Mean\tRPKM\n' > {output.covstats:q}
			cp {output.covstats:q} {output.covstats_unique:q}
			for basecov in {output.basecov:q} {output.unique_basecov:q}; do
				printf '' | gzip -c > "$basecov"
			done
		fi
		"""

rule buildBowtieDB_reference:
	input:
		contaminants=REFERENCE_DIR+ "/" + REFERENCE + ".fasta",
	output:
		contigs_bt2=REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2",
	params:
		prefix=REFERENCE_DIR+ "/" + REFERENCE + "",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_reference/" + REFERENCE + "/tot.tsv"
	threads: 32
	shell:
		"""
		bowtie2-build {input.contaminants} {params.prefix} --threads {threads}
		"""

rule buildBowtieDB_reference_long:
	input:
		contaminants=REFERENCE_DIR+ "/" + REFERENCE + ".fasta",
	output:
		contigs_bt2=REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2l",
	params:
		prefix=REFERENCE_DIR+ "/" + REFERENCE + "",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_reference_long/" + REFERENCE + "/tot.tsv"
	threads: 32
	shell:
		"""
		bowtie2-build {input.contaminants} {params.prefix} --threads {threads} --large-index
		"""

def input_bowtie_reference(wildcards):
	if LONG_INDEX:
		return(REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2l")
	else:
		return(REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2")

rule mapReads_reference:
	input:
		contigs_bt2=input_bowtie_reference,
		forward_paired=lambda wc: mapping_read_path(wc, "forward"),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse"),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot.sam"),
		bam=(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_" + REFERENCE + "_{sample}.tot.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_filtered_" + REFERENCE + "_{sample}.tot.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_unique_" + REFERENCE + "_{sample}.tot.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique_sorted.bam"),
		unique_sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique_sorted.bam.bai"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique_covstats.txt",
		basecov=(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		unique_basecov=(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_tot_unique_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unique_flags="-F 2308" if ISOLATES else "-f 0x2",
		unpaired=mapping_orphan_option,
		prefix=REFERENCE_DIR+ "/" + REFERENCE + "",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapReads_reference/" + REFERENCE + "/sample={sample}.tsv"
	threads: 16
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --fast
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ {threads} -h {params.unique_flags} {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ {threads} -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ {threads} {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		for basecov in {output.basecov:q}; do
			bedtools genomecov -dz -ibam {output.filtered_bam} | gzip -c > "$basecov"
		done
		for basecov in {output.unique_basecov:q}; do
			bedtools genomecov -dz -ibam {output.unique_sorted_bam} | gzip -c > "$basecov"
		done
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""

rule mapReads_raw_reference:
	input:
		contigs_bt2=input_bowtie_reference,
		forward_paired=(dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config["forward_tag"]) + ".fastq.gz"),
		reverse_paired=(dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config["reverse_tag"]) + ".fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_" + REFERENCE + "_{sample}_raw.sam"),
		bam=dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_" + REFERENCE + "_{sample}_raw.bam",
		sorted_bam=temp(dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_" + REFERENCE + "_{sample}_raw_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_" + REFERENCE + "_{sample}_raw_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_" + REFERENCE + "_{sample}_raw_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_flagstats_" + REFERENCE + "_{sample}.raw.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_flagstats_filtered_" + REFERENCE + "_{sample}.raw.txt",
		covstats=dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_" + REFERENCE + "_{sample}_raw_covstats.txt",
		basecov=(dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW/bowtie2_" + REFERENCE + "_{sample}_raw_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		prefix=REFERENCE_DIR + "/" + REFERENCE,
		outdir=dirs_dict["MAPPING_DIR"] + "/REFERENCES/RAW",
	message:
		"Mapping raw paired-end reads to the reference"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapReads_raw_reference/" + REFERENCE + "/sample={sample}.tsv"
	threads: 16
	shell:
		"""
		mkdir -p {params.outdir}
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --fast
		samtools view -@ {threads} -bS {output.sam} > {output.bam}
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm -o {output.covstats}
		for basecov in {output.basecov:q}; do
			bedtools genomecov -dz -ibam {output.filtered_bam:q} | gzip -c > "$basecov"
		done
		"""

rule gene_Abundance:
	input:
		NR_bt2_150=expand(dirs_dict["ANNOTATION"] + "/predicted_genes_NR_95_85_150bp_tot.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"]),
		forward_paired=lambda wc: mapping_read_path(wc, "forward"),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse"),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/GENES/bowtie2_predicted_genes_NR_95_85_150bp_{sample}_tot.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/GENES/bowtie2_predicted_genes_NR_95_85_150bp_{sample}_tot.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/GENES/bowtie2_predicted_genes_NR_95_85_150bp_{sample}_tot_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/GENES/bowtie2_predicted_genes_NR_95_85_150bp_{sample}_tot_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/GENES/bowtie2_predicted_genes_NR_95_85_150bp_{sample}_tot_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/GENES/bowtie2_flagstats_predicted_genes_NR_95_85_150bp_{sample}.tot.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/GENES/bowtie2_flagstats_filtered_predicted_genes_NR_95_85_150bp_{sample}.tot.txt",
		covstats=dirs_dict["MAPPING_DIR"]+ "/GENES/bowtie2_predicted_genes_NR_95_85_150bp_{sample}_tot_covstats.txt",
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unpaired=mapping_orphan_option,
		prefix=dirs_dict["ANNOTATION"]+ "/predicted_genes_NR_95_85_150bp_tot",
	message:
		"Mapping reads to NR genes"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/gene_Abundance/sample={sample}.tsv"
	threads: 16
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --very-sensitive
		samtools view -@ {threads} -h -F 0x900 -bS {output.sam} > {output.bam}
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm -o {output.covstats}
		"""

rule mapReads_reference_sub:
	input:
		contigs_bt2=REFERENCE_DIR+ "/" + REFERENCE + ".1.bt2",
		forward_paired=lambda wc: mapping_read_path(wc, "forward", subset=True),
		reverse_paired=lambda wc: mapping_read_path(wc, "reverse", subset=True),
		unpaired=mapping_orphan_reads,
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub.sam"),
		bam=(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub.bam"),
		sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_" + REFERENCE + "_{sample}.sub.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_filtered_" + REFERENCE + "_{sample}.sub.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_flagstats_unique_" + REFERENCE + "_{sample}.sub.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique_sorted.bam"),
		unique_sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique_sorted.bam.bai"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique_covstats.txt",
		basecov=(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		unique_basecov=(dirs_dict["MAPPING_DIR"]+ "/REFERENCES/bowtie2_" + REFERENCE + "_{sample}_sub_unique_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		max_insert_bp=MAPPING_MAX_INSERT_BP,
		unique_flags="-F 2308" if ISOLATES else "-f 0x2",
		unpaired=mapping_orphan_option,
		prefix=REFERENCE_DIR+ "/" + REFERENCE + "",
	message:
		"Mapping reads to contigs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapReads_reference_sub/" + REFERENCE + "/sample={sample}.tsv"
	threads: 16
	shell:
		"""
		bowtie2 -X {params.max_insert_bp} -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} {params.unpaired} -S {output.sam} --threads {threads} --no-unal --fast
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ {threads} -h {params.unique_flags} {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ {threads} -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ {threads} {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		for basecov in {output.basecov:q}; do
			bedtools genomecov -dz -ibam {output.filtered_bam} | gzip -c > "$basecov"
		done
		for basecov in {output.unique_basecov:q}; do
			bedtools genomecov -dz -ibam {output.unique_sorted_bam} | gzip -c > "$basecov"
		done
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""

rule extract_mapped_reads:
	input:
		bam=(dirs_dict["MAPPING_DIR"]+ "/bowtie2_{sample}_tot.bam"),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
	output:
		mapped_reads=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_mapped_reads.txt"),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_mapped.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_mapped.fastq.gz"),
	message:
		"Extracting mapped reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/extract_mapped_reads/sample={sample}.tsv"
	threads: 16
	shell:
		"""
 		samtools view -F 4 {input.bam} | cut -f1 | sort | uniq > {output.mapped_reads}
		seqtk subseq {input.forward_paired} {output.mapped_reads} | gzip > {output.forward_paired}
		seqtk subseq {input.reverse_paired} {output.mapped_reads} | gzip > {output.reverse_paired}
		"""

def input_long_read_coverage_assembly(wildcards):
	if NANOPORE & NANOPORE_ONLY:
		return(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta")
	if NANOPORE & (not NANOPORE_ONLY) & PAIRED:
		return(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER +"_corrected_scaffolds_pilon.{sampling}.fasta")
	if PACBIO & PACBIO_ONLY:
		return(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta")
	if PACBIO & PACBIO_HYBRID:
		return(dirs_dict["ASSEMBLY_DIR"] + "/polypolish_{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta")
		
def input_long_read_coverage_reads(wildcards):
	if NANOPORE:
		return(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq.gz")
	if PACBIO:
		return(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_pacbio_clean.{sampling}.fastq.gz")

rule long_read_contig_coverage:
	input:
		assembly=input_long_read_coverage_assembly,
		reads=input_long_read_coverage_reads
	output:
		bam=temp(dirs_dict["MAPPING_DIR"] + "/{sample}_{sampling}_long_read_sorted.bam"),
		bai=temp(dirs_dict["MAPPING_DIR"] + "/{sample}_{sampling}_long_read_sorted.bam.bai"),
		coverage=dirs_dict["MAPPING_DIR"] + "/{sample}_{sampling}_long_read_contig_coverage.tsv"
	params:
		preset=lambda wildcards: "map-hifi" if PACBIO else "map-ont"
	message:
		"Calculating long-read per-contig coverage with CoverM"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/long_read_contig_coverage/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		minimap2 -ax {params.preset} -t {threads} {input.assembly} {input.reads} | samtools sort -@ {threads} -o {output.bam}
		samtools index {output.bam}
		coverm contig -b {output.bam} -m mean length covered_bases count variance trimmed_mean rpkm -o {output.coverage} -t {threads}
		"""
