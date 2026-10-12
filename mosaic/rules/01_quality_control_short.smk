#ruleorder: trim_adapters_quality_illumina_PE > trim_adapters_quality_illumina_SE
#ruleorder: listContaminants_PE > listContaminants_SE
#ruleorder: removeContaminants_PE > removeContaminants_SE
#ruleorder: subsampleReadsIllumina_PE > subsampleReadsIllumina_SE
#ruleorder: normalizeReads_PE > normalizeReads_SE
#ruleorder: postQualityCheckIlluminaPE > postQualityCheckIlluminaSE

# Prefer one counting job per paired-end sample/stage; retain single-file fallbacks.
ruleorder: countReads_raw > countReads_trimmed > countReads_noEuk > countReads_clean > countReads_norm > countReads_gz > countReads

rule download_SRA:
	input:
		sratoolkit="tools/sratoolkit.2.10.0-ubuntu64"
	output:
		forward_file=(dirs_dict["RAW_DATA_DIR"] + "/{SRA}_pass_1.fastq"),
		reverse_file=(dirs_dict["RAW_DATA_DIR"] + "/{SRA}_pass_2.fastq"),
	params:
		SRA_dir=dirs_dict["RAW_DATA_DIR"],
	message:
		"Downloading SRA run"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/download_SRA/SRA={SRA}.tsv"
#	threads: 1
	shell:
		"""
		{input.sratoolkit}/bin/fastq-dump --outdir {params.SRA_dir} --skip-technical --readids --read-filter pass \\
		--dumpbase --split-files --clip -N 0 -M 0 {wildcards.SRA}
		"""

rule countReads_raw:
	input:
		forward_paired=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq.gz",
		reverse_paired=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq.gz",
	output:
		forward_paired=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + "_read_count.txt",
		reverse_paired=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + "_read_count.txt",
	message:
		"Counting raw paired-end reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/countReads_raw/sample={sample}.tsv"
	wildcard_constraints:
		sample="|".join(re.escape(sample) for sample in SAMPLES) or "(?!)",
	resources:
		runtime_min=10,
		mem_mb=1000,
	shell:
		"""
		gzip -cd -- {input.forward_paired:q} | awk 'END {{print int(NR / 4)}}' > {output.forward_paired:q}
		gzip -cd -- {input.reverse_paired:q} | awk 'END {{print int(NR / 4)}}' > {output.reverse_paired:q}
		"""

rule countReads_trimmed:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.tot.fastq.gz",
	output:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_read_count.txt",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_read_count.txt",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.tot_read_count.txt",
	message:
		"Counting trimmed paired-end and unpaired reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/countReads_trimmed/sample={sample}.tsv"
	wildcard_constraints:
		sample="|".join(re.escape(sample) for sample in SAMPLES) or "(?!)",
	resources:
		runtime_min=15,
		mem_mb=1000,
	shell:
		"""
		gzip -cd -- {input.forward_paired:q} | awk 'END {{print int(NR / 4)}}' > {output.forward_paired:q}
		gzip -cd -- {input.reverse_paired:q} | awk 'END {{print int(NR / 4)}}' > {output.reverse_paired:q}
		gzip -cd -- {input.unpaired:q} | awk 'END {{print int(NR / 4)}}' > {output.unpaired:q}
		"""

rule countReads_noEuk:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_noEuk.tot.fastq",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_noEuk.tot.fastq",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_noEuk.tot.fastq",
	output:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_noEuk.tot_read_count.txt",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_noEuk.tot_read_count.txt",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_noEuk.tot_read_count.txt",
	message:
		"Counting paired-end and unpaired reads after eukaryote removal"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/countReads_noEuk/sample={sample}.tsv"
	wildcard_constraints:
		sample="|".join(re.escape(sample) for sample in SAMPLES) or "(?!)",
	resources:
		runtime_min=15,
		mem_mb=1000,
	shell:
		"""
		awk 'END {{print int(NR / 4)}}' {input.forward_paired:q} > {output.forward_paired:q}
		awk 'END {{print int(NR / 4)}}' {input.reverse_paired:q} > {output.reverse_paired:q}
		awk 'END {{print int(NR / 4)}}' {input.unpaired:q} > {output.unpaired:q}
		"""

rule countReads_clean:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz",
	output:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot_read_count.txt",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot_read_count.txt",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot_read_count.txt",
	message:
		"Counting clean paired-end and unpaired reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/countReads_clean/sample={sample}.tsv"
	wildcard_constraints:
		sample="|".join(re.escape(sample) for sample in SAMPLES) or "(?!)",
	resources:
		runtime_min=15,
		mem_mb=1000,
	shell:
		"""
		gzip -cd -- {input.forward_paired:q} | awk 'END {{print int(NR / 4)}}' > {output.forward_paired:q}
		gzip -cd -- {input.reverse_paired:q} | awk 'END {{print int(NR / 4)}}' > {output.reverse_paired:q}
		gzip -cd -- {input.unpaired:q} | awk 'END {{print int(NR / 4)}}' > {output.unpaired:q}
		"""

rule countReads_norm:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq.gz",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq.gz",
	output:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}_read_count.txt",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}_read_count.txt",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}_read_count.txt",
	message:
		"Counting normalized paired-end and unpaired reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/countReads_norm/sample={sample}__sampling={sampling}.tsv"
	wildcard_constraints:
		sample="|".join(re.escape(sample) for sample in SAMPLES) or "(?!)",
	resources:
		runtime_min=15,
		mem_mb=1000,
	shell:
		"""
		gzip -cd -- {input.forward_paired:q} | awk 'END {{print int(NR / 4)}}' > {output.forward_paired:q}
		gzip -cd -- {input.reverse_paired:q} | awk 'END {{print int(NR / 4)}}' > {output.reverse_paired:q}
		gzip -cd -- {input.unpaired:q} | awk 'END {{print int(NR / 4)}}' > {output.unpaired:q}
		"""

rule countReads_gz:
	input:
		fastq="{fastq_name}.fastq.gz",
	output:
		counts="{fastq_name}_read_count.txt",
	message:
		"Counting reads on fastq file"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/countReads_gz/fastq_name={fastq_name}.tsv"
	group:
		"read_counts_gz"
	resources:
		runtime_min= 5,
		mem_mb= 1000,
	shell:
		"""
		echo $(( $(zgrep -Ec "$" {input.fastq}) / 4 )) > {output.counts} 
		"""

rule countReads:
	input:
		fastq="{fastq_name}.fastq",
	output:
		counts="{fastq_name}_read_count.txt",
	message:
		"Counting reads on fastq file"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/countReads/fastq_name={fastq_name}.tsv"
	group:
		"read_counts"
	resources:
		runtime_min= 5,
		mem_mb= 1000,
	shell:
		"""
		echo $(( $(grep -Ec "$" {input.fastq}) / 4 )) > {output.counts} 
		"""

rule fastQC_pre:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"] + "/{fastq_name}.fastq.gz"
	output:
		html=temp(dirs_dict["RAW_DATA_DIR"] + "/{fastq_name}_fastqc.html"),
		zipped=(dirs_dict["RAW_DATA_DIR"] + "/{fastq_name}_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/fastQC_pre/fastq_name={fastq_name}.tsv"
	resources:
		runtime_min= 412,
		mem_mb= 500,
	shell:
		"""
		fastqc {input}
		"""

rule fastQC_post:
	input:
		raw_fastq=dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}.fastq.gz"
	output:
		html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}_fastqc.html"),
		zipped=(dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/fastQC_post/fastq_name={fastq_name}.tsv"
	resources:
		runtime_min= 412,
		mem_mb= 500,
	shell:
		"""
		fastqc {input}
		"""


rule superDeduper_pcr:
	input:
		forward_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq.gz",
		reverse_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq.gz",
	output:
		duplicate_stats=(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.log"),
		deduplicate=temp(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.out"),
	message:
		"Detect PCR duplicates"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/superDeduper_pcr/sample={sample}.tsv"
	resources:
		runtime_min= 30,
		mem_mb= 1000,
	shell:
		"""
		hts_SuperDeduper --log_freq 50000 -L {output.duplicate_stats} -1 {input.forward_file} -2 {input.reverse_file} > {output.deduplicate}
		"""

rule trim_adapters_quality_illumina_PE:
	input:
		forward_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq.gz",
		reverse_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq.gz",
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		forward_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq.gz"),
		reverse_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq.gz"),
		# Keep every trimmed orphan in one file; the two component files are temporary.
		merged_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.tot.fastq.gz",
		html=dirs_dict["QC_DIR"] + "/{sample}_fastp.html",
		json=dirs_dict["QC_DIR"] + "/{sample}_fastp.json",
	message:
		"Trimming Illumina adapters, poly-G tails and low-quality sequence with fastp"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/trim_adapters_quality_illumina_PE/sample={sample}.tsv"
	resources:
		runtime_min= 350,
		mem_mb= 4000,
	threads: 8
	log:
		dirs_dict["QC_DIR"] + "/{sample}_fastp.log"
	shell:
		"""
		fastp --thread {threads} --in1 {input.forward_file:q} --in2 {input.reverse_file:q} \
			--out1 {output.forward_paired:q} --out2 {output.reverse_paired:q} \
			--unpaired1 {output.forward_unpaired:q} --unpaired2 {output.reverse_unpaired:q} \
			--detect_adapter_for_pe --trim_poly_g --poly_g_min_len {config[fastp_poly_g_min_len]} \
			--cut_front --cut_front_window_size {config[fastp_cut_front_window_size]} \
			--cut_front_mean_quality {config[fastp_cut_front_mean_quality]} \
			--cut_right --cut_right_window_size {config[fastp_cut_right_window_size]} \
			--cut_right_mean_quality {config[fastp_cut_right_mean_quality]} \
			--length_required {config[fastp_length_required]} \
			--qualified_quality_phred {config[fastp_qualified_quality_phred]} \
			--unqualified_percent_limit {config[fastp_unqualified_percent_limit]} \
			--n_base_limit {config[fastp_n_base_limit]} \
			--html {output.html:q} --json {output.json:q} > {log:q} 2>&1
		cat {output.forward_unpaired:q} {output.reverse_unpaired:q} > {output.merged_unpaired:q}
		"""

rule sourmash_sketch_trim:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
	output:
		manysketch_csv=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_manysketch.csv"),
		sketch=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_sourmash.sig.zip"),
	params:
		sample="{sample}"
	message:
		"Building sketches with sourmash"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_sketch_trim/sample={sample}.tsv"
	threads: 8
	shell:
		"""
		echo name,read1,read2 > {output.manysketch_csv}
		echo {params.sample},{input.forward_paired},{input.reverse_paired} >> {output.manysketch_csv}
		sourmash scripts manysketch {output.manysketch_csv} -p k=31,k=51,abund,scaled=1000,DNA -o {output.sketch} -c {threads}
		"""

rule sourmash_sketch_clean_reads:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz",
	output:
		sketch=SOURMASH_CLEAN_DIR + "/SAMPLES/{sample}_sourmash.sig.zip",
	message:
		"Sketching all final clean paired and unpaired reads with Sourmash"
	conda:
		dirs_dict["ENVS_DIR"] + "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_sketch_clean_reads/sample={sample}.tsv"
	threads: 1
	shell:
		"""
		rm -f -- {output.sketch:q}
		sourmash sketch dna {input.forward_paired:q} {input.reverse_paired:q} {input.unpaired:q} \
			-p k=31,scaled=1000,abund --merge {wildcards.sample:q} -o {output.sketch:q}
		if [ ! -s {output.sketch:q} ]; then
			python - {output.sketch:q} {wildcards.sample:q} <<-'PYTHON'
			import sys
			import sourmash
			from sourmash.sourmash_args import SaveSignaturesToLocation
			minhash = sourmash.MinHash(n=0, ksize=31, scaled=1000, track_abundance=True)
			with SaveSignaturesToLocation(sys.argv[1]) as output:
			    output.add(sourmash.SourmashSignature(minhash, name=sys.argv[2]))
			PYTHON
		fi
		"""

rule sourmash_pool_clean_reads:
	input:
		sketches=expand(SOURMASH_CLEAN_DIR + "/SAMPLES/{sample}_sourmash.sig.zip", sample=SAMPLES),
	output:
		sketch=SOURMASH_CLEAN_DIR + "/POOLED/pooled_sourmash.sig.zip",
	message:
		"Pooling clean-read sketches from every sample, including the negative control"
	conda:
		dirs_dict["ENVS_DIR"] + "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_pool_clean_reads/tot.tsv"
	threads: 1
	shell:
		"""
		rm -f -- {output.sketch:q}
		if [ -n "{input.sketches}" ]; then
			sourmash sig merge {input.sketches:q} -k 31 --dna --name pooled -o {output.sketch:q}
		else
			python - {output.sketch:q} <<-'PYTHON'
			import sys
			import sourmash
			from sourmash.sourmash_args import SaveSignaturesToLocation
			minhash = sourmash.MinHash(n=0, ksize=31, scaled=1000, track_abundance=True)
			with SaveSignaturesToLocation(sys.argv[1]) as output:
			    output.add(sourmash.SourmashSignature(minhash, name="pooled"))
			PYTHON
		fi
		"""

rule sourmash_gather:
    benchmark:
        dirs_dict["BENCHMARKS"] + "/sourmash_gather/basedir={basedir}__sample={sample}.tsv"
    input:
        sketch="{basedir}/{sample}_sourmash.sig.zip",
        sourmash_rocksdb=config["sourmash_rocksdb"],
    output:
        gather="{basedir}/{sample}_gather_sourmash.csv",
    wildcard_constraints:
        basedir="(?:" + re.escape(dirs_dict["CLEAN_DATA_DIR"]) + "|" + re.escape(SOURMASH_CLEAN_DIR + "/SAMPLES") + "|" + re.escape(SOURMASH_CLEAN_DIR + "/POOLED") + ")",
    params:
        threshold_bp=lambda wildcards: int(config.get("sourmash_clean_min_shared_bp", 50000)) if wildcards.basedir != dirs_dict["CLEAN_DATA_DIR"] else 50000,
        abundances=lambda wildcards: wildcards.basedir != dirs_dict["CLEAN_DATA_DIR"],
    message:
        "Metagenome containment with Sourmash"
    conda:
        dirs_dict["ENVS_DIR"] + "/sourmash.yaml"
    threads: lambda wildcards: 1 if wildcards.basedir != dirs_dict["CLEAN_DATA_DIR"] else 8
    shell:
        """
        printf '' > {output.gather:q}
        sourmash_hashes=$(python - {input.sketch:q} <<-'PYTHON'
		import sys
		import sourmash
		print(sum(len(signature.minhash) for signature in sourmash.load_file_as_signatures(sys.argv[1])))
		PYTHON
        )
        if [ "$sourmash_hashes" -gt 0 ]; then
            if [ "{params.abundances}" = "True" ]; then
                sourmash gather {input.sketch:q} {input.sourmash_rocksdb:q} \
                    -k 31 --scaled 1000 --threshold-bp {params.threshold_bp} -o {output.gather:q} --create-empty-results
            else
                sourmash scripts fastmultigather \
                    {input.sketch:q} {input.sourmash_rocksdb:q} \
                    -c {threads} -o {output.gather:q} -t {params.threshold_bp} -s 1000
            fi
        fi
        if [ ! -s {output.gather:q} ]; then
            printf 'query_name,name,f_unique_weighted,f_unique_to_query,unique_intersect_bp,remaining_bp,query_md5,query_filename,query_bp,ksize,scaled,query_n_hashes\n' > {output.gather:q}
        fi
        """

rule sourmash_tax:
	input:
		gather="{basedir}/{sample}_gather_sourmash.csv",
		sourmash_tax=config['sourmash_tax'],
	output:
		kreport="{basedir}/{sample}_sourmash.kreport.txt",
		profile="{basedir}/{sample}_sourmash.summarized.csv",
		annotated_gather="{basedir}/{sample}_gather_sourmash.with-lineages.csv",
	wildcard_constraints:
		basedir="(?:" + re.escape(dirs_dict["CLEAN_DATA_DIR"]) + "|" + re.escape(SOURMASH_CLEAN_DIR + "/SAMPLES") + "|" + re.escape(SOURMASH_CLEAN_DIR + "/POOLED") + ")",
	params:
		sample="{sample}_sourmash",
		outdir="{basedir}",
	message:
		"Assigning taxonomy with sourmash tax"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_tax/basedir={basedir}__sample={sample}.tsv"
	threads: 1
	shell:
		"""
		if [ "$(wc -l < {input.gather:q})" -gt 1 ]; then
			# In v4, kreport and the CSV f_weighted_at_rank column retain abundance.
			# Explicit --use-abundances incorrectly rejects equal weighted/unweighted fractions.
			sourmash tax metagenome --gather-csv {input.gather:q} -t {input.sourmash_tax:q} -o {params.sample:q} \
				--v4 --output-format kreport csv_summary --rank species -f --output-dir {params.outdir:q}
			sourmash tax annotate --gather-csv {input.gather:q} -t {input.sourmash_tax:q} --output-dir {params.outdir:q}
		else
			printf '' > {output.kreport:q}
			printf 'query_name,rank,fraction,lineage,f_weighted_at_rank,bp_match_at_rank\n' > {output.profile:q}
			printf 'query_name,name,f_unique_weighted,f_unique_to_query,unique_intersect_bp,f_match_orig,average_abund,median_abund,scaled,lineage\n' > {output.annotated_gather:q}
		fi
		"""

rule contaminants_KRAKEN:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		merged_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.tot.fastq.gz"),
		kraken_db=(config['kraken_db']),
		kraken_tools=(config['kraken_tools']),
	output:
		kraken_output_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_tot.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_tot.csv"),
		kraken_domain=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_domains_tot.csv"),
		kraken_output_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_unpaired_tot.csv"),
		kraken_report_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_unpaired_tot.csv"),
	params:
		kraken_db=config['kraken_db'],
	message:
		"Assesing contamination with kraken2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/contaminants_KRAKEN/sample={sample}.tsv"
	threads: 16
	resources:
		runtime_min= 15,
		mem_mb= 18000,
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired}
		grep -P 'D\t' {output.kraken_report_paired} | sort -r > {output.kraken_domain}
		#UNPAIRED
		kraken2 --db {params.kraken_db} --threads {threads} {input.merged_unpaired}  \
			--output {output.kraken_output_unpaired} --report {output.kraken_report_unpaired}
		"""
		# python {input.kraken_tools}/combine_kreports.py \
		# 	-r {output.kraken_report_paired} {output.kraken_report_unpaired} \
		# 	-o {output.kraken_report_combined}
		# ""

rule contaminants_KRAKEN_microbial:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		kraken_db=(config['kraken_db_nt']),
	output:
		kraken_output_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_microbial.tot.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_microbial.tot.csv"),
	params:
		kraken_db=config['kraken_db_nt'],
	message:
		"Assesing contamination with kraken2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/contaminants_KRAKEN_microbial/sample={sample}.tsv"
	threads: 32
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired} \
			--report-minimizer-data
		"""

rule remove_euk:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		merged_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.tot.fastq.gz"),
		kraken_output_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_tot.csv"),
		kraken_output_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_unpaired_tot.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_tot.csv"),
		kraken_report_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_unpaired_tot.csv"),
		kraken_tools=(config['kraken_tools']),
	output:
		forward_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_noEuk.tot.fastq"),
		reverse_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_noEuk.tot.fastq"),
		unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_noEuk.tot.fastq"),
	message:
		"Removing eukaryotic reads with Kraken"
	params:
		# unclassified_name_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R#.tot.fastq",
		host_taxid=config["contaminants_taxid"] 
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/remove_euk/sample={sample}.tsv"
	threads: 4
	resources:
		mem_mb=40000,
		runtime_min= 1100,
	shell:
		"""
		python {input.kraken_tools}/extract_kraken_reads.py -k {input.kraken_output_paired} \
			-s1 {input.forward_paired} -s2 {input.reverse_paired} \
			-o {output.forward_paired} -o2 {output.reverse_paired} \
			--exclude -t {params.host_taxid} --include-children -r {input.kraken_report_paired} --fastq-output
		python {input.kraken_tools}/extract_kraken_reads.py -k {input.kraken_output_unpaired} \
			-s {input.merged_unpaired} -o {output.unpaired} --exclude -t {params.host_taxid} --include-children \
			-r {input.kraken_report_unpaired} --fastq-output
		"""

def remove_user_contaminants_forward(wildcards):
	if REMOVE_EUK:
		return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_forward_paired_noEuk.tot.fastq"
	return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_forward_paired.fastq.gz"

def remove_user_contaminants_reverse(wildcards):
	if REMOVE_EUK:
		return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_reverse_paired_noEuk.tot.fastq"
	return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_reverse_paired.fastq.gz"

def remove_user_contaminants_unpaired(wildcards):
	if REMOVE_EUK:
		return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_unpaired_noEuk.tot.fastq"
	return dirs_dict["CLEAN_DATA_DIR"] + f"/{wildcards.sample}_merged_unpaired.tot.fastq.gz"

rule remove_user_contaminants_PE:
	input:
		forward_paired=remove_user_contaminants_forward,
		reverse_paired=remove_user_contaminants_reverse,
		unpaired=remove_user_contaminants_unpaired,
		contaminants_fasta=expand(dirs_dict["CONTAMINANTS_DIR_DB"] +"/{contaminants}.fasta",contaminants=CONTAMINANTS) if not ISOLATES else [],
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz"),
		phix_contaminants_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{sample}_contaminants.fasta",
		stats=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_contaminant_stats_bbduk.tot.txt"),
	params:
		has_contaminants=bool(CONTAMINANTS) and not ISOLATES,
	message:
		"Removing configured contaminants or passing reads through unchanged"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/remove_user_contaminants_PE/sample={sample}.tsv"
	threads: 4
	resources:
		mem_mb=40000,
		runtime_min= 15,
	shell:
		"""
		if [ "{params.has_contaminants}" = "True" ]; then
			cat {input.contaminants_fasta:q} > {output.phix_contaminants_fasta:q}
			bbduk.sh -Xmx{resources.mem_mb}m in1={input.forward_paired:q} in2={input.reverse_paired:q} out1={output.forward_paired:q} out2={output.reverse_paired:q} \
				ref={output.phix_contaminants_fasta:q} k=31 hdist=1 threads={threads} stats={output.stats:q}
			bbduk.sh -Xmx{resources.mem_mb}m in={input.unpaired:q} out={output.unpaired:q} ref={output.phix_contaminants_fasta:q} k=31 hdist=1 threads={threads}
		else
			: > {output.phix_contaminants_fasta:q}
			if [[ {input.forward_paired:q} == *.gz ]]; then
				ln -sfnr -- {input.forward_paired:q} {output.forward_paired:q}
			else
				gzip -c -- {input.forward_paired:q} > {output.forward_paired:q}
			fi
			if [[ {input.reverse_paired:q} == *.gz ]]; then
				ln -sfnr -- {input.reverse_paired:q} {output.reverse_paired:q}
			else
				gzip -c -- {input.reverse_paired:q} > {output.reverse_paired:q}
			fi
			if [[ {input.unpaired:q} == *.gz ]]; then
				ln -sfnr -- {input.unpaired:q} {output.unpaired:q}
			else
				gzip -c -- {input.unpaired:q} > {output.unpaired:q}
			fi
			printf 'No contaminant references configured; read filtering skipped.\n' > {output.stats:q}
		fi
		"""

rule contaminants_KRAKEN_clean:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz",
		kraken_db=(config['kraken_db']),
		kraken_tools=(config['kraken_tools']),
	output:
		kraken_output_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_clean_tot.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_clean_tot.csv"),
	params:
		kraken_db=config['kraken_db'],
	message:
		"Assesing taxonomy with kraken2 on clean reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/contaminants_KRAKEN_clean/sample={sample}.tsv"
	priority: 1
	threads: 8
	resources:
		runtime_min= 15,
		mem_mb= 18000,
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired} \
			--report-minimizer-data
		"""

rule preMultiQC:
	input:
		#html=expand(dirs_dict["RAW_DATA_DIR"]+"/{sample}_{reads}_fastqc.html", sample=SAMPLES, reads=READ_TYPES),
		zipped=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{reads}_fastqc.zip", sample=SAMPLES, reads=READ_TYPES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
		multiqc_txt=dirs_dict["QC_DIR"]+ "/preQC_illumina_report_data/multiqc_fastqc.txt",
	params:
		fastqc_dir=dirs_dict["RAW_DATA_DIR"],
		html_name="preQC_illumina_report.html",
		multiqc_dir=dirs_dict["QC_DIR"],
	message:
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/preMultiQC/tot.tsv"
	resources:
		runtime_min= 5,
		mem_mb= 4000,
	shell:
		"""
		multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule postMultiQC:
	input:
		# html_forward=expand(dirs_dict["CLEAN_DATA_DIR"]  + "/{sample}_forward_paired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot_fastqc.zip", sample=SAMPLES),
		# html_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot_fastqc.zip", sample=SAMPLES),
		# html_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot_fastqc.html", sample=SAMPLES),
		zipped_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"]  + "/{sample}_unpaired_clean.tot_fastqc.zip", sample=SAMPLES),
		fastp_json=expand(dirs_dict["QC_DIR"] + "/{sample}_fastp.json", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/postQC_illumina_report.html",
		multiqc_txt=dirs_dict["QC_DIR"]+ "/postQC_illumina_report_data/multiqc_fastqc.txt",
	params:
		html_name="postQC_illumina_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/postMultiQC/tot.tsv"
	priority: 1
	resources:
		runtime_min= 5,
		mem_mb= 4000,
	shell:
		"""
		multiqc -f {input:q} -o {params.multiqc_dir:q} -n {params.html_name:q}
		"""

rule prekrakenMultiQC:
	input:
		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_tot.csv", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/pre_decontamination_kraken_multiqc_report.html"
		# 		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="pre_decontamination_kraken_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report kraken pre"
	priority: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/prekrakenMultiQC/tot.tsv"
	resources:
		runtime_min= 5,
		mem_mb= 4000,
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule postkrakenMultiQC:
	input:
		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_clean_tot.csv", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/post_decontamination_kraken_multiqc_report.html"
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="post_decontamination_kraken_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report kraken post"
	priority: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/postkrakenMultiQC/tot.tsv"
	resources:
		runtime_min= 5,
		mem_mb= 4000,
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule krakenMicrobialMultiQC:
	input:
		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_microbial.tot.csv", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/microbial_kraken_multiqc_report.html"
		# 		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="microbial_kraken_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report kraken pre"
	priority: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/krakenMicrobialMultiQC/tot.tsv"
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule normalizeReads_PE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.{sampling}.fastq.gz",
	output:
		forward_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.{sampling}.fastq.gz"),
		reverse_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.{sampling}.fastq.gz"),
		unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq.gz"),
		histogram_pre=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_histogram_pre.{sampling}.txt"),
		histogram_post=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_histogram_post.{sampling}.txt"),
		peaks=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_count_peaks.{sampling}.txt"),
	message:
		"Normalizing reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/normalizeReads_PE/sample={sample}__sampling={sampling}.tsv"
	params:
		min_depth=config['min_norm'],
		max_depth=config['max_norm'],
		heap_mb=MEMORY_ECORR
	threads: 16
	priority: 1
	wildcard_constraints:
		sampling="tot|sub"  
	resources:
		mem_mb=BBTOOLS_MEM_MB,
		runtime_min= 125,
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{params.heap_mb}m in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
			target={params.max_depth} mindepth={params.min_depth} t={threads} khist={output.histogram_pre} peaks={output.peaks} khistout={output.histogram_post}
		#unpaired
		bbnorm.sh -Xmx{params.heap_mb}m in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth} threads={threads}
		"""

rule concatenate_subassembly:
	input:
		forward_paired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",sample=SAMPLES),
		reverse_paired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",sample=SAMPLES),
		unpaired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz",sample=SAMPLES),
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/ALL_forward_paired_clean.tot.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/ALL_reverse_paired_clean.tot.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/ALL_unpaired_clean.tot.fastq.gz",
	message:
		"Concatenating clean reads for cross assembly"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/concatenate_subassembly/tot.tsv"
	shell:
		"""
		cat {input.forward_paired} > {output.forward_paired}
		cat {input.reverse_paired} > {output.reverse_paired}
		cat {input.unpaired} > {output.unpaired}
		"""
rule kmer_rarefraction:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		histogram=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.{sampling}.csv"),
	message:
		"Counting unique reads with BBtools"
	params:
		heap_mb=MEMORY_ECORR
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/kmer_rarefraction/sample={sample}__sampling={sampling}.tsv"
	threads: 1
	resources:
		mem_mb=BBTOOLS_MEM_MB,
		runtime_min= 412,
	shell:
		"""
		bbcountunique.sh -Xmx{params.heap_mb}m in1={input.forward_paired} in2={input.reverse_paired} out={output.histogram} interval={config[kmer_window]}
		"""
