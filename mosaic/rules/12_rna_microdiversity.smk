# Additional RNA discovery and per-sample variation against the shared filtered vOTU catalogue.

rule rna_assemble_spades:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",
	output:
		fasta=temp(RNA_DIR + "/{sample}_rnaviralspades.unrenamed.fasta"),
	params:
		work_prefix=RNA_DIR + "/{sample}/spades_work_",
		mem_gb=RNA_SPADES_MEM_GB,
	message:
		"Assembling paired-end reads with RNAviralSPAdes"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rna_assemble_spades/sample={sample}.tsv"
	threads: int(config.get("rna_assembly_threads", 16))
	resources:
		mem_mb=RNA_SPADES_MEM_MB,
	log:
		RNA_DIR + "/{sample}/rnaviralspades.log"
	shell:
		r"""
		work_dir=$(mktemp -d {params.work_prefix:q}XXXXXX)
		trap 'rm -rf -- "$work_dir"' EXIT
		(
			spades.py --rnaviral -1 {input.forward_paired:q} -2 {input.reverse_paired:q} \
				-o "$work_dir/spades" -t {threads} -m {params.mem_gb}
			cp "$work_dir/spades/contigs.fasta" {output.fasta:q}
		) > {log:q} 2>&1
		"""

rule rna_assemble_megahit:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",
	output:
		fasta=temp(RNA_DIR + "/{sample}_megahit.unrenamed.fasta"),
	params:
		work_prefix=RNA_DIR + "/{sample}/megahit_work_",
		mem_bytes=RNA_MEGAHIT_MEM_GB * 1024 ** 3,
		min_length=int(config.get("rna_min_contig_length", 500)),
	message:
		"Assembling paired-end reads with MEGAHIT"
	conda:
		dirs_dict["ENVS_DIR"] + "/rna_megahit.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rna_assemble_megahit/sample={sample}.tsv"
	threads: int(config.get("rna_assembly_threads", 16))
	resources:
		mem_mb=RNA_MEGAHIT_MEM_MB,
	log:
		RNA_DIR + "/{sample}/megahit.log"
	shell:
		r"""
		work_dir=$(mktemp -d {params.work_prefix:q}XXXXXX)
		trap 'rm -rf -- "$work_dir"' EXIT
		(
			megahit -1 {input.forward_paired:q} -2 {input.reverse_paired:q} -o "$work_dir/megahit" \
				-t {threads} -m {params.mem_bytes} --min-contig-len {params.min_length}
			cp "$work_dir/megahit/final.contigs.fa" {output.fasta:q}
		) > {log:q} 2>&1
		"""

rule rna_assemble_trinity:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",
	output:
		fasta=temp(RNA_DIR + "/{sample}_trinity.unrenamed.fasta"),
	params:
		assembly_dir=lambda wildcards: os.path.abspath(RNA_DIR + "/" + wildcards.sample + "/trinity_out"),
		assembled_fasta=lambda wildcards: os.path.abspath(RNA_DIR + "/" + wildcards.sample + "/trinity_out.Trinity.fasta"),
		mem_gb=RNA_TRINITY_MEM_GB,
		bfly_heap_gb=RNA_TRINITY_BFLY_HEAP_GB,
		min_length=int(config.get("rna_min_contig_length", 500)),
		strand_flag="--SS_lib_type " + config["rna_trinity_strandedness"] if config.get("rna_trinity_strandedness", "") else "",
	message:
		"Assembling paired-end reads with Trinity"
	conda:
		dirs_dict["ENVS_DIR"] + "/rna_trinity.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rna_assemble_trinity/sample={sample}.tsv"
	threads: int(config.get("rna_trinity_threads", 8))
	resources:
		mem_mb=lambda wc, threads: memory_reservation("rna_trinity_mem_mb",
			(RNA_TRINITY_MEM_GB + threads * RNA_TRINITY_BFLY_HEAP_GB) * 1024, 224 * 1024, 280 * 1024),
	log:
		RNA_DIR + "/{sample}/trinity.log"
	shell:
		r"""
		(
			Trinity --seqType fq --left {input.forward_paired:q} --right {input.reverse_paired:q} \
				--output {params.assembly_dir:q} --CPU {threads} --max_memory {params.mem_gb}G \
				--bflyHeapSpaceMax {params.bfly_heap_gb}G \
				--no_normalize_reads --no_salmon --full_cleanup --min_contig_length {params.min_length} {params.strand_flag}
			cp {params.assembled_fasta:q} {output.fasta:q}
			test -s {output.fasta:q}
			rm -rf -- {params.assembly_dir:q}
		) > {log:q} 2>&1
		"""

rule rna_genomad_assembler:
	input:
		fasta=RNA_DIR + "/{sample}_{assembler}.fasta",
		db=config["genomad_db"],
	output:
		outdir=directory(RNA_DIR + "/{sample}/{assembler}_genomad"),
	params:
		empty_summary=RNA_DIR + "/{sample}/{assembler}_genomad/{sample}_{assembler}_summary/{sample}_{assembler}_virus_summary.tsv",
	message:
		"Identifying viral contigs in the named RNA assembly with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rna_genomad_assembler/sample={sample}__assembler={assembler}.tsv"
	threads: 8
	wildcard_constraints:
		assembler="rnaviralspades|megahit|trinity"
	log:
		RNA_DIR + "/{sample}/{assembler}_genomad.log"
	shell:
		r"""
		if [ -s {input.fasta:q} ]; then
			genomad end-to-end --restart --cleanup --splits 8 -t {threads} \
				{input.fasta:q} {output.outdir:q} {input.db:q} --relaxed > {log:q} 2>&1
		else
			mkdir -p "$(dirname {params.empty_summary:q})"
			printf 'seq_name\tlength\ttopology\tvirus_score\ttaxonomy\n' > {params.empty_summary:q}
			printf 'No assembled contigs to classify.\n' > {log:q}
		fi
		"""

rule rna_combine_assemblies:
	input:
		fastas=lambda wc: expand(RNA_DIR + "/{sample}_{assembler}.fasta", sample=wc.sample, assembler=RNA_ASSEMBLERS),
	output:
		fasta=RNA_DIR + "/{sample}_combined_derreplicated.fasta",
		provenance=RNA_DIR + "/{sample}_assembly_provenance.tsv",
	params:
		names=RNA_ASSEMBLERS,
		min_length=int(config.get("rna_min_contig_length", 500)),
	message:
		"Combining RNA assemblies and recording exact-sequence provenance"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rna_combine_assemblies/sample={sample}.tsv"
	threads: 1
	run:
		import csv
		import hashlib

		COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")


		def records(path):
			name, chunks = None, []
			with open(path) as handle:
				for line in handle:
					line = line.strip()
					if not line:
						continue
					if line.startswith(">"):
						if name is not None:
							yield name, "".join(chunks).upper().replace("U", "T")
						name, chunks = line[1:].split()[0], []
					elif name is None:
						raise ValueError(f"Invalid FASTA: {path}")
					else:
						chunks.append(line)
			if name is not None:
				yield name, "".join(chunks).upper().replace("U", "T")


		if len(input.fastas) != len(params.names):
			raise ValueError("Each input FASTA must have a provenance source name")
		seen = {}
		with open(output.fasta, "w") as fasta, open(output.provenance, "w") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow(["representative", "source", "contig_id", "length", "status"])
			for path, source in zip(input.fastas, params.names):
				for name, seq in records(path):
					if not seq or len(seq) < params.min_length:
						writer.writerow(["", source, name, len(seq), "below_min_length"])
						continue
					if set(seq) - set("ACGTRYMKSWBDHVN"):
						raise ValueError(f"Invalid nucleotide sequence in {path}: {name}")
					canonical = min(seq, seq.translate(COMPLEMENT)[::-1])
					key = hashlib.sha256(canonical.encode()).digest()
					if key in seen:
						writer.writerow([seen[key], source, name, len(seq), "exact_duplicate"])
						continue
					identifier = name
					seen[key] = identifier
					fasta.write(f">{identifier}\n{seq}\n")
					writer.writerow([identifier, source, name, len(seq), "retained"])

rule assemblyStats_RNA:
	input:
		fastas=expand(RNA_DIR + "/{sample}_{assembler}.fasta", sample=SAMPLES, assembler=RNA_ASSEMBLERS + ["combined_derreplicated"]),
	output:
		quast_report_dir=directory(RNA_DIR + "/statistics_quast_tot"),
		quast_txt=RNA_DIR + "/assembly_quast_report.tot.txt",
		quast_tsv=RNA_DIR + "/statistics_quast_tot/transposed_report.tsv",
	message:
		"Creating QUAST statistics for each RNA assembler and the combined assemblies"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/assemblyStats_RNA/tot.tsv"
	log:
		RNA_DIR + "/quast.tot.log"
	threads: 4
	resources:
		mem_mb=QUAST_MEM_MB
	shell:
		"""
		rna_quast_fastas=()
		rna_quast_labels=""
		for rna_quast_fasta in {input.fastas:q}; do
			if [ -s "$rna_quast_fasta" ]; then
				rna_quast_fastas+=("$rna_quast_fasta")
				rna_quast_label=$(basename "$rna_quast_fasta" .fasta)
				rna_quast_labels+="${{rna_quast_label}},"
			fi
		done
		if [ "${{#rna_quast_fastas[@]}}" -gt 0 ]; then
			# Assess every saved contig, matching the RNA notebook's length-filtered FASTA statistics.
			quast.py "${{rna_quast_fastas[@]}}" -o {output.quast_report_dir:q} \
				--labels "${{rna_quast_labels%,}}" --min-contig 1 --threads {threads} > {log:q} 2>&1
			cp {output.quast_report_dir:q}/report.txt {output.quast_txt:q}
		else
			mkdir -p {output.quast_report_dir:q}
			printf 'Assembly\n' > {output.quast_tsv:q}
			printf 'No RNA contigs to assess.\n' > {output.quast_report_dir:q}/report.txt
			cp {output.quast_report_dir:q}/report.txt {output.quast_txt:q}
			printf 'No RNA contigs to assess.\n' > {log:q}
		fi
		"""

rule rna_identify_candidates:
	input:
		fasta=RNA_DIR + "/{sample}_combined_derreplicated.fasta",
		db=config["virSorter_db"],
	output:
		fasta=RNA_DIR + "/{sample}/virsorter/final-viral-combined.fa",
		scores=RNA_DIR + "/{sample}/virsorter/final-viral-score.tsv",
	params:
		outdir=RNA_DIR + "/{sample}/virsorter",
		groups=config.get("rna_viral_groups", "RNA"),
	message:
		"Identifying RNA-branch viral candidates with VirSorter2"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rna_identify_candidates/sample={sample}.tsv"
	threads: 8
	log:
		RNA_DIR + "/{sample}/virsorter.log"
	shell:
		r"""
		if [ -s {input.fasta:q} ]; then
			virsorter run -w {params.outdir:q} -i {input.fasta:q} -j {threads} --db-dir {input.db:q} \
				--include-groups {params.groups:q} --seqname-suffix-off --provirus-off \
				--keep-original-seq --min-length 0 all > {log:q} 2>&1
		else
			mkdir -p {params.outdir:q}
			printf '' > {output.fasta:q}
			printf 'seqname\tmax_score\tmax_score_group\n' > {output.scores:q}
			printf 'No contigs passed the assembly length filter.\n' > {log:q}
		fi
		test -f {output.fasta:q}
		test -f {output.scores:q}
		find {params.outdir:q} -mindepth 1 -maxdepth 1 \
			! -name 'final-viral-combined.fa' ! -name 'final-viral-score.tsv' \
			! -name 'final-viral-boundary.tsv' -exec rm -rf -- {{}} +
		rm -f -- {log:q}
		"""

rule rna_viral_positive_fasta:
	input:
		fasta=RNA_DIR + "/{sample}/virsorter/final-viral-combined.fa",
	output:
		fasta=RNA_DIR + "/{sample}_virsorter2_RNA_viral_positive.fasta",
	params:
		target="{sample}/virsorter/final-viral-combined.fa",
	message:
		"Linking RNA VirSorter2-positive contigs into the RNA assembly directory"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rna_viral_positive_fasta/sample={sample}.tsv"
	threads: 1
	shell:
		"""
		ln -sfn {params.target:q} {output.fasta:q}
		"""

rule rna_checkv:
	input:
		fasta=RNA_DIR + "/{sample}_virsorter2_RNA_viral_positive.fasta",
		db=config["checkv_db"],
	output:
		quality=RNA_DIR + "/{sample}/checkv/quality_summary.tsv",
		completeness=RNA_DIR + "/{sample}/checkv/completeness.tsv",
		contamination=RNA_DIR + "/{sample}/checkv/contamination.tsv",
	params:
		outdir=RNA_DIR + "/{sample}/checkv",
	message:
		"Assessing RNA-branch viral candidates with CheckV"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/rna_checkv/sample={sample}.tsv"
	threads: 4
	resources:
		mem_mb=8000,
	log:
		RNA_DIR + "/{sample}/checkv.log"
	shell:
		r"""
		mkdir -p {params.outdir:q}
		(
			if [ -s {input.fasta:q} ]; then
				checkv contamination {input.fasta:q} {params.outdir:q} -t {threads} -d {input.db:q}
				checkv completeness {input.fasta:q} {params.outdir:q} -t {threads} -d {input.db:q}
				checkv complete_genomes {input.fasta:q} {params.outdir:q}
				checkv quality_summary {input.fasta:q} {params.outdir:q}
			else
				printf 'contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n' > {output.quality:q}
				printf '' > {output.completeness:q}
				printf '' > {output.contamination:q}
				printf 'No RNA-branch viral candidates to assess.\n'
			fi
		) > {log:q} 2>&1
		"""

rule microdiversity_reference:
	input:
		fasta=dirs_dict["vOUT_DIR"] + "/filtered_" + REPRESENTATIVE_CONTIGS_BASE + ".tot.fasta",
	output:
		fasta=MICRO_DIR + "/reference/catalogue.fasta",
		fai=MICRO_DIR + "/reference/catalogue.fasta.fai",
		bowtie_index=temp(directory(MICRO_DIR + "/reference/bowtie2")),
	params:
		prefix=MICRO_DIR + "/reference/bowtie2/catalogue",
	message:
		"Indexing the shared filtered vOTU catalogue for microdiversity"
	conda:
		dirs_dict["ENVS_DIR"] + "/microdiversity.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/microdiversity_reference/tot.tsv"
	threads: 8
	shell:
		r"""
		test -s {input.fasta:q} || {{ echo 'The final vOTU catalogue is empty; no microdiversity reference.' >&2; exit 1; }}
		cp {input.fasta:q} {output.fasta:q}
		samtools faidx {output.fasta:q}
		mkdir -p {output.bowtie_index:q}
		bowtie2-build --threads {threads} {output.fasta:q} {params.prefix:q}
		"""

rule microdiversity_map:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",
		bowtie_index=MICRO_DIR + "/reference/bowtie2",
	output:
		bam=MICRO_DIR + "/{sample}/aligned.sorted.bam",
		bai=MICRO_DIR + "/{sample}/aligned.sorted.bam.bai",
		flagstat=MICRO_DIR + "/{sample}/flagstat.txt",
	params:
		prefix=MICRO_DIR + "/reference/bowtie2/catalogue",
	message:
		"Mapping paired-end reads to the shared catalogue for microdiversity"
	conda:
		dirs_dict["ENVS_DIR"] + "/microdiversity.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/microdiversity_map/sample={sample}.tsv"
	threads: int(config.get("microdiversity_threads", 8))
	resources:
		mem_mb=8000,
	log:
		MICRO_DIR + "/{sample}/bowtie2.log"
	shell:
		r"""
		bowtie2 --very-sensitive --seed 1 -x {params.prefix:q} -1 {input.forward_paired:q} -2 {input.reverse_paired:q} \
			-p {threads} 2> {log:q} | samtools sort -@ {threads} -m 256M -o {output.bam:q} -
		samtools index {output.bam:q}
		samtools flagstat {output.bam:q} > {output.flagstat:q}
		"""

rule microdiversity_statistics:
	input:
		bam=MICRO_DIR + "/{sample}/aligned.sorted.bam",
		bai=MICRO_DIR + "/{sample}/aligned.sorted.bam.bai",
		fasta=MICRO_DIR + "/reference/catalogue.fasta",
		fai=MICRO_DIR + "/reference/catalogue.fasta.fai",
	output:
		positions=MICRO_DIR + "/{sample}/positions.tsv.gz",
		variants=MICRO_DIR + "/{sample}/snv_candidates.tsv",
		summary=MICRO_DIR + "/{sample}/summary.tsv",
	params:
		min_depth=int(config.get("microdiversity_min_depth", 100)),
		baseq=int(config.get("microdiversity_min_base_quality", 30)),
		mapq=int(config.get("microdiversity_min_mapping_quality", 30)),
		min_af=float(config.get("microdiversity_min_alt_frequency", 0.03)),
		min_count=int(config.get("microdiversity_min_alt_count", 5)),
		max_depth=int(config.get("microdiversity_max_depth", 100000)),
	message:
		"Calculating callable-site nucleotide diversity against the final vOTU catalogue"
	conda:
		dirs_dict["ENVS_DIR"] + "/microdiversity.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/microdiversity_statistics/sample={sample}.tsv"
	threads: 1
	resources:
		mem_mb=4000,
	shell:
		r"""
		python - {input.bam:q} {input.fasta:q} {output.positions:q} {output.variants:q} {output.summary:q} \
			{wildcards.sample:q} {params.min_depth} {params.baseq} {params.mapq} {params.min_count} {params.max_depth} {params.min_af} <<-'PYTHON'
		import csv
		import gzip
		import math
		from collections import Counter


		def diversity(counts):
		    depth = sum(counts)
		    if depth < 2:
		        return float("nan"), float("nan")
		    frequencies = [count / depth for count in counts if count]
		    pi = (depth * depth - sum(count * count for count in counts)) / (depth * (depth - 1))
		    entropy = -sum(p * math.log2(p) for p in frequencies)
		    return pi, entropy


		import sys
		import pysam

		bam_path, fasta_path, positions_path, variants_path, summary_path, sample = sys.argv[1:7]
		min_depth, baseq, mapq, min_count, max_depth = map(int, sys.argv[7:12])
		min_af = float(sys.argv[12])

		if not (2 <= min_depth < max_depth and min_count >= 1 and 0 < min_af <= 0.5):
		    raise ValueError("Require 2 <= min_depth < max_depth, min_alt_count >= 1, and 0 < min_alt_frequency <= 0.5")
		if min(baseq, mapq) < 0:
		    raise ValueError("Quality thresholds cannot be negative")
		bases = "ACGT"
		with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(fasta_path) as ref, \
		        gzip.open(positions_path, "wt") as pos_handle, \
		        open(variants_path, "w") as var_handle, open(summary_path, "w") as summary_handle:
		    positions = csv.writer(pos_handle, delimiter="\t", lineterminator="\n")
		    variants = csv.writer(var_handle, delimiter="\t", lineterminator="\n")
		    summary = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
		    positions.writerow(["sample", "contig", "position", "reference", "depth", "A", "C", "G", "T",
		                        "A_reverse", "C_reverse", "G_reverse", "T_reverse", "callable", "depth_capped", "pi", "entropy_bits"])
		    variants.writerow(["sample", "contig", "position", "reference", "alternate", "depth", "alt_count",
		                       "alt_frequency", "alt_forward", "alt_reverse", "type"])
		    summary.writerow(["sample", "contig", "length", "acgt_reference_sites", "covered_sites", "callable_sites",
		                      "callable_fraction", "depth_capped_sites", "mean_depth", "polymorphic_sites",
		                      "pi_callable", "entropy_callable_bits", "min_depth", "min_base_quality",
		                      "min_mapping_quality", "min_alt_count", "min_alt_frequency", "max_depth"])
		    for contig, length in zip(ref.references, ref.lengths):
		        sequence = ref.fetch(contig).upper()
		        acgt_sites = sum(base in bases for base in sequence)
		        covered = callable_sites = capped_sites = depth_sum = polymorphic = 0
		        pi_sum = entropy_sum = 0.0
		        for col in bam.pileup(contig, 0, length, truncate=True, stepper="samtools", fastafile=ref,
		                              min_base_quality=baseq, min_mapping_quality=mapq,
		                              ignore_overlaps=True, ignore_orphans=True, compute_baq=False,
		                              flag_filter=3844, max_depth=max_depth):
		            pos = col.reference_pos
		            counts, reverse = Counter(), Counter()
		            # Access pileups explicitly: count only A/C/G/T, excluding deletions and skipped bases.
		            for read in col.pileups:
		                if read.is_del or read.is_refskip:
		                    continue
		                alignment = read.alignment
		                base = alignment.query_sequence[read.query_position].upper()
		                if base in bases:
		                    counts[base] += 1
		                    if alignment.is_reverse:
		                        reverse[base] += 1
		            depth = sum(counts.values())
		            depth_sum += depth
		            covered += depth > 0
		            capped = col.nsegments >= max_depth
		            capped_sites += capped
		            callable_site = depth >= min_depth and sequence[pos] in bases and not capped
		            pi = entropy = "NA"
		            if callable_site:
		                callable_sites += 1
		                pi, entropy = diversity([counts[b] for b in bases])
		                pi_sum += pi
		                entropy_sum += entropy
		                ordered = sorted(counts.values(), reverse=True)
		                is_polymorphic = len(ordered) > 1 and ordered[1] >= min_count and ordered[1] / depth >= min_af
		                polymorphic += is_polymorphic
		                for base in bases:
		                    n = counts[base]
		                    if base != sequence[pos] and n >= min_count and n / depth >= min_af:
		                        variants.writerow([sample, contig, pos + 1, sequence[pos], base, depth, n,
		                                           n / depth, n - reverse[base], reverse[base],
		                                           "polymorphic" if is_polymorphic else "reference_difference"])
		            positions.writerow([sample, contig, pos + 1, sequence[pos], depth,
		                                *[counts[b] for b in bases], *[reverse[b] for b in bases],
		                                int(callable_site), int(capped), pi, entropy])
		        summary.writerow([sample, contig, length, acgt_sites, covered, callable_sites,
		                          callable_sites / acgt_sites if acgt_sites else "NA", capped_sites,
		                          depth_sum / length if length else "NA", polymorphic,
		                          pi_sum / callable_sites if callable_sites else "NA",
		                          entropy_sum / callable_sites if callable_sites else "NA",
		                          min_depth, baseq, mapq, min_count, min_af, max_depth])
		PYTHON
		"""
