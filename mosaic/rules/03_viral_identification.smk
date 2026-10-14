rule satellite_finder:
	input:
		genomad_outdir=dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}/",
	output:
		satellite_finder_outdir=(dirs_dict["VIRAL_DIR"] + "/{sample}_{sampling}_satellite_finder_{model}/best_solution_summary.tsv"),
		faa_temp=(dirs_dict["VIRAL_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}_proteins.faa_fixed_{model}"),
		faa_idx_temp=temp(dirs_dict["VIRAL_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}_proteins.faa_fixed_{model}.idx"),
	params:
		faa=dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}/{sample}_spades_filtered_scaffolds.{sampling}_annotate/{sample}_spades_filtered_scaffolds.{sampling}_proteins.faa",
		model="{model}",
		satellite_finder_dir=os.path.join(workflow.basedir, "tools"),
		satellite_finder_outdir=(dirs_dict["VIRAL_DIR"] + "/{sample}_{sampling}_satellite_finder_{model}/"),
	message:
		"Identifying viral satellites with satellite_finder"
	conda:
		dirs_dict["ENVS_DIR"] + "/satellite_finder.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/satellite_finder/model={model}__sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		python scripts/process_fasta_satellite_finder.py {params.faa} {output.faa_temp}
		apptainer run -H ${{HOME}} {params.satellite_finder_dir}/fixed_satellite_finder.sif --db-type gembase --models {params.model} --sequence-db {output.faa_temp} -w {threads} -o {params.satellite_finder_outdir} --mute
		"""

rule satellite_finder_get_fasta:
	input:
		scaffolds_spades=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
		summary=(dirs_dict["VIRAL_DIR"] + "/{sample}_{sampling}_satellite_finder_{model}/best_solution_summary.tsv"),
	output:
		summary_positive=(dirs_dict["VIRAL_DIR"] + "/{sample}_{sampling}_satellite_finder_{model}/positive_satellites.tsv"),
		fasta_positive=(dirs_dict["VIRAL_DIR"] + "/{sample}_{sampling}_satellite_finder_{model}/positive_satellites.fasta"),
	params:
		model="{model}",
		satellite_finder_dir=os.path.join(workflow.basedir, "tools"),
	message:
		"Identifying viral satellites with satellite_finder"
	conda:
		dirs_dict["ENVS_DIR"] + "/satellite_finder.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/satellite_finder_get_fasta/model={model}__sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		grep -v "#" {input.summary} | grep -v "replicon" | awk '$2>0' | cut -f1 > {output.summary_positive}
		sed -i "s/-/_/g" {output.summary_positive} 
		seqtk subseq {input.scaffolds_spades} {output.summary_positive} > {output.fasta_positive}
		"""

rule combine_satellite_finder:
	input:
		satellite_finder_positive=expand(dirs_dict["VIRAL_DIR"] + "/{sample}_{{sampling}}_satellite_finder_{{model}}/positive_satellites.fasta" ,sample=SAMPLES),
	output:
		satellite_finder_all=(dirs_dict["VIRAL_DIR"] + "/satellite_finder_{sampling}_{model}_positive.fasta"),
	message:
		"Identifying viral satellites with satellite_finder"
	conda:
		dirs_dict["ENVS_DIR"] + "/satellite_finder.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/combine_satellite_finder/model={model}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		cat {input.satellite_finder_positive} > {output.satellite_finder_all}
		"""
		
rule genomad_viral_id:
	input:
		scaffolds_spades=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
		genomad_db=(config['genomad_db']),
	output:
		genomad_outdir=directory(dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_{sampling}/"),
		positive_contigs=dirs_dict["VIRAL_DIR"]+ "/{sample}_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		viral_fasta=lambda wc, input: dirs_dict["VIRAL_DIR"] + "/" + wc.sample + "_geNomad_" + wc.sampling + "/" + os.path.basename(input.scaffolds_spades).removesuffix(".fasta") + "_summary/" + os.path.basename(input.scaffolds_spades).removesuffix(".fasta") + "_virus.fna",
	message:
		"Identifying viral contigs with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/genomad_viral_id/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		if [ -s {input.scaffolds_spades:q} ]; then
			genomad end-to-end --restart --cleanup --splits 8 -t {threads} {input.scaffolds_spades:q} {output.genomad_outdir:q} {input.genomad_db:q} --relaxed
			cat {params.viral_fasta:q} | sed "s/|/_/g" > {output.positive_contigs:q}
		else
			mkdir -p "$(dirname {params.viral_fasta:q})"
			printf 'seq_name\\tlength\\ttopology\\tvirus_score\\ttaxonomy\\n' > "$(dirname {params.viral_fasta:q})/$(basename {input.scaffolds_spades:q} .fasta)_virus_summary.tsv"
			printf '' > {output.positive_contigs:q}
		fi
		"""

rule report_assembly_circularity:
	input:
		assembly=lambda wc: (
			dirs_dict["HOST_DIR"] + f"/{wc.sample}.fasta" if wc.assembler == "host" else dirs_dict["ASSEMBLY_DIR"] + f"/{wc.sample}_spades_filtered_scaffolds.tot.fasta"
			if wc.assembler == "spades" else RNA_DIR + f"/{wc.sample}_{wc.assembler}.fasta"
		),
		genomad=lambda wc: (
			dirs_dict["HOST_DIR"] + f"/{wc.sample}_geNomad" if wc.assembler == "host" else dirs_dict["VIRAL_DIR"] + f"/{wc.sample}_geNomad_tot/"
			if wc.assembler == "spades" else RNA_DIR + f"/{wc.sample}/{wc.assembler}_genomad"
		),
	output:
		tsv=dirs_dict["VIRAL_DIR"] + "/{sample}_{assembler}_circularity.tot.tsv",
	params:
		min_repeat=int(config.get("circularity_min_repeat_bp", 30)),
		max_repeat=int(config.get("circularity_max_repeat_bp", 5000)),
	message:
		"Reporting terminal-repeat evidence for every assembled contig"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/report_assembly_circularity/sample={sample}__assembler={assembler}.tsv"
	threads: 1
	wildcard_constraints:
		assembler="spades|rnaviralspades|megahit|trinity|host"
	run:
		import csv

		def fasta_records(path):
			name, chunks = None, []
			with open(path) as handle:
				for line in handle:
					if line.startswith(">"):
						if name is not None:
							yield name, "".join(chunks).upper().replace("U", "T")
						name, chunks = line[1:].split()[0], []
					else:
						chunks.append(line.strip())
			if name is not None:
				yield name, "".join(chunks).upper().replace("U", "T")

		prefix = os.path.splitext(os.path.basename(input.assembly))[0]
		summary = os.path.join(input.genomad, prefix + "_summary", prefix + "_virus_summary.tsv")
		viruses = {}
		with open(summary) as handle:
			for row in csv.DictReader(handle, delimiter="\t"):
				name = row["seq_name"].split("|provirus_", 1)[0]
				viruses.setdefault(name, []).append(row)

		complement = str.maketrans("ACGTN", "TGCAN")
		os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
		with open(output.tsv, "w") as handle:
			writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
			writer.writerow([
				"sample", "assembler", "contig_id", "contig_length_bp", "dtr_bp", "itr_bp",
				"dtr_sequence", "dtr_left_start", "dtr_left_end", "dtr_right_start", "dtr_right_end",
				"itr_left_sequence", "itr_right_sequence", "itr_left_start", "itr_left_end",
				"itr_right_start", "itr_right_end",
				"terminal_repeat_type", "repeat_warning", "genomad_virus_call",
				"genomad_topology", "genomad_virus_score", "genomad_taxonomy"
			])
			for name, seq in fasta_records(input.assembly):
				hits = viruses.get(name, [])
				window = min(params.max_repeat, len(seq) // 2)
				if any("DTR" in hit.get("topology", "").upper() for hit in hits):
					window = len(seq) // 2
				dtr, itr = 0, 0
				if window:
					start, end = seq[:window], seq[-window:]
					text = start + "#" + end
					border = [0] * len(text)
					for i in range(1, len(text)):
						j = border[i - 1]
						while j and text[i] != text[j]:
							j = border[j - 1]
						if text[i] == text[j]:
							j += 1
						border[i] = j
					dtr = border[-1]
					reverse_end = end.translate(complement)[::-1]
					while itr < window and start[itr] == reverse_end[itr]:
						itr += 1
				repeat_types = []
				if dtr >= params.min_repeat:
					repeat_types.append("DTR")
				if itr >= params.min_repeat:
					repeat_types.append("ITR")
				repeats = [seq[:n] for n in (dtr, itr) if n >= params.min_repeat]
				warnings = []
				if any(set(repeat) - set("ACGT") for repeat in repeats):
					warnings.append("ambiguous_bases")
				if any(len(set(repeat)) == 1 for repeat in repeats):
					warnings.append("homopolymer")
				def coordinates(length):
					return (1, length, len(seq) - length + 1, len(seq)) if length else ("", "", "", "")

				dtr_left_start, dtr_left_end, dtr_right_start, dtr_right_end = coordinates(dtr)
				itr_left_start, itr_left_end, itr_right_start, itr_right_end = coordinates(itr)
				writer.writerow([
					wildcards.sample, wildcards.assembler, name, len(seq), dtr, itr,
					seq[:dtr], dtr_left_start, dtr_left_end, dtr_right_start, dtr_right_end,
					seq[:itr], seq[-itr:] if itr else "", itr_left_start, itr_left_end,
					itr_right_start, itr_right_end,
					"+".join(repeat_types) if repeat_types else "none",
					";".join(warnings), "yes" if hits else "no",
					" | ".join(hit.get("topology", "") for hit in hits),
					" | ".join(hit.get("virus_score", "") for hit in hits),
					" | ".join(hit.get("taxonomy", "") for hit in hits)
				])

def input_genomad_viral_id_long(wildcards):
	if NANOPORE and NANOPORE_ONLY:
		return dirs_dict["ASSEMBLY_DIR"] + "/medaka_polished_{sample}_contigs_" + LONG_ASSEMBLER + ".{sampling}.fasta"
	if NANOPORE and (not NANOPORE_ONLY) and PAIRED:
		return dirs_dict["ASSEMBLY_DIR"] + "/{sample}_" + LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}.fasta"
	if PACBIO and PACBIO_ONLY:
		return dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_" + LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta"
	if PACBIO and PACBIO_HYBRID:
		return dirs_dict["ASSEMBLY_DIR"] + "/polypolish_{sample}_contigs_" + LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta"
	raise ValueError("No supported long-read mode detected for geNomad")


rule genomad_viral_id_long:
	input:
		scaffolds=input_genomad_viral_id_long,
		genomad_db=config["genomad_db"]
	output:
		genomad_outdir=directory(dirs_dict["VIRAL_DIR"] + "/{sample}_long_geNomad_{sampling}/"),
		positive_contigs=dirs_dict["VIRAL_DIR"] + "/{sample}_long_" + VIRAL_CONTIGS_BASE + ".{sampling}.fasta",
	params:
		viral_fasta=dirs_dict["VIRAL_DIR"] + "/{sample}_long_geNomad_{sampling}/*summary/*_virus.fna",
		genomad_filter="--" + config["genomad_filter"],
	message:
		"Identifying viral contigs with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/genomad_viral_id_long/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.scaffolds} {output.genomad_outdir} {input.genomad_db} {params.genomad_filter}
		cat {params.viral_fasta} 2>/dev/null | sed "s/|/_/g" > {output.positive_contigs}
		"""

# VIRAL FILTERING vOTUS
rule virSorter2:
	input:
		representatives=lambda wc: annotation_fasta_path(wc.sequence + "." + wc.sampling),
		virSorter_db=config['virSorter_db'],
	output:
		positive_fasta=dirs_dict["vOUT_DIR"] + "/VirSorter2_{sequence}_{sampling}/final-viral-combined.fa",
		table_virsorter=dirs_dict["vOUT_DIR"] + "/VirSorter2_{sequence}_{sampling}/final-viral-score.tsv",
		positive_list=dirs_dict["vOUT_DIR"] + "/VirSorter2_{sequence}_{sampling}/positive_VS_list_{sampling}.txt",
		# DRAM_tab=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/for-dramv/viral-affi-contigs-for-dramv.tab",
		# DRAM_fasta=dirs_dict["vOUT_DIR"] + "/VirSorter2_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/for-dramv/final-viral-combined-for-dramv.fa",
		iter=directory(dirs_dict["vOUT_DIR"] + "/VirSorter2_{sequence}_{sampling}/iter-0"),
	params:
		out_folder=dirs_dict["vOUT_DIR"] + "/VirSorter2_{sequence}_{sampling}"
	message:
		"Classifing contigs with VirSorter"
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/virSorter2/sequence={sequence}__sampling={sampling}.tsv"
	threads: 64
	wildcard_constraints:
		sequence="[^/]+",
		sampling="tot|sub"
	shell:
		"""
		if [ -s {input.representatives:q} ]; then
		virsorter run -w {params.out_folder:q} -i {input.representatives:q} -j {threads} --db-dir {input.virSorter_db:q} \
				--include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --seqname-suffix-off  --provirus-off --min-length 0
		else
			mkdir -p {output.iter:q}
			: > {output.positive_fasta:q}
			printf 'seqname\tmax_score\tmax_score_group\n' > {output.table_virsorter:q}
		fi
		grep ">" {output.positive_fasta:q} | cut -f1 -d\| | sed "s/>//g" > {output.positive_list:q} || true
		"""

rule genomad_annotations:
	input:
		fasta=lambda wc: annotation_fasta_path(wc.sequence),
		database=config["genomad_db"],
	output:
		results=directory(dirs_dict["ANNOTATION"] + "/geNomad_{sequence}"),
	params:
		prefix="{sequence}",
	message:
		"Collecting geNomad virus/plasmid evidence for {wildcards.sequence}"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/genomad_annotations/sequence={sequence}.tsv"
	threads: 16
	resources:
		mem_mb=16000,
	wildcard_constraints:
		sequence="[^/]+",
	shell:
		r"""
		if [ -s {input.fasta:q} ]; then
			genomad end-to-end --relaxed --cleanup --splits 8 -t {threads} \
				{input.fasta:q} {output.results:q} {input.database:q}
		else
			mkdir -p {output.results:q}/{params.prefix:q}_summary
			printf 'seq_name\tvirus_score\ttaxonomy\ttopology\n' > {output.results:q}/{params.prefix:q}_summary/{params.prefix:q}_virus_summary.tsv
			printf 'seq_name\tplasmid_score\ttopology\n' > {output.results:q}/{params.prefix:q}_summary/{params.prefix:q}_plasmid_summary.tsv
		fi
		"""

rule genomad_vOTUs:
	input:
		representatives=dirs_dict["vOUT_DIR"]+ "/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		genomad_db=(config['genomad_db']),
	output:
		virus_summary=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_virus_summary.tsv",
		plasmid_summary=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_plasmid_summary.tsv",
		viral_fasta=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_virus.fna",																																																
		positive_contigs=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/formatted_viral_" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}.fasta",
		positive_contigs_conservative=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/" + REPRESENTATIVE_CONTIGS_BASE + ".{sampling}_summary/formatted_viral_" + REPRESENTATIVE_CONTIGS_BASE + "_conservative.{sampling}.fasta",
	params:
		genomad_outdir=dirs_dict["vOUT_DIR"] + "/geNomad_" + REPRESENTATIVE_CONTIGS_BASE + "_{sampling}/",
	message:
		"Identifying viral contigs with geNomad"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/genomad_vOTUs/sampling={sampling}.tsv"
	threads: 32
	shell:
		"""
		rm -rf {params.genomad_outdir}
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.representatives} {params.genomad_outdir} {input.genomad_db} --conservative  
		cat {output.viral_fasta} | sed "s/|/_/g" > {output.positive_contigs_conservative}
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.representatives} {params.genomad_outdir} {input.genomad_db}  
		cat {output.viral_fasta} | sed "s/|/_/g" > {output.positive_contigs}
		# leave genomad folder as --relaxed to have all contigs in the summary file
		genomad end-to-end --cleanup --splits 8 -t {threads} {input.representatives} {params.genomad_outdir} {input.genomad_db}  --relaxed
		"""

rule annotate_VIBRANT:
	input:
		representatives=lambda wc: annotation_fasta_path(wc.sequence + "." + wc.sampling),
		VIBRANT_dir=os.path.join(workflow.basedir, config['vibrant_dir']),
	output:
		vibrant_circular=dirs_dict["vOUT_DIR"] + "/VIBRANT_{sequence}_circular.{sampling}.csv",
		vibrant_positive=dirs_dict["vOUT_DIR"] + "/VIBRANT_{sequence}_positive_list.{sampling}.csv",
		vibrant_quality=dirs_dict["vOUT_DIR"] + "/VIBRANT_{sequence}_positive_quality.{sampling}.csv",
		vibrant_summary=dirs_dict["vOUT_DIR"] + "/VIBRANT_{sequence}_summary_results.{sampling}.csv",
	params:
		vibrant_outdir=dirs_dict["vOUT_DIR"] + "/VIBRANT_{sequence}.{sampling}",
	conda:
		dirs_dict["ENVS_DIR"] + "/vibrant.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/annotate_VIBRANT/sequence={sequence}__sampling={sampling}.tsv"
	message:
		"Annotating viral contigs with VIBRANT"
	threads: 16
	wildcard_constraints:
		sequence="[^/]+",
		sampling="tot|sub"
	shell:
		"""
		rm -rf -- {params.vibrant_outdir:q}
		mkdir -p {params.vibrant_outdir:q}
		if [ -s {input.representatives:q} ]; then
			{input.VIBRANT_dir:q}/VIBRANT_run.py -i {input.representatives:q} -t {threads} -virome -folder {params.vibrant_outdir:q}
		fi
		: > {output.vibrant_circular:q}
		: > {output.vibrant_positive:q}
		printf 'scaffold\ttype\tQuality\n' > {output.vibrant_quality:q}
		printf 'scaffold\n' > {output.vibrant_summary:q}
		shopt -s nullglob
		for path in {params.vibrant_outdir:q}/VIBRANT_*/VIBRANT_results*/*complete_circular*tsv; do
			cut -f1 "$path" > {output.vibrant_circular:q}
		done
		for path in {params.vibrant_outdir:q}/VIBRANT_*/VIBRANT_phages_*/*phages_combined.txt; do
			cp "$path" {output.vibrant_positive:q}
		done
		for path in {params.vibrant_outdir:q}/VIBRANT_*/VIBRANT_results*/VIBRANT_genome_quality*.tsv; do
			cp "$path" {output.vibrant_quality:q}
		done
		for path in {params.vibrant_outdir:q}/VIBRANT_*/VIBRANT_results*/VIBRANT_summary_results*.tsv; do
			cp "$path" {output.vibrant_summary:q}
		done
		"""
