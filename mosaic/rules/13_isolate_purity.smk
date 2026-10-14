# Full-read isolate purity: selection, shared reference pools and one reusable mapping stage.
def isolate_stage_reference(wildcards):
	stage=wildcards.stage
	own={"01_own_retained": "retained", "04_own_host_excluded": "excluded_host", "05_own_threshold_excluded": "excluded_threshold"}
	if stage in own:
		return ISOLATE_CONTIG_DIR + "/" + wildcards.sample + "_" + own[stage] + ".tot.fasta"
	pool={"02_host_chromosomes": "host_chromosomes", "03_host_viral": "host_viral", "06_other_retained": "retained",
		"07_other_host_excluded": "excluded_host", "08_other_threshold_excluded": "excluded_threshold"}[stage]
	scope=("HOST/" + ISOLATE_HOST_ASSIGNMENTS[wildcards.sample] + "/" if wildcards.sample in ISOLATE_HOST_ASSIGNMENTS else "UNASSIGNED/") if stage in ["02_host_chromosomes", "03_host_viral"] else ""
	return ISOLATE_CONTIG_DIR + "/REFERENCES/" + scope + pool + ".fasta"

def isolate_stage_reads(wildcards, mate):
	index=ISOLATE_STAGES.index(wildcards.stage)
	if index == 0:
		name={"R1": "forward_paired", "R2": "reverse_paired", "U": "unpaired"}[mate]
		return dirs_dict["CLEAN_DATA_DIR"] + "/" + wildcards.sample + "_" + name + "_clean.tot.fastq.gz"
	return ISOLATE_MAPPING_DIR + "/" + wildcards.sample + "/" + ISOLATE_STAGES[index - 1] + ".remaining_" + mate + ".fastq.gz"

def isolate_stage_index(wildcards):
	if wildcards.stage == "01_own_retained":
		return []
	return [isolate_stage_reference(wildcards) + "." + part + ".bt2l" for part in ["1", "2", "3", "4", "rev.1", "rev.2"]]

def isolate_reference_hosts(wildcards):
	if wildcards.reference_set in ["retained", "excluded", "excluded_host", "excluded_threshold"] or wildcards.reference_scope == "UNASSIGNED/":
		return []
	if wildcards.reference_scope:
		return [wildcards.reference_scope.split("/")[1]]
	return HOSTS

rule select_isolate_contigs:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.tot.fasta",
		assembled_covstats=dirs_dict["MAPPING_DIR"] + "/STATS_FILES/bowtie2_{sample}_assembled_contigs.tot_covstats.txt",
		host_blast=lambda wc: [dirs_dict["vOUT_DIR"] + "/ISOLATES/" + wc.sample + "/blastn_out_assembly_" + ISOLATE_HOST_ASSIGNMENTS[wc.sample] + ".tot.csv"] if wc.sample in ISOLATE_HOST_ASSIGNMENTS else [],
		host_viral_blast=lambda wc: [dirs_dict["vOUT_DIR"] + "/ISOLATES/" + wc.sample + "/blastn_out_assembly_" + ISOLATE_HOST_ASSIGNMENTS[wc.sample] + "_viral_regions.tot.csv"] if wc.sample in ISOLATE_HOST_ASSIGNMENTS else [],
		genomad=dirs_dict["VIRAL_DIR"] + "/{sample}_geNomad_tot/",
		host_genomad=lambda wc: [dirs_dict["HOST_DIR"] + "/" + ISOLATE_HOST_ASSIGNMENTS[wc.sample] + "_geNomad"] if wc.sample in ISOLATE_HOST_ASSIGNMENTS else [],
	output:
		selection=ISOLATE_CONTIG_DIR + "/SELECTION/{sample}.tot.tsv",
	params:
		phase="selection",
		samples=lambda wc: [wc.sample],
		hosts=lambda wc: [ISOLATE_HOST_ASSIGNMENTS[wc.sample]] if wc.sample in ISOLATE_HOST_ASSIGNMENTS else [],
		expected_host=lambda wc: ISOLATE_HOST_ASSIGNMENTS.get(wc.sample, "not reported"),
		sampling="tot",
		results_dir=RESULTS_DIR,
		min_depth=float(config.get("isolate_min_depth", 5)),
		min_length=int(config.get("isolate_min_length_bp", 4000)),
		short_depth=float(config.get("isolate_short_min_depth", 10)),
		host_min_identity=float(config.get("isolate_host_min_identity", 90)),
		host_min_query_coverage=float(config.get("isolate_host_min_query_coverage", 90)),
		host_plasmid_min_identity=float(config.get("isolate_host_plasmid_min_identity", 95)),
		host_plasmid_min_query_coverage=float(config.get("isolate_host_plasmid_min_query_coverage", 75)),
		host_plasmid_min_subject_coverage=float(config.get("isolate_host_plasmid_min_subject_coverage", 90)),
	benchmark:
		dirs_dict["BENCHMARKS"] + "/select_isolate_contigs/sample={sample}.tsv"
	# Reuse selection code without saving a per-sample report notebook.
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_isolate_contig_catalogue.py.ipynb"

rule extract_isolate_contig_sets:
	input:
		metadata=ISOLATE_CONTIG_DIR + "/SELECTION/{sample}.tot.tsv",
		fasta=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.tot.fasta",
	output:
		retained=ISOLATE_CONTIG_DIR + "/{sample}_retained.tot.fasta",
		excluded=ISOLATE_CONTIG_DIR + "/{sample}_excluded.tot.fasta",
		excluded_host=ISOLATE_CONTIG_DIR + "/{sample}_excluded_host.tot.fasta",
		excluded_threshold=ISOLATE_CONTIG_DIR + "/{sample}_excluded_threshold.tot.fasta",
		groups=ISOLATE_CONTIG_DIR + "/SELECTION/{sample}.tot.exclusion_groups.tsv",
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/extract_isolate_contig_sets/sample={sample}.tsv"
	threads: 1
	shell:
		r"""
		python - {input.metadata:q} {input.fasta:q} {output.retained:q} {output.excluded:q} {output.excluded_host:q} {output.excluded_threshold:q} {output.groups:q} <<-'PYTHON'
		import sys
		import numpy as np
		import pandas as pd
		from Bio import SeqIO

		metadata_path, fasta_path, retained_path, excluded_path, host_path, threshold_path, groups_path=sys.argv[1:]
		metadata=pd.read_csv(metadata_path, sep="\t", dtype={{"best_host_blast_contig": str, "host_plasmid_contig": str}})
		retained=set(metadata.loc[metadata["contig_set"].eq("retained"), "original_id"])
		prophage=metadata.get("host_prophage_match", pd.Series(False, index=metadata.index)).isin([True, "True"])
		host=metadata.get("host_match", pd.Series(False, index=metadata.index)).isin([True, "True"]) | prophage
		excluded=metadata["contig_set"].eq("excluded")
		metadata["exclusion_group"]=np.select([excluded & host, excluded], ["excluded_host", "excluded_threshold"], default="retained")
		reasons=metadata.get("exclusion_reason", pd.Series("not reported", index=metadata.index))
		metadata["primary_exclusion_reason"]=np.select([excluded & prophage, excluded & host,
		    excluded & reasons.str.contains("below_minimum_depth", regex=False, na=False),
		    excluded & reasons.str.contains("short_contig_below_required_depth", regex=False, na=False)],
		    ["host_prophage_associated", "host_chromosome_match", "below_minimum_depth", "short_contig_below_required_depth"],
		    default="not excluded")
		metadata.loc[metadata["exclusion_group"].eq("excluded_threshold") & metadata["primary_exclusion_reason"].eq("not excluded"),
		    "primary_exclusion_reason"]="length_depth_threshold"
		host_ids=set(metadata.loc[metadata["exclusion_group"].eq("excluded_host"), "original_id"])
		metadata.to_csv(groups_path, sep="\t", index=False, na_rep="not reported")
		with open(retained_path, "w") as kept, open(excluded_path, "w") as removed, open(host_path, "w") as host_removed, open(threshold_path, "w") as threshold_removed:
		    for record in SeqIO.parse(fasta_path, "fasta"):
		        if record.id in retained:
		            SeqIO.write(record, kept, "fasta")
		        else:
		            SeqIO.write(record, removed, "fasta")
		            SeqIO.write(record, host_removed if record.id in host_ids else threshold_removed, "fasta")
		PYTHON
		"""

rule pool_isolate_references:
	input:
		contigs=lambda wc: expand(ISOLATE_CONTIG_DIR + "/{sample}_" + wc.reference_set + ".tot.fasta", sample=SAMPLES) if wc.reference_set in ["retained", "excluded", "excluded_host", "excluded_threshold"] else [],
		fasta=lambda wc: [ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_rep_seq.tot.fasta"] if wc.reference_set in ["retained", "excluded", "excluded_host", "excluded_threshold"] else [],
		exact=lambda wc: [ALL_ASSEMBLED_DIR + "/phage_isolates_contigs_derreplicated_cluster.tot.tsv"] if wc.reference_set in ["retained", "excluded", "excluded_host", "excluded_threshold"] else [],
		hosts=lambda wc: expand(dirs_dict["HOST_DIR"] + "/host_masked_prophages/{host}_masked_prophages.fasta", host=isolate_reference_hosts(wc)) if wc.reference_set == "host_chromosomes" else [],
		viral=lambda wc: expand(dirs_dict["HOST_DIR"] + "/prophages/{host}_prophages.fasta", host=isolate_reference_hosts(wc)) if wc.reference_set == "host_viral" else [],
	output:
		fasta=ISOLATE_CONTIG_DIR + "/REFERENCES/{reference_scope}{reference_set}.fasta",
		membership=ISOLATE_CONTIG_DIR + "/REFERENCES/{reference_scope}{reference_set}.membership.tsv",
	wildcard_constraints:
		reference_set="host_chromosomes|host_viral|retained|excluded|excluded_host|excluded_threshold",
		reference_scope="|HOST/[^/]+/|UNASSIGNED/",
	params:
		context=lambda wc, input: json.dumps(dict(samples=SAMPLES, hosts=isolate_reference_hosts(wc),
			contigs=list(input.contigs), fasta=list(input.fasta), exact=list(input.exact),
			host_fastas=list(input.hosts), viral_fastas=list(input.viral), sample_hosts=ISOLATE_HOST_ASSIGNMENTS)),
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/pool_isolate_references/{reference_scope}reference={reference_set}.tsv"
	threads: 1
	shell:
		r"""
		python - {wildcards.reference_set:q} {output.fasta:q} {output.membership:q} {params.context:q} <<-'PYTHON'
		import csv
		import json
		import sys
		import pandas as pd
		from Bio import SeqIO

		reference_set, fasta_path, membership_path, context=sys.argv[1:]
		context=json.loads(context)
		fields=["reference_id", "sample", "original_id", "host", "evidence_scope", "source_expected_host"]
		with open(fasta_path, "w") as fasta, open(membership_path, "w") as table:
		    writer=csv.DictWriter(table, fieldnames=fields, delimiter="\t", lineterminator="\n")
		    writer.writeheader()
		    if reference_set in ["retained", "excluded", "excluded_host", "excluded_threshold"]:
		        members=[]
		        for sample, path in zip(context["samples"], context["contigs"]):
		            with open(path) as handle:
		                members.extend({{"sample": sample, "original_id": record.id}} for record in SeqIO.parse(handle, "fasta"))
		        members=pd.DataFrame(members, columns=["sample", "original_id"])
		        exact=pd.read_csv(context["exact"][0], sep="\t", names=["exact_rep", "original_id"])
		        members=members.merge(exact, on="original_id", how="left")
		        groups={{rep: group for rep, group in members.groupby("exact_rep")}}
		        with open(context["fasta"][0]) as handle:
		            for record in SeqIO.parse(handle, "fasta"):
		                if record.id not in groups:
		                    continue
		                SeqIO.write(record, fasta, "fasta")
		                for row in groups[record.id].itertuples():
		                    writer.writerow(dict(reference_id=record.id, sample=row.sample,
		                        original_id=row.original_id, evidence_scope=reference_set,
		                        source_expected_host=context.get("sample_hosts", {{}}).get(row.sample, "not reported")))
		    else:
		        for host, path in zip(context["hosts"], context["host_fastas"] or context["viral_fastas"]):
		            with open(path) as handle:
		                for record in SeqIO.parse(handle, "fasta"):
		                    original=record.id
		                    if reference_set == "host_chromosomes":
		                        record.id=host + "_" + original
		                    record.description=""
		                    SeqIO.write(record, fasta, "fasta")
		                    writer.writerow(dict(reference_id=record.id, original_id=original, host=host,
		                        evidence_scope=reference_set))
		PYTHON
		"""

rule buildBowtieDB_isolate_stage:
	input:
		fasta="{reference}.fasta",
	output:
		index=temp(expand("{{reference}}.fasta.{part}.bt2l", part=["1", "2", "3", "4", "rev.1", "rev.2"])),
	wildcard_constraints:
		reference=re.escape(ISOLATE_CONTIG_DIR) + "/(?:(?:REFERENCES/(?:HOST/[^/]+/|UNASSIGNED/)?(?:host_chromosomes|host_viral|retained|excluded|excluded_host|excluded_threshold))|(?:[^/]+_excluded(?:_host|_threshold)?.tot))",
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_isolate_stage/reference={reference}.tsv"
	threads: int(config.get("isolate_mapping_threads", 8))
	shell:
		r"""
		if [ -s {input.fasta:q} ]; then
			bowtie2-build --large-index --threads {threads} {input.fasta:q} {input.fasta:q}
		else
			touch {output.index:q}
		fi
		"""

rule isolate_read_accounting:
	input:
		fasta=isolate_stage_reference,
		index=isolate_stage_index,
		forward=lambda wc: isolate_stage_reads(wc, "R1"),
		rev_reads=lambda wc: isolate_stage_reads(wc, "R2"),
		unpaired=lambda wc: isolate_stage_reads(wc, "U"),
		assembly_bam=lambda wc: [dirs_dict["MAPPING_DIR"] + "/STATS_FILES/bowtie2_" + wc.sample + "_assembled_contigs_tot_filtered.bam"] if wc.stage == "01_own_retained" else [],
		membership=lambda wc: [isolate_stage_reference(wc).removesuffix(".fasta") + ".membership.tsv"] if wc.stage in ["02_host_chromosomes", "03_host_viral", "06_other_retained", "07_other_host_excluded", "08_other_threshold_excluded"] else [],
		previous=lambda wc: [ISOLATE_MAPPING_DIR + "/" + wc.sample + "/" + ISOLATE_STAGES[ISOLATE_STAGES.index(wc.stage)-1] + ".summary.tsv"] if wc.stage != "01_own_retained" else [],
	output:
		bam=ISOLATE_MAPPING_DIR + "/{sample}/{stage}.bam",
		bam_index=ISOLATE_MAPPING_DIR + "/{sample}/{stage}.bam.bai",
		reads=ISOLATE_MAPPING_DIR + "/{sample}/{stage}.reads.tsv.gz",
		summary=ISOLATE_MAPPING_DIR + "/{sample}/{stage}.summary.tsv",
		reference_counts=ISOLATE_MAPPING_DIR + "/{sample}/{stage}.reference_reads.tsv",
		covstats=ISOLATE_MAPPING_DIR + "/{sample}/{stage}.covstats.tsv",
		basecov=(ISOLATE_MAPPING_DIR + "/{sample}/{stage}.basecov.tsv.gz") if PER_BASE_COVERAGE else [],
	wildcard_constraints:
		stage="|".join(ISOLATE_STAGES),
		sample="|".join(re.escape(sample) for sample in SAMPLES) or "(?!)",
	params:
		expected_host=lambda wc: ISOLATE_HOST_ASSIGNMENTS.get(wc.sample, "not reported"),
		assembly_bam=lambda wc, input: shlex.quote(input.assembly_bam[0] if input.assembly_bam else ""),
		membership=lambda wc, input: shlex.quote(input.membership[0] if input.membership else ""),
		previous=lambda wc, input: shlex.quote(input.previous[0] if input.previous else ""),
		basecov=lambda wc, output: shlex.quote(str(output.basecov) if output.basecov else ""),
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/isolate_read_accounting/sample={sample}__stage={stage}.tsv"
	log:
		ISOLATE_MAPPING_DIR + "/{sample}/{stage}.log"
	threads: int(config.get("isolate_mapping_threads", 8))
	resources:
		mem_mb=int(config.get("isolate_mapping_mem_mb", 16000)),
	shell:
		r"""
		python - {wildcards.sample:q} {wildcards.stage:q} {params.expected_host:q} {input.fasta:q} {input.forward:q} {input.rev_reads:q} {input.unpaired:q} \
			{params.assembly_bam} {params.membership} {params.previous} {output.bam:q} {output.reads:q} \
			{output.summary:q} {output.covstats:q} {params.basecov} {log:q} {threads} <<-'PYTHON'
		import csv
		import gzip
		import os
		import re
		import subprocess
		import sys
		import tempfile
		from collections import Counter

		(sample, stage, expected_host, fasta, r1, r2, orphan, assembly_bam, membership, previous,
		 bam, reads_path, summary_path, covstats, basecov, log_path, threads) = sys.argv[1:]
		threads = int(threads)

		def fastq(path):
		    with gzip.open(path, "rt") as handle:
		        for header in handle:
		            yield header, handle.readline(), handle.readline(), handle.readline()

		def read_name(header):
		    name = header.split()[0].lstrip("@")
		    return re.sub(r"/[12]$", "", name)

		def original_identity(name, mate):
		    match = re.fullmatch(r"(.+)__MOSAIC_mate([12])", name)
		    return (match.group(1), int(match.group(2))) if match else (name, mate)

		lengths = {{}}
		with open(fasta) as handle:
		    name = None
		    for line in handle:
		        if line.startswith(">"):
		            name = line[1:].split()[0]
		            lengths[name] = 0
		        elif name is not None:
		            lengths[name] += len(line.strip())

		owners = {{}}
		if membership:
		    with open(membership) as handle:
		        for row in csv.DictReader(handle, delimiter="\t"):
		            owners.setdefault(row["reference_id"], []).append(row)
		allowed = set(lengths)
		host_scope = "not_applicable"
		if stage in ["02_host_chromosomes", "03_host_viral"]:
		    host_scope = "expected_host" if expected_host != "not reported" else "not_assessed"
		    allowed = {{name for name in allowed if any(row["host"] == expected_host for row in owners.get(name, []))}}
		if stage in ["06_other_retained", "07_other_host_excluded", "08_other_threshold_excluded"]:
		    # A shared exact representative can belong to this sample AND another sample.
		    allowed = {{name for name in allowed if any(row["sample"] != sample for row in owners.get(name, []))}}

		def nonempty(path):
		    with gzip.open(path, "rt") as handle:
		        return bool(handle.read(1))

		accepted = {{}}
		with tempfile.TemporaryDirectory(prefix=".accounting_", dir=os.path.dirname(bam)) as workspace, open(log_path, "w") as log:
		    alignment = assembly_bam
		    if not alignment and allowed and (nonempty(r1) or nonempty(orphan)):
		        sam = os.path.join(workspace, "bowtie2.sam")
		        command = ["bowtie2", "-x", fasta, "--threads", str(threads), "--very-sensitive", "--all", "-S", sam]
		        if nonempty(r1):
		            command += ["-1", r1, "-2", r2]
		        if nonempty(orphan):
		            command += ["-U", orphan]
		        subprocess.run(command, check=True, stderr=log)
		        sorted_bam = os.path.join(workspace, "sorted.bam")
		        subprocess.run(["samtools", "sort", "-@", str(threads), "-o", sorted_bam, sam], check=True, stderr=log)
		        alignment = os.path.join(workspace, "accepted.bam")
		        subprocess.run(["coverm", "filter", "-b", sorted_bam, "-o", alignment,
		                        "--min-read-percent-identity", "95", "--min-read-aligned-percent", "85",
		                        "--include-secondary", "--exclude-supplementary",
		                        "-t", str(threads)], check=True, stderr=log)
		    headers = []
		    if alignment:
		        process = subprocess.Popen(["samtools", "view", "-h", alignment], stdout=subprocess.PIPE, text=True, stderr=log)
		        for line in process.stdout:
		            if line.startswith("@"):
		                headers.append(line)
		                continue
		            fields = line.rstrip("\n").split("\t")
		            flag = int(fields[1])
		            if flag & 4 or fields[2] not in allowed:
		                continue
		            mate = 1 if flag & 64 else 2 if flag & 128 else 0
		            key = (read_name(fields[0]), mate)
		            tags = {{tag.split(":", 2)[0]: tag.split(":", 2)[-1] for tag in fields[11:]}}
		            score = (not bool(flag & (256 | 2048)), int(tags.get("AS", 0)), -int(tags.get("NM", 0)),
		                     int(fields[4]), fields[2], -int(fields[3]))
		            evidence = accepted.setdefault(key, {{"targets": set(), "best": None, "score": None, "xs": False}})
		            evidence["targets"].add(fields[2])
		            evidence["xs"] = evidence["xs"] or "XS" in tags
		            if evidence["score"] is None or score > evidence["score"]:
		                evidence["best"], evidence["score"] = fields, score
		        process.stdout.close()
		        if process.wait():
		            sys.exit("samtools view failed; see " + log_path)
		    else:
		        headers = ["@HD\tVN:1.6\tSO:unsorted\n"] + [f"@SQ\tSN:{{name}}\tLN:{{length}}\n" for name, length in lengths.items()]
		        log.write("No reference sequence or no entering reads; mapping skipped.\n")

		    # One deterministic accepted placement per READ contributes to coverage.
		    # Alternative accepted references remain in the read-level evidence table.
		    selected_sam = os.path.join(workspace, "selected.sam")
		    with open(selected_sam, "w") as output:
		        output.writelines(headers)
		        for (name, mate), evidence in accepted.items():
		            fields = list(evidence["best"])
		            flag = int(fields[1]) & ~(256 | 2048)
		            fields[0] = name
		            other = accepted.get((name, 3 - mate)) if mate else None
		            if mate:
		                flag &= ~(2 | 8 | 32)
		                if other:
		                    partner = other["best"]
		                    proper = (int(evidence["best"][1]) & 2 and int(partner[1]) & 2
		                              and fields[2] == partner[2]
		                              and int(evidence["best"][7]) == int(partner[3])
		                              and int(partner[7]) == int(fields[3]))
		                    if proper:
		                        flag |= 2
		                    if int(partner[1]) & 16:
		                        flag |= 32
		                    fields[6] = "=" if fields[2] == partner[2] else partner[2]
		                    fields[7] = partner[3]
		                    if not proper:
		                        fields[8] = "0"
		                else:
		                    flag |= 8
		                    fields[6:9] = ["*", "0", "0"]
		            fields[1] = str(flag)
		            output.write("\t".join(fields) + "\n")
		    subprocess.run(["samtools", "sort", "-@", str(threads), "-o", bam, selected_sam], check=True, stderr=log)
		    subprocess.run(["samtools", "index", bam], check=True, stderr=log)
		    if accepted:
		        subprocess.run(["coverm", "contig", "-b", bam, "-m", "mean", "length", "covered_bases",
		                        "count", "variance", "trimmed_mean", "rpkm", "-t", str(threads),
		                        "-o", covstats], check=True, stderr=log)
		        # Initial own-assembly headers can also contain excluded references.
		        with open(covstats) as handle:
		            rows = list(csv.reader(handle, delimiter="\t"))
		        with open(covstats, "w") as handle:
		            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
		            writer.writerow(rows[0])
		            writer.writerows(row for row in rows[1:] if row[0] in lengths)
		        if basecov:
		            with gzip.open(basecov, "wt") as output:
		                process = subprocess.Popen(["bedtools", "genomecov", "-dz", "-ibam", bam], stdout=subprocess.PIPE, text=True, stderr=log)
		                for line in process.stdout:
		                    output.write(line)
		                process.stdout.close()
		                if process.wait():
		                    sys.exit("bedtools genomecov failed; see " + log_path)
		    else:
		        with open(covstats, "w") as output:
		            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
		            writer.writerow(["Contig", "Mean", "Length", "Covered Bases", "Read Count", "Variance", "Trimmed Mean", "RPKM"])
		            writer.writerows([name, 0, length, 0, 0, 0, 0, 0] for name, length in lengths.items())
		        if basecov:
		            with gzip.open(basecov, "wt") as output:
		                pass

		    totals = Counter()
		    reference_counts = Counter()
		    columns = ["sample", "stage", "read_name", "original_read_name", "original_mate", "alignment_mate",
		               "selected_reference", "accepted_references", "number_accepted_references",
		               "non_unique", "pairing_scope", "matching_samples", "matching_hosts"]
		    with gzip.open(reads_path, "wt") as table:
		        writer = csv.DictWriter(table, fieldnames=columns, delimiter="\t", lineterminator="\n")
		        writer.writeheader()

		        def account(record, mate):
		            name = read_name(record[0])
		            evidence = accepted.get((name, mate))
		            totals["entering_reads"] += 1
		            if evidence is None:
		                totals["remaining_reads"] += 1
		                return False
		            totals["assigned_reads"] += 1
		            reference_counts[evidence["best"][2]] += 1
		            targets = sorted(evidence["targets"])
		            non_unique = len(targets) > 1 or evidence["xs"]
		            totals["non_unique_reads"] += non_unique
		            other = accepted.get((name, 3 - mate)) if mate else None
		            if mate:
		                own_fields = evidence["best"]
		                proper = (other is not None and int(own_fields[1]) & 2 and int(other["best"][1]) & 2
		                          and own_fields[2] == other["best"][2]
		                          and int(own_fields[7]) == int(other["best"][3])
		                          and int(other["best"][7]) == int(own_fields[3]))
		                pairing = "proper_pair" if proper else "discordant" if other else "single_mate"
		            else:
		                pairing = "residual_mate" if "__MOSAIC_mate" in name else "unpaired"
		            totals[pairing + "_reads"] += 1
		            original_name, original_mate = original_identity(name, mate)
		            matches = [row for target in targets for row in owners.get(target, [])]
		            writer.writerow(dict(sample=sample, stage=stage, read_name=name, original_read_name=original_name,
		                original_mate=original_mate, alignment_mate=mate, selected_reference=evidence["best"][2],
		                accepted_references=";".join(targets), number_accepted_references=len(targets),
		                non_unique=non_unique, pairing_scope=pairing,
		                matching_samples=";".join(sorted({{row["sample"] for row in matches if row["sample"]}})),
		                matching_hosts=";".join(sorted({{row["host"] for row in matches if row["host"]}}))))
		            return True

		        for first, second in zip(fastq(r1), fastq(r2)):
		            account(first, 1)
		            account(second, 2)
		        for record in fastq(orphan):
		            account(record, 0)

		    original_total = totals["entering_reads"]
		    if previous:
		        with open(previous) as handle:
		            original_total = int(next(csv.DictReader(handle, delimiter="\t"))["original_reads"])
		    summary = dict(sample=sample, stage=stage, reference_assessment="assessed" if allowed else "no_eligible_reference",
		                   expected_host=expected_host, host_reference_scope=host_scope, original_reads=original_total)
		    summary.update({{name: totals[name] for name in ["entering_reads", "assigned_reads", "remaining_reads", "non_unique_reads",
		                    "proper_pair_reads", "discordant_reads", "single_mate_reads", "residual_mate_reads", "unpaired_reads"]}})
		    summary["percent_original"] = 100 * totals["assigned_reads"] / original_total if original_total else 0
		    summary["percent_entering"] = 100 * totals["assigned_reads"] / totals["entering_reads"] if totals["entering_reads"] else 0
		    summary["remaining_percent_original"] = 100 * totals["remaining_reads"] / original_total if original_total else 0
		    with open(summary_path, "w") as output:
		        writer = csv.DictWriter(output, fieldnames=list(summary), delimiter="\t", lineterminator="\n")
		        writer.writeheader()
		        writer.writerow(summary)
		    with open(summary_path.replace(".summary.tsv", ".reference_reads.tsv"), "w") as output:
		        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
		        writer.writerow(["sample", "stage", "reference", "assigned_reads", "genome_length"])
		        writer.writerows([sample, stage, name, reference_counts[name], length] for name, length in lengths.items())
		PYTHON
		"""

rule extract_isolate_remaining_reads:
	input:
		forward=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.tot.fastq.gz",
		rev_reads=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.tot.fastq.gz",
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.tot.fastq.gz",
		assigned=lambda wc: [ISOLATE_MAPPING_DIR + "/" + wc.sample + "/" + stage + ".reads.tsv.gz" for stage in ISOLATE_STAGES[:ISOLATE_STAGES.index(wc.stage) + 1]],
	output:
		forward=temp(ISOLATE_MAPPING_DIR + "/{sample}/{stage}.remaining_R1.fastq.gz"),
		rev_reads=temp(ISOLATE_MAPPING_DIR + "/{sample}/{stage}.remaining_R2.fastq.gz"),
		unpaired=temp(ISOLATE_MAPPING_DIR + "/{sample}/{stage}.remaining_U.fastq.gz"),
	wildcard_constraints:
		stage="|".join(ISOLATE_STAGES),
		sample="|".join(re.escape(sample) for sample in SAMPLES) or "(?!)",
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/extract_isolate_remaining_reads/sample={sample}__stage={stage}.tsv"
	threads: 1
	shell:
		r"""
		python - {input.forward:q} {input.rev_reads:q} {input.unpaired:q} \
			{output.forward:q} {output.rev_reads:q} {output.unpaired:q} {input.assigned:q} <<-'PYTHON'
		import csv
		import gzip
		import re
		import sys

		r1, r2, orphan, remaining_r1, remaining_r2, remaining_u=sys.argv[1:7]
		assigned=set()
		for path in sys.argv[7:]:
		    with gzip.open(path, "rt") as handle:
		        for row in csv.DictReader(handle, delimiter="\t"):
		            assigned.add((row["original_read_name"], int(row["original_mate"])))

		def fastq(path):
		    with gzip.open(path, "rt") as handle:
		        for header in handle:
		            yield header, handle.readline(), handle.readline(), handle.readline()

		def read_name(header):
		    return re.sub(r"/[12]$", "", header.split()[0].lstrip("@"))

		with gzip.open(remaining_r1, "wt") as forward, gzip.open(remaining_r2, "wt") as rev_reads, gzip.open(remaining_u, "wt") as unpaired:
		    for first, second in zip(fastq(r1), fastq(r2)):
		        remaining_first=(read_name(first[0]), 1) not in assigned
		        remaining_second=(read_name(second[0]), 2) not in assigned
		        if remaining_first and remaining_second:
		            forward.writelines(first)
		            rev_reads.writelines(second)
		        else:
		            for record, remaining, mate in [(first, remaining_first, 1), (second, remaining_second, 2)]:
		                if remaining:
		                    unpaired.writelines(("@" + read_name(record[0]) + "__MOSAIC_mate" + str(mate) + "\n", record[1], "+\n", record[3]))
		    for record in fastq(orphan):
		        name=read_name(record[0])
		        original=re.fullmatch(r"(.+)__MOSAIC_mate([12])", name)
		        identity=(original.group(1), int(original.group(2))) if original else (name, 0)
		        if identity not in assigned:
		            unpaired.writelines(record)
		PYTHON
		"""

rule isolate_unexplained_reads:
	input:
		forward=ISOLATE_MAPPING_DIR + "/{sample}/" + ISOLATE_STAGES[-1] + ".remaining_R1.fastq.gz",
		rev_reads=ISOLATE_MAPPING_DIR + "/{sample}/" + ISOLATE_STAGES[-1] + ".remaining_R2.fastq.gz",
		unpaired=ISOLATE_MAPPING_DIR + "/{sample}/" + ISOLATE_STAGES[-1] + ".remaining_U.fastq.gz",
	output:
		forward=ISOLATE_MAPPING_DIR + "/{sample}/unexplained_R1.fastq.gz",
		rev_reads=ISOLATE_MAPPING_DIR + "/{sample}/unexplained_R2.fastq.gz",
		unpaired=ISOLATE_MAPPING_DIR + "/{sample}/unexplained_unpaired.fastq.gz",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/isolate_unexplained_reads/sample={sample}.tsv"
	threads: 1
	shell:
		r"""
		ln -f -- {input.forward:q} {output.forward:q}
		ln -f -- {input.rev_reads:q} {output.rev_reads:q}
		ln -f -- {input.unpaired:q} {output.unpaired:q}
		"""

rule select_host_bacphlip_genomes:
	input:
		fasta=dirs_dict["HOST_DIR"] + "/prophages/{host}_prophages.fasta",
		evidence=dirs_dict["HOST_DIR"] + "/prophages/{host}_viral_evidence.tsv",
		checkv=dirs_dict["HOST_DIR"] + "/prophages/{host}_checkV/quality_summary.tsv",
	output:
		fasta=dirs_dict["HOST_DIR"] + "/prophages/{host}_bacphlip_input.fasta",
		eligibility=dirs_dict["HOST_DIR"] + "/prophages/{host}_bacphlip_eligibility.tsv",
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/select_host_bacphlip_genomes/host={host}.tsv"
	threads: 1
	shell:
		r"""
		python - {input.fasta:q} {input.evidence:q} {input.checkv:q} {output.fasta:q} {output.eligibility:q} <<-'PYTHON'
		import os
		import sys
		import pandas as pd
		from Bio import SeqIO

		fasta_path, evidence_path, checkv_path, output_fasta, eligibility_path=sys.argv[1:]
		evidence=pd.read_csv(evidence_path, sep="\t")
		checkv=pd.read_csv(checkv_path, sep="\t") if os.path.getsize(checkv_path) else pd.DataFrame(columns=["contig_id", "checkv_quality"])
		evidence=evidence.merge(checkv[["contig_id", "checkv_quality"]], left_on="viral_id", right_on="contig_id", how="left")
		phage=evidence["taxonomy"].fillna("").str.contains("Caudoviricetes", regex=False)
		complete=evidence["checkv_quality"].eq("Complete")
		evidence["bacphlip_assessment"]="eligible"
		evidence.loc[~phage, "bacphlip_assessment"]="not_assessed_outside_model_scope"
		evidence.loc[phage & ~complete, "bacphlip_assessment"]="not_assessed_incomplete_genome"
		eligible=set(evidence.loc[phage & complete, "viral_id"])
		with open(output_fasta, "w") as handle:
		    for record in SeqIO.parse(fasta_path, "fasta"):
		        if record.id in eligible:
		            SeqIO.write(record, handle, "fasta")
		evidence[["viral_id", "bacphlip_assessment"]].to_csv(eligibility_path, sep="\t", index=False)
		PYTHON
		"""


rule hosts_summary:
	input:
		checkm=dirs_dict["HOST_DIR"] + "/checkM_summary.csv",
		host_fastas=expand(dirs_dict["HOST_DIR"] + "/{host}.fasta", host=HOSTS),
		evidence=expand(dirs_dict["HOST_DIR"] + "/prophages/{host}_viral_evidence.tsv", host=HOSTS),
		checkv=expand(dirs_dict["HOST_DIR"] + "/prophages/{host}_checkV/quality_summary.tsv", host=HOSTS),
		eligibility=expand(dirs_dict["HOST_DIR"] + "/prophages/{host}_bacphlip_eligibility.tsv", host=HOSTS),
		bacphlip=expand(dirs_dict["ANNOTATION"] + "/host_viral_{host}_bacphlip.csv", host=HOSTS),
		circularity=expand(dirs_dict["VIRAL_DIR"] + "/{sample}_host_circularity.tot.tsv", sample=HOSTS),
	output:
		html=dirs_dict["PLOTS_DIR"] + "/08_hosts_summary.tot.html",
		hosts=dirs_dict["PLOTS_DIR"] + "/08_hosts_summary.tot.tsv",
		viral_regions=dirs_dict["PLOTS_DIR"] + "/08_host_viral_regions.tot.tsv",
		png=dirs_dict["PLOTS_DIR"] + "/08_hosts_prophages.tot.png",
		svg=dirs_dict["PLOTS_DIR"] + "/08_hosts_prophages.tot.svg",
	params:
		hosts=HOSTS,
	benchmark:
		dirs_dict["BENCHMARKS"] + "/hosts_summary/tot.tsv"
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/08_hosts_summary.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_hosts_summary.py.ipynb"


# Reuse full-read host mapping, not the sequential residual host/prophage stages.
rule host_prophage_activity_depth:
	input:
		bam=dirs_dict["MAPPING_DIR"] + "/HOST/bowtie2_{sample}_vs_{host}_filtered.bam",
	output:
		depth=dirs_dict["MAPPING_DIR"] + "/HOST/ACTIVITY/{sample}_vs_{host}.basecov.tsv.gz",
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/host_prophage_activity_depth/host={host}__sample={sample}.tsv"
	threads: 2
	shell:
		r"""
		# One primary placement per read; secondary/supplementary records do not add depth.
		samtools view -@ {threads} -u -F 2308 {input.bam:q} | \
			bedtools genomecov -dz -ibam stdin | gzip -c > {output.depth:q}
		"""


rule host_prophage_activity:
	input:
		regions=dirs_dict["PLOTS_DIR"] + "/08_host_viral_regions.tot.tsv",
		assignments=ALL_ASSEMBLED_DIR + "/phage_isolates.tot/sample_host_assignments.tsv",
		host_fastas=expand(dirs_dict["HOST_DIR"] + "/{host}.fasta", host=HOSTS),
		depth=lambda wc: [dirs_dict["MAPPING_DIR"] + "/HOST/ACTIVITY/" + sample + "_vs_" + ISOLATE_HOST_ASSIGNMENTS[sample] + ".basecov.tsv.gz"
			for sample in SAMPLES if sample in ISOLATE_HOST_ASSIGNMENTS],
	output:
		table=dirs_dict["PLOTS_DIR"] + "/08_host_prophage_activity.tot.tsv",
		png=dirs_dict["PLOTS_DIR"] + "/08_host_prophage_activity.tot.png",
		svg=dirs_dict["PLOTS_DIR"] + "/08_host_prophage_activity.tot.svg",
		figures=directory(dirs_dict["PLOTS_DIR"] + "/08_host_prophage_activity.tot"),
	params:
		samples=[sample for sample in SAMPLES if sample in ISOLATE_HOST_ASSIGNMENTS],
		hosts=HOSTS,
		min_ratio=float(config.get("host_prophage_activity_min_ratio", 2.0)),
		min_cohen_d=float(config.get("host_prophage_activity_min_cohen_d", 0.70)),
		min_mean_depth=float(config.get("host_prophage_activity_min_mean_depth", 1.0)),
		min_breadth_percent=float(config.get("host_prophage_activity_min_breadth_percent", 50)),
		mask_bp=int(config.get("host_prophage_activity_mask_bp", 150)),
		min_length_bp=int(config.get("host_prophage_activity_min_length_bp", 1000)),
		plot_flank_bp=int(config.get("host_prophage_activity_plot_flank_bp", 50000)),
	message:
		"Checking host prophage coverage enrichment with PropagAtE activity defaults"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/host_prophage_activity/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/08_host_prophage_activity.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_host_prophage_activity.py.ipynb"


# Coverage-first exploration outside geNomad predictions; no new mapping or viral calls.
rule host_enriched_regions:
	input:
		regions=dirs_dict["PLOTS_DIR"] + "/08_host_viral_regions.tot.tsv",
		assignments=ALL_ASSEMBLED_DIR + "/phage_isolates.tot/sample_host_assignments.tsv",
		host_fastas=expand(dirs_dict["HOST_DIR"] + "/{host}.fasta", host=HOSTS),
		depth=lambda wc: [dirs_dict["MAPPING_DIR"] + "/HOST/ACTIVITY/" + sample + "_vs_" + ISOLATE_HOST_ASSIGNMENTS[sample] + ".basecov.tsv.gz"
			for sample in SAMPLES if sample in ISOLATE_HOST_ASSIGNMENTS],
	output:
		table=dirs_dict["PLOTS_DIR"] + "/08_host_enriched_regions.tot.tsv",
		summary=dirs_dict["PLOTS_DIR"] + "/08_host_enriched_regions_summary.tot.tsv",
		figures=directory(dirs_dict["PLOTS_DIR"] + "/08_host_enriched_regions.tot"),
	params:
		samples=[sample for sample in SAMPLES if sample in ISOLATE_HOST_ASSIGNMENTS],
		hosts=HOSTS,
		min_ratio=float(config.get("host_enriched_regions_min_ratio", 2.0)),
		min_cohen_d=float(config.get("host_enriched_regions_min_cohen_d", 0.70)),
		min_mean_depth=float(config.get("host_enriched_regions_min_mean_depth", 5.0)),
		breadth_min_depth=float(config.get("host_enriched_regions_breadth_min_depth", 1.0)),
		min_breadth_percent=float(config.get("host_enriched_regions_min_breadth_percent", 80)),
		min_local_ratio=float(config.get("host_enriched_regions_min_local_ratio", 2.0)),
		mask_bp=int(config.get("host_prophage_activity_mask_bp", 150)),
		bin_bp=int(config.get("host_enriched_regions_bin_bp", 1000)),
		min_length_bp=int(config.get("host_enriched_regions_min_length_bp", 5000)),
		local_flank_bp=int(config.get("host_enriched_regions_local_flank_bp", 10000)),
		plot_flank_bp=int(config.get("host_enriched_regions_plot_flank_bp", 200000)),
	message:
		"Exploring host coverage enrichment outside geNomad predictions"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/host_enriched_regions/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/08_host_enriched_regions.tot.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/08_host_enriched_regions.py.ipynb"


rule extract_host_enriched_regions:
	input:
		regions=dirs_dict["PLOTS_DIR"] + "/08_host_enriched_regions.tot.tsv",
		viral_regions=dirs_dict["PLOTS_DIR"] + "/08_host_viral_regions.tot.tsv",
		host_fastas=expand(dirs_dict["HOST_DIR"] + "/{host}.fasta", host=HOSTS),
	output:
		fasta=dirs_dict["HOST_DIR"] + "/enriched_regions/host_enriched_regions.tot.fasta",
		table=dirs_dict["HOST_DIR"] + "/enriched_regions/host_enriched_regions.tot.tsv",
	params:
		hosts=HOSTS,
		adjacency_bp=int(config.get("host_enriched_regions_prophage_adjacency_bp", 1000)),
	message:
		"Exporting unique enriched host loci and their supporting samples"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/extract_host_enriched_regions/tot.tsv"
	threads: 1
	run:
		from pathlib import Path
		import pandas as pd
		from Bio import SeqIO
		from Bio.SeqRecord import SeqRecord

		calls=pd.read_csv(input.regions, sep="\t", dtype={"host": str, "parent_contig": str, "sample": str})
		viral_regions=pd.read_csv(input.viral_regions, sep="\t", dtype={"host": str, "parent_contig": str})
		loci=[]
		# Merge overlapping calls across samples, not their much larger plotting windows.
		for (host, parent), group in calls.groupby(["host", "parent_contig"], sort=False):
			current=None
			for row in group.sort_values(["start", "end", "sample"]).itertuples():
				if current is None or row.start > current["end"]:
					current=dict(host=host, parent_contig=parent, start=int(row.start), end=int(row.end), calls=[])
					loci.append(current)
				current["end"]=max(current["end"], int(row.end))
				current["calls"].append(row)
		columns=["region_id", "host", "parent_contig", "start", "end", "length_bp", "n_samples", "samples", "n_sample_calls", "source_intervals",
			"next_to_prophage", "nearest_prophage_id", "nearest_prophage_distance_bp", "prophage_adjacency_threshold_bp"]
		rows=[]
		Path(output.fasta).parent.mkdir(parents=True, exist_ok=True)
		with open(output.fasta, "w") as handle:
			for host, path in zip(params.hosts, input.host_fastas):
				wanted={}
				for locus in loci:
					if locus["host"] == host:
						wanted.setdefault(locus["parent_contig"], []).append(locus)
				with open(path) as source:
					records=list(SeqIO.parse(source, "fasta"))
				for record in records:
					for locus in wanted.get(record.id, []):
						start, end=locus["start"], locus["end"]
						region_id=host + "_" + record.id + "_enriched_" + str(start) + "_" + str(end)
						samples=sorted({call.sample for call in locus["calls"]})
						# Reported coordinates are 1-based inclusive; no context flanks are exported.
						sequence=record.seq[start-1:end]
						neighbours=[]
						for viral in viral_regions[viral_regions["host"].eq(host) & viral_regions["parent_contig"].eq(record.id) &
							viral_regions["evidence_scope"].eq("embedded_prophage")].itertuples():
							left, right=sorted([int(viral.start), int(viral.end)])
							neighbours.append((max(left - end - 1, start - right - 1, 0), viral.viral_id))
						nearest=min(neighbours) if neighbours else None
						distance=nearest[0] if nearest else "not reported"
						nearest_id=nearest[1] if nearest else "not reported"
						adjacent=nearest is not None and nearest[0] <= params.adjacency_bp
						SeqIO.write(SeqRecord(sequence, id=region_id,
							description="host=" + host + " contig=" + record.id + " length=" + str(len(sequence)) + " samples=" + ";".join(samples) +
							" next_to_prophage=" + str(adjacent) + " nearest_prophage=" + nearest_id.replace(" ", "_") +
							" prophage_gap_bp=" + str(distance).replace(" ", "_") + " adjacency_threshold_bp=" + str(params.adjacency_bp)), handle, "fasta")
						rows.append(dict(region_id=region_id, host=host, parent_contig=record.id,
							start=start, end=end, length_bp=len(sequence), n_samples=len(samples), samples=";".join(samples),
							n_sample_calls=len(locus["calls"]), source_intervals=";".join(str(a) + "-" + str(b)
								for a, b in sorted({(call.start, call.end) for call in locus["calls"]})),
							next_to_prophage=adjacent, nearest_prophage_id=nearest_id, nearest_prophage_distance_bp=distance,
							prophage_adjacency_threshold_bp=params.adjacency_bp))
		pd.DataFrame(rows, columns=columns).to_csv(output.table, sep="\t", index=False)


rule host_enriched_region_metadata:
	input:
		regions=dirs_dict["HOST_DIR"] + "/enriched_regions/host_enriched_regions.tot.tsv",
		genomad=dirs_dict["ANNOTATION"] + "/geNomad_host_enriched_regions.tot",
		checkv=dirs_dict["HOST_DIR"] + "/enriched_regions/host_enriched_regions.tot_checkV/quality_summary.tsv",
		phold=dirs_dict["ANNOTATION"] + "/phold_host_enriched_regions.tot",
		taxmyphage=dirs_dict["ANNOTATION"] + "/taxmyphage_host_enriched_regions.tot",
		refseq=dirs_dict["ANNOTATION"] + "/blast_output_ViralRefSeq_host_enriched_regions.tot.csv",
		metavr=[dirs_dict["ANNOTATION"] + "/blast_output_METAVR_host_enriched_regions.tot.csv"] if METAVR_blast else [],
		metavr_metadata=[os.path.join(config["METAVR_db"], "METAVR_main_table.parquet")] if METAVR_blast else [],
	output:
		metadata=dirs_dict["HOST_DIR"] + "/enriched_regions/host_enriched_regions_metadata.tot.tsv",
	message:
		"Combining enriched-locus provenance with Phold, geNomad, CheckV and reference taxonomy"
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/host_enriched_region_metadata/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	shell:
		r"""
		python - {input.regions:q} {input.genomad:q} {input.checkv:q} {input.phold:q} {input.taxmyphage:q} {input.refseq:q} {output.metadata:q} {input.metavr:q} {input.metavr_metadata:q} <<-'PYTHON'
		import sys
		from types import SimpleNamespace
		input=SimpleNamespace(regions=sys.argv[1], genomad=sys.argv[2], checkv=sys.argv[3], phold=sys.argv[4],
		    taxmyphage=sys.argv[5], refseq=sys.argv[6], metavr=sys.argv[8:9], metavr_metadata=sys.argv[9:10])
		output=SimpleNamespace(metadata=sys.argv[7])
		from pathlib import Path
		import pandas as pd

		def table(path):
		    path=Path(path)
		    return pd.read_csv(path, sep="\t") if path.is_file() and path.stat().st_size else pd.DataFrame()

		metadata=table(input.regions).set_index("region_id")
		prefix="host_enriched_regions.tot"
		folder=Path(input.genomad)
		virus=table(folder / (prefix + "_summary") / (prefix + "_virus_summary.tsv"))
		plasmid=table(folder / (prefix + "_summary") / (prefix + "_plasmid_summary.tsv"))
		scores=table(folder / (prefix + "_aggregated_classification") / (prefix + "_aggregated_classification.tsv"))
		if "seq_name" in scores:
		    scores=scores.drop_duplicates("seq_name").set_index("seq_name")
		    metadata=metadata.join(scores.reindex(metadata.index).add_prefix("geNomad_contig_"))
		metadata["geNomad_classification"]="no viral/plasmid call"
		metadata["geNomad_evidence_scope"]="not reported"
		for biology, frame in [("plasmid", plasmid), ("virus", virus)]:
		    if "seq_name" not in frame:
		        continue
		    frame=frame.copy()
		    frame["region_id"]=frame["seq_name"].str.split("|provirus_", regex=False).str[0]
		    for identifier, hits in frame.groupby("region_id", sort=False):
		        if identifier not in metadata.index:
		            continue
		        whole=hits["seq_name"].eq(identifier).any()
		        metadata.loc[identifier, "geNomad_classification"]=biology if whole else "viral region"
		        metadata.loc[identifier, "geNomad_evidence_scope"]="whole contig" if whole else "provirus region"
		        for column in frame.columns.drop("region_id"):
		            values=hits[column].dropna().astype(str).drop_duplicates()
		            field="geNomad_" + column if column in ["virus_score", "plasmid_score", "taxonomy"] else "geNomad_" + biology + "_" + column
		            metadata.loc[identifier, field]=" | ".join(values) or "not reported"
		for biology in ["virus", "plasmid"]:
		    field="geNomad_" + biology + "_score"
		    contig_field="geNomad_contig_" + biology + "_score"
		    if field not in metadata:
		        metadata[field]=pd.NA
		    if contig_field in metadata:
		        metadata[field]=metadata[field].where(metadata[field].notna(), metadata[contig_field])

		for path, key, name in [(input.checkv, "contig_id", "CheckV"),
		    (Path(input.taxmyphage) / "Summary_taxonomy.tsv", "Genome", "TaxMyPhage")]:
		    frame=table(path)
		    if key in frame:
		        metadata=metadata.join(frame.drop_duplicates(key).set_index(key).reindex(metadata.index).add_prefix(name + "_"))
		functions=table(Path(input.phold) / (prefix + "_all_cds_functions.tsv"))
		if all(field in functions for field in ["Contig", "Description", "Count"]):
		    functions=functions.pivot(index="Contig", columns="Description", values="Count")
		    functions.columns=["Phold_" + str(name).replace(" ", "_") for name in functions.columns]
		    metadata=metadata.join(functions.reindex(metadata.index))
		cds=table(Path(input.phold) / (prefix + "_per_cds_predictions.tsv"))
		if "contig_id" in cds:
		    metadata["Phold_CDS_count"]=cds.groupby("contig_id").size().reindex(metadata.index, fill_value=0)
		    for field in ["function", "product", "annotation_confidence"]:
		        if field in cds:
		            metadata["Phold_" + field + "s"]=cds.groupby("contig_id")[field].agg(
		                lambda values: " | ".join(sorted(set(values.dropna().astype(str)))))

		blast_columns=["qseqid", "sseqid", "description", "qstart", "qend", "qlen", "slen", "qcovs", "evalue", "alignment_length", "pident"]
		for paths, name in [([input.refseq], "RefSeq"), (input.metavr, "METAVR")]:
		    metadata[name + "_assessment"]="run" if paths else "not requested"
		    for path in paths:
		        frame=pd.read_csv(path, sep="\t", header=None, names=blast_columns) if Path(path).stat().st_size else pd.DataFrame(columns=blast_columns)
		        frame=frame[frame["qseqid"].ne("qseqid")].copy()
		        for field in blast_columns[3:]:
		            frame[field]=pd.to_numeric(frame[field], errors="coerce")
		        # Use the same best-hit ordering as the existing all-assembled metadata.
		        frame=frame.sort_values(["evalue", "qcovs", "pident", "alignment_length", "sseqid"],
		            ascending=[True, False, False, False, True], kind="stable")
		        metadata=metadata.join(frame.drop_duplicates("qseqid").set_index("qseqid").reindex(metadata.index).add_prefix(name + "_"))
		if input.metavr_metadata:
		    accessions=metadata["METAVR_sseqid"].dropna().str.split("|").str[0].unique().tolist()
		    metavr=pd.read_parquet(input.metavr_metadata[0], columns=["uvig", "quality", "genome_type", "ictv_taxonomy"],
		        filters=[("uvig", "in", accessions)]).drop_duplicates("uvig").set_index("uvig") if accessions else pd.DataFrame()
		    for field in ["quality", "genome_type", "ictv_taxonomy"]:
		        if field in metavr:
		            metadata["METAVR_" + field]=metadata["METAVR_sseqid"].fillna("").str.split("|").str[0].map(metavr[field])
		metadata.reset_index().to_csv(output.metadata, sep="\t", index=False, na_rep="not reported")
		PYTHON
		"""
