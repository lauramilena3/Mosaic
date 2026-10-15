# Optional retained-contig network. Existing retention, clustering and mapping are unchanged.
rule isolate_network:
	input:
		dirs_dict["PLOTS_DIR"] + "/08_isolate_network.tot.html",


rule select_isolate_network_contigs:
	input:
		metadata=ALL_ASSEMBLED_DIR + "/phage_isolates.tot/all_contig_metadata.tsv",
		fasta=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs.tot.fasta",
		summary=dirs_dict["PLOTS_DIR"] + "/08_phage_isolates_summary.tot.csv",
	output:
		fasta=ISOLATE_NETWORK_PREFIX + ".own.fasta",
		nodes=ISOLATE_NETWORK_PREFIX + ".own_nodes.tsv",
	message:
		"Collecting retained isolate contigs and existing evidence for the network"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/select_isolate_network_contigs/tot.tsv"
	threads: 1
	run:
		import hashlib
		import pandas as pd
		from Bio import SeqIO
		from Bio.SeqRecord import SeqRecord

		metadata=pd.read_csv(input.metadata, sep="\t")
		nodes=metadata.loc[metadata.retained.astype(str).str.lower().eq("true")].copy()
		summary=pd.read_csv(input.summary)
		nodes["decision"]=nodes["sample"].map(summary.set_index("sample")["decision"])
		nodes["node_id"]="MOSAIC:" + nodes["original_id"]
		nodes["source"]="MOSAIC"
		nodes["accession"]=nodes["original_id"]
		nodes["description"]=nodes["original_id"]
		wanted=set(nodes.original_id)
		comparison_ids={}
		sequences={}
		with open(input.fasta) as handle:
			for record in SeqIO.parse(handle, "fasta"):
				if record.id not in wanted:
					continue
				sequence=record.seq.upper()
				canonical=min(str(sequence), str(sequence.reverse_complement()))
				comparison_id="own_" + hashlib.sha256(canonical.encode()).hexdigest()[:24]
				comparison_ids[record.id]=comparison_id
				sequences[comparison_id]=SeqRecord(sequence, id=comparison_id, description="")
		nodes["comparison_id"]=nodes["original_id"].map(comparison_ids)
		# Exact sequence hashing avoids treating shorter MMseqs-contained members as identical.
		SeqIO.write([sequences[key] for key in sorted(sequences)], output.fasta, "fasta-2line")
		nodes.sort_values(["sample", "original_id"]).to_csv(output.nodes, sep="\t", index=False)


rule select_isolate_network_references:
	input:
		nodes=ISOLATE_NETWORK_PREFIX + ".own_nodes.tsv",
		parsers=os.path.join(workflow.basedir, "notebooks/08_phage_isolates_summary.py.ipynb"),
		blast=lambda wc: [dirs_dict["ANNOTATION"] + "/blast_output_" + {"RefSeq": "ViralRefSeq", "METAVR": "METAVR"}[source] + "_phage_isolates_cluster_representatives.tot.csv" for source in ISOLATE_NETWORK_SOURCES],
	output:
		hits=ISOLATE_NETWORK_PREFIX + ".reference_hits.tsv",
	params:
		sources=ISOLATE_NETWORK_SOURCES,
		max_references=int(config.get("isolate_network_max_references_per_query", 3)),
		min_identity=float(config.get("isolate_network_min_reference_identity", 70)),
		min_query_coverage=float(config.get("isolate_network_min_reference_query_coverage", 50)),
	message:
		"Selecting a bounded set of RefSeq and optional METAVR relatives"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/select_isolate_network_references/tot.tsv"
	threads: 1
	run:
		import ast
		import json
		import re
		import numpy as np
		import pandas as pd
		from pathlib import Path

		# Reuse the summary's overlap-aware BLAST parser without running a notebook or worker.
		functions=["_read_table", "_to_numeric", "_merge_intervals", "_uncovered_intervals",
			"_blast_alignment_summary", "_read_blast_alignments", "_group_blast_alignments"]
		state={"pd": pd, "np": np, "Path": Path, "re": re}
		for cell in json.loads(Path(input.parsers).read_text())["cells"]:
			if cell["cell_type"] != "code":
				continue
			source="".join(cell["source"])
			if "BLAST_FIELDS =" not in source and not any("def " + name + "(" in source for name in functions):
				continue
			for node in ast.parse(source).body:
				if (isinstance(node, ast.FunctionDef) and node.name in functions or
					isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == "BLAST_FIELDS" for target in node.targets)):
					exec(compile(ast.Module(body=[node], type_ignores=[]), input.parsers, "exec"), state)
		queries=set(pd.read_csv(input.nodes, sep="\t").cluster_rep)
		selected=[]
		for source, path in zip(params.sources, input.blast):
			frame=state["_read_blast_alignments"](path)
			frame=frame.loc[frame["query"].isin(queries)].copy()
			if frame.empty or params.max_references <= 0:
				continue
			pairs=state["_group_blast_alignments"](frame, ["query", "subject", "subject_title", "q_len", "s_len"])
			pairs=pairs.loc[pairs.pident.ge(params.min_identity) & pairs.query_coverage_percent.ge(params.min_query_coverage)].copy()
			pairs=pairs.sort_values(["query", "query_coverage_percent", "pident", "aligned_bases", "subject"], ascending=[True, False, False, False, True])
			pairs=pairs.groupby("query", sort=False).head(params.max_references).copy()
			pairs["source"]=source
			pairs["accession"]=pairs["subject"].str.replace(r"^(?:ref|gb|emb|dbj)\|([^|]+)\|$", r"\1", regex=True)
			pairs["selection_scope"]="existing cluster-representative BLAST; discovery only, not evidence inherited by members"
			selected.append(pairs)
		columns=["query", "subject", "subject_title", "q_len", "s_len", "pident", "query_coverage_percent", "subject_coverage_percent", "source", "accession", "selection_scope"]
		hits=pd.concat(selected, ignore_index=True) if selected else pd.DataFrame(columns=columns)
		hits.to_csv(output.hits, sep="\t", index=False)


rule prepare_isolate_network_genomes:
	input:
		fasta=ISOLATE_NETWORK_PREFIX + ".own.fasta",
		nodes=ISOLATE_NETWORK_PREFIX + ".own_nodes.tsv",
		hits=ISOLATE_NETWORK_PREFIX + ".reference_hits.tsv",
		refseq=lambda wc: [config["RefSeqViral_db"]] if "RefSeq" in ISOLATE_NETWORK_SOURCES else [],
		metavr=lambda wc: [os.path.join(config["METAVR_db"], "METAVR_UViG_blastdb")] if "METAVR" in ISOLATE_NETWORK_SOURCES else [],
	output:
		fasta=ISOLATE_NETWORK_PREFIX + ".fasta",
		nodes=ISOLATE_NETWORK_PREFIX + ".nodes.tsv",
		missing=ISOLATE_NETWORK_PREFIX + ".missing_references.tsv",
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	message:
		"Preparing the retained genomes and selected reference sequences for Mash"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/prepare_isolate_network_genomes/tot.tsv"
	threads: 1
	shell:
		r"""
		python - {input.fasta:q} {input.nodes:q} {input.hits:q} {output.fasta:q} {output.nodes:q} {output.missing:q} \
			'{input.refseq}' '{input.metavr}' <<-'PYTHON'
		import io
		import re
		import sys
		import subprocess
		from pathlib import Path
		import pandas as pd
		from Bio import SeqIO
		from Bio.SeqRecord import SeqRecord

		own_fasta, own_nodes, hit_path, out_fasta, out_nodes, missing_path, refseq, metavr=sys.argv[1:]
		nodes=pd.read_csv(own_nodes, sep="\t")
		hits=pd.read_csv(hit_path, sep="\t")
		references=hits.drop_duplicates(["source", "accession"])
		with open(own_fasta) as handle:
		    sequences=list(SeqIO.parse(handle, "fasta"))
		found={{}}
		if refseq:
		    wanted=set(references.loc[references.source.eq("RefSeq"), "accession"])
		    with open(refseq) as handle:
		        for record in SeqIO.parse(handle, "fasta"):
		            accession=re.sub(r"^(?:ref|gb|emb|dbj)\|([^|]+)\|$", r"\1", record.id)
		            if accession in wanted:
		                found[("RefSeq", accession)]=record
		                wanted.remove(accession)
		            if not wanted:
		                break
		if metavr:
		    database=str(Path(metavr) / "METAVR_UViG.blastdb")
		    for row in references.loc[references.source.eq("METAVR")].itertuples():
		        result=subprocess.run(["blastdbcmd", "-db", database, "-entry", row.subject, "-outfmt", "%f"], capture_output=True, text=True)
		        records=list(SeqIO.parse(io.StringIO(result.stdout), "fasta")) if result.returncode == 0 else []
		        if records:
		            found[("METAVR", row.accession)]=records[0]
		extra=[]
		missing=[]
		for row in references.itertuples():
		    key=(row.source, row.accession)
		    if key not in found:
		        missing.append(dict(source=row.source, accession=row.accession, subject=row.subject, reason="sequence not retrieved from the configured database"))
		        continue
		    record=found[key]
		    node_id=row.source + ":" + row.accession
		    comparison_id="ref_" + row.source + "_" + re.sub(r"[^A-Za-z0-9_.-]", "_", row.accession)
		    sequences.append(SeqRecord(record.seq, id=comparison_id, description=""))
		    extra.append(dict(node_id=node_id, comparison_id=comparison_id, original_id=row.accession,
		        source=row.source, accession=row.accession, description=row.subject_title, length_bp=len(record),
		        retained=False, genomad_classification="not reported", CheckV_checkv_quality="not reported",
		        genomad_evidence_scope="not reported", selection_queries=";".join(sorted(set(hits.loc[hits.source.eq(row.source) & hits.accession.eq(row.accession), "query"])))))
		SeqIO.write(sequences, out_fasta, "fasta-2line")
		pd.concat([nodes, pd.DataFrame(extra)], ignore_index=True).to_csv(out_nodes, sep="\t", index=False)
		pd.DataFrame(missing, columns=["source", "accession", "subject", "reason"]).to_csv(missing_path, sep="\t", index=False)
		PYTHON
		"""


rule mash_isolate_network:
	input:
		fasta=ISOLATE_NETWORK_PREFIX + ".fasta",
	output:
		sketch=ISOLATE_NETWORK_PREFIX + ".mash.msh",
		distances=ISOLATE_NETWORK_PREFIX + ".mash_distances.tsv",
		version=ISOLATE_NETWORK_PREFIX + ".mash_version.txt",
	params:
		prefix=ISOLATE_NETWORK_PREFIX + ".mash",
		kmer_size=int(config.get("isolate_network_kmer_size", 15)),
		sketch_size=int(config.get("isolate_network_sketch_size", 25000)),
		seed=int(config.get("isolate_network_seed", 42)),
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	message:
		"Calculating sketch-based Mash distances for the isolate network"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mash_isolate_network/tot.tsv"
	threads: int(config.get("isolate_network_threads", 8))
	shell:
		r"""
		mash --version > {output.version:q}
		printf 'reference\tquery\tmash_distance\tp_value\tmatching_hashes\n' > {output.distances:q}
		if [ -s {input.fasta:q} ]; then
			mash sketch -i -k {params.kmer_size} -s {params.sketch_size} -S {params.seed} -p {threads} \
				-o {params.prefix:q} {input.fasta:q}
			mash dist -p {threads} {output.sketch:q} {output.sketch:q} >> {output.distances:q}
		else
			: > {output.sketch:q}
		fi
		"""


rule isolate_network_report:
	input:
		nodes=ISOLATE_NETWORK_PREFIX + ".nodes.tsv",
		hits=ISOLATE_NETWORK_PREFIX + ".reference_hits.tsv",
		missing=ISOLATE_NETWORK_PREFIX + ".missing_references.tsv",
		distances=ISOLATE_NETWORK_PREFIX + ".mash_distances.tsv",
		version=ISOLATE_NETWORK_PREFIX + ".mash_version.txt",
		template=os.path.join(workflow.basedir, "templates/08_isolate_network.html"),
		cytoscape=os.path.join(workflow.basedir, "tools/cytoscape-3.33.1/cytoscape.min.js"),
	output:
		html=dirs_dict["PLOTS_DIR"] + "/08_isolate_network.tot.html",
		nodes=dirs_dict["PLOTS_DIR"] + "/08_isolate_network.tot/nodes.tsv",
		edges=dirs_dict["PLOTS_DIR"] + "/08_isolate_network.tot/edges.tsv",
		hits=dirs_dict["PLOTS_DIR"] + "/08_isolate_network.tot/reference_hits.tsv",
		missing=dirs_dict["PLOTS_DIR"] + "/08_isolate_network.tot/missing_references.tsv",
		manifest=dirs_dict["PLOTS_DIR"] + "/08_isolate_network.tot/provenance.json",
	params:
		default_distance=float(config.get("isolate_network_default_distance", 0.15)),
		kmer_size=int(config.get("isolate_network_kmer_size", 15)),
		sketch_size=int(config.get("isolate_network_sketch_size", 25000)),
		seed=int(config.get("isolate_network_seed", 42)),
		sources=ISOLATE_NETWORK_SOURCES,
		max_references=int(config.get("isolate_network_max_references_per_query", 3)),
		min_identity=float(config.get("isolate_network_min_reference_identity", 70)),
		min_query_coverage=float(config.get("isolate_network_min_reference_query_coverage", 50)),
	message:
		"Writing the standalone interactive isolate genome network"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/isolate_network_report/tot.tsv"
	threads: 1
	run:
		import hashlib
		import itertools
		import json
		import shutil
		from datetime import datetime, timezone
		from pathlib import Path
		import pandas as pd

		nodes=pd.read_csv(input.nodes, sep="\t")
		comparisons=pd.read_csv(input.distances, sep="\t")
		self_comparisons=comparisons.loc[comparisons.reference.eq(comparisons["query"])].set_index("reference")
		comparisons=comparisons.loc[comparisons.mash_distance.lt(1) & comparisons.matching_hashes.str.split("/").str[0].astype(int).gt(0)]
		members=nodes.groupby("comparison_id")["node_id"].agg(list).to_dict()
		edges=[]
		seen=set()
		for row in comparisons.itertuples():
			if row.reference == row.query:
				continue
			key=tuple(sorted([row.reference, row.query]))
			if key in seen:
				continue
			seen.add(key)
			shared, total=map(int, row.matching_hashes.split("/"))
			for source, target in itertools.product(members.get(row.reference, []), members.get(row.query, [])):
				edges.append(dict(source=source, target=target, mash_distance=row.mash_distance,
					p_value=row.p_value, matching_hashes=row.matching_hashes,
					shared_hashes=shared, sketch_hashes=total, evidence="Mash"))
		for comparison_id, group in members.items():
			measurement=self_comparisons.loc[comparison_id] if comparison_id in self_comparisons.index else None
			matching=measurement.matching_hashes if measurement is not None else None
			for source, target in itertools.combinations(group, 2):
				edges.append(dict(source=source, target=target, mash_distance=0,
					p_value=measurement.p_value if measurement is not None else None, matching_hashes=matching,
					shared_hashes=int(matching.split("/")[0]) if matching else None,
					sketch_hashes=int(matching.split("/")[1]) if matching else None,
					evidence="identical sequence or reverse complement; compared once"))
		columns=["source", "target", "mash_distance", "p_value", "matching_hashes", "shared_hashes", "sketch_hashes", "evidence"]
		edges=pd.DataFrame(edges, columns=columns)
		Path(output.nodes).parent.mkdir(parents=True, exist_ok=True)
		nodes.to_csv(output.nodes, sep="\t", index=False)
		edges.to_csv(output.edges, sep="\t", index=False)
		shutil.copyfile(input.hits, output.hits)
		shutil.copyfile(input.missing, output.missing)
		missing=pd.read_csv(input.missing, sep="\t")
		manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),
			distance_metric="Mash distance (sketch-based; lower is closer)", mash_version=Path(input.version).read_text().strip(),
			kmer_size=params.kmer_size, sketch_size=params.sketch_size, seed=params.seed,
			existing_clustering="95% ANI, 85% target coverage; unchanged",
			retention="existing retained flag; no geNomad filtering", sources=list(params.sources),
			max_references_per_query=params.max_references, min_reference_identity=params.min_identity,
			min_reference_query_coverage=params.min_query_coverage, default_distance=params.default_distance,
			own_contigs=int(nodes.source.eq("MOSAIC").sum()), reference_sequences=int(nodes.source.ne("MOSAIC").sum()),
			distinct_comparison_sequences=int(nodes.comparison_id.nunique()), missing_references=len(missing),
			mash_distances=str(input.distances), cytoscape_version="3.33.1",
			input_sha256={str(path): hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in [input.nodes, input.hits, input.missing, input.distances, input.version]})
		Path(output.manifest).write_text(json.dumps(manifest, indent=2) + "\n")
		payload=dict(nodes=json.loads(nodes.to_json(orient="records")), edges=json.loads(edges.to_json(orient="records")), manifest=manifest)
		data=json.dumps(payload, allow_nan=False).replace("<", "\\u003c")
		html=Path(input.template).read_text().replace("__MOSAIC_NETWORK_DATA__", data)
		html=html.replace("__CYTOSCAPE_LIBRARY__", Path(input.cytoscape).read_text())
		Path(output.html).write_text(html)
