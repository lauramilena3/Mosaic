# Within-cluster VIRIDIC comparisons reuse the existing retention and fragment evidence.
# The checkpoint creates one job per eligible cluster, not a whole-catalogue comparison.
checkpoint select_isolate_viridic_clusters:
	input:
		metadata=ALL_ASSEMBLED_DIR + "/phage_isolates.{sampling}/all_contig_metadata.tsv",
		fasta=ALL_ASSEMBLED_DIR + "/phage_isolates_contigs.{sampling}.fasta",
	output:
		fasta_dir=directory(ISOLATE_VIRIDIC_INPUTS + "/INPUTS"),
		clusters=ISOLATE_VIRIDIC_INPUTS + "/clusters.tsv",
		members=ISOLATE_VIRIDIC_INPUTS + "/members.tsv",
	params:
		complete_only=config_bool("isolate_viridic_complete_only", False),
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	message:
		"Selecting retained non-fragment genomes for within-cluster VIRIDIC"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/select_isolate_viridic_clusters/sampling={sampling}.tsv"
	threads: 1
	shell:
		r"""
		python - {input.metadata:q} {input.fasta:q} {output.fasta_dir:q} {output.clusters:q} {output.members:q} \
			{params.complete_only} <<-'PYTHON'
		import hashlib
		import sys
		from pathlib import Path
		import pandas as pd
		from Bio import SeqIO

		metadata_path, fasta_path, fasta_dir, clusters_path, members_path, complete_only=sys.argv[1:]
		metadata=pd.read_csv(metadata_path, sep="\t")
		retained=metadata.loc[metadata.retained.astype(str).str.lower().eq("true")].copy()
		selected=retained.loc[~retained.fragment_candidate.astype(str).str.lower().eq("true")].copy()
		wanted=set(selected.original_id)
		records={{}}
		hashes={{}}
		with open(fasta_path) as handle:
		    for record in SeqIO.parse(handle, "fasta"):
		        if record.id in wanted:
		            sequence=record.seq.upper()
		            canonical=min(str(sequence), str(sequence.reverse_complement()))
		            hashes[record.id]=hashlib.sha256(canonical.encode()).hexdigest()
		            records[record.id]=record
		selected["sequence_sha256"]=[hashes[contig] for contig in selected.original_id]
		selected["is_cluster_representative"]=selected.original_id.eq(selected.cluster_rep)
		quality_order={{"Complete": 0, "High-quality": 1, "Medium-quality": 2, "Low-quality": 3, "Not-determined": 4}}
		selected["quality_order"]=selected.CheckV_checkv_quality.map(quality_order).fillna(5)
		selected["numeric_completeness"]=pd.to_numeric(selected.CheckV_completeness, errors="coerce")
		selected=selected.sort_values(["cluster_rep", "is_cluster_representative", "quality_order",
		    "numeric_completeness", "sample", "original_id"], ascending=[True, False, True, False, True, True])
		comparisons=selected.drop_duplicates(["cluster_rep", "sequence_sha256"])
		comparison_ids={{(row.cluster_rep, row.sequence_sha256): row.original_id for row in comparisons.itertuples()}}
		selected["comparison_id"]=[comparison_ids[(row.cluster_rep, row.sequence_sha256)] for row in selected.itertuples()]
		selected["is_comparison_record"]=selected.original_id.eq(selected.comparison_id)
		selected=selected.drop(columns=["quality_order", "numeric_completeness"])

		clusters=pd.DataFrame({{"cluster_rep": sorted(metadata.cluster_rep.unique())}})
		qualities=metadata.set_index("original_id").CheckV_checkv_quality
		clusters["representative_CheckV_quality"]=clusters.cluster_rep.map(qualities)
		for column, counts in [
		    ("number_retained_contigs", retained.groupby("cluster_rep").size()),
		    ("number_nonfragment_contigs", selected.groupby("cluster_rep").size()),
		    ("number_distinct_sequences", selected.groupby("cluster_rep").sequence_sha256.nunique()),
		    ("number_samples", selected.groupby("cluster_rep")["sample"].nunique())]:
		    clusters[column]=clusters.cluster_rep.map(counts).fillna(0).astype(int)
		clusters["selected"]=clusters.number_distinct_sequences.ge(2)
		clusters["selection_reason"]="selected"
		clusters.loc[~clusters.selected, "selection_reason"]="fewer_than_two_distinct_retained_nonfragment_sequences"
		clusters.loc[clusters.number_retained_contigs.eq(0), "selection_reason"]="no_retained_contigs"
		if complete_only.lower() == "true":
		    incomplete=~clusters.representative_CheckV_quality.eq("Complete") & clusters.selected
		    clusters.loc[incomplete, "selected"]=False
		    clusters.loc[incomplete, "selection_reason"]="representative_not_Complete"
		selected=selected.loc[selected.cluster_rep.isin(clusters.loc[clusters.selected, "cluster_rep"])].copy()
		fasta_dir=Path(fasta_dir)
		fasta_dir.mkdir(parents=True, exist_ok=True)
		for cluster, group in selected.loc[selected.is_comparison_record].groupby("cluster_rep", sort=True):
		    SeqIO.write([records[contig] for contig in group.original_id], fasta_dir / (cluster + ".fasta"), "fasta-2line")
		clusters.to_csv(clusters_path, sep="\t", index=False, na_rep="not reported")
		selected.to_csv(members_path, sep="\t", index=False, na_rep="not reported")

		PYTHON
		"""


def input_isolate_viridic_selection(wildcards):
	if not ISOLATE_VIRIDIC:
		return []
	return [checkpoints.select_isolate_viridic_clusters.get(sampling=wildcards.sampling).output.clusters]


def input_isolate_viridic_members(wildcards):
	if not ISOLATE_VIRIDIC:
		return []
	return [checkpoints.select_isolate_viridic_clusters.get(sampling=wildcards.sampling).output.members]


def input_isolate_viridic_results(wildcards):
	if not ISOLATE_VIRIDIC:
		return []
	ckpt=checkpoints.select_isolate_viridic_clusters.get(sampling=wildcards.sampling)
	with open(ckpt.output.clusters) as handle:
		clusters=[row["cluster_rep"] for row in csv.DictReader(handle, delimiter="\t") if row["selected"] == "True"]
	return expand(ISOLATE_VIRIDIC_RESULTS, sampling=wildcards.sampling, cluster=clusters)


# Same VIRIDIC command and output format as the existing relatives rule.
# Skip its stock heatmap: the summary notebook explicitly draws annotated plots.
use rule viridic_relatives_phages as viridic_isolate_cluster with:
	input:
		cat_isolates_relatives=lambda wc: os.path.join(checkpoints.select_isolate_viridic_clusters.get(sampling=wc.sampling).output.fasta_dir, wc.cluster + ".fasta"),
		viridic_singularity_folder=config["viridic_folder"],
	output:
		viridic_out=directory(ISOLATE_VIRIDIC_RESULTS),
	params:
		steps="sim_clust",
		run_log=dirs_dict["ANNOTATION"] + "/VIRIDIC_CLUSTERS/{sampling}/LOGS/{cluster}.log",
	log:
		dirs_dict["ANNOTATION"] + "/VIRIDIC_CLUSTERS/{sampling}/LOGS/{cluster}.log",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/viridic_isolate_cluster/sampling={sampling}/cluster={cluster}.tsv"
	threads: int(config.get("isolate_viridic_threads", 8))
