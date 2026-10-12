# Cenote-Taker3 on abundance-selected cluster representatives

Run from the `mosaic/` directory with your usual configuration, adding:

```bash
snakemake --use-conda -p runWorkflow --config \
  input_dir=/path/to/00_RAW_DATA \
  map_to_all_assembled=True run_cenote=True \
  all_assembled_cluster_top_n=1000 all_assembled_top_per_sample=100 \
  -j 32 --rerun-incomplete
```

The automatic Cenote branch requires both flags. `run_cenote=False` is the
default and leaves the existing top-contig analysis unchanged. RNA enrichment
is independent: when enabled, RNA assemblies can contribute to the same
abundance-ranked selection. These are selected cluster representatives, not
necessarily complete or viral genomes.

## Selection, clustering and mapping

`map_to_all_assembled=True` maps the same existing subset of up to 2 million
cleaned read pairs per sample once, to the full MMseqs-dereplicated catalogue.
Selection takes the union of the top `all_assembled_top_per_sample` contigs
per non-NC sample (100 by default) and the top `all_assembled_cluster_top_n`
contigs by mean raw RPKM across non-NC samples (1,000 by default). Only non-zero
abundances qualify, and shared contigs are counted once. There is no overall
selection cap: the union can contain more than 1,000 contigs. Ties are resolved
by contig ID. NC values remain in the reported columns but do not affect selection.

The selection is written to
`03_CONTIGS/ALL_ASSEMBLED/all_assembled_selected_contigs.tot.fasta`, with ranking
and MMseqs membership tables alongside it. The existing `vOUTclustering` rule
compares only this selection, using `anicalc_checkv.py` and CheckV's original `aniclust`
with `--min_ani 95 --min_tcov 85 --min_qcov 0`: independently at least 95% ANI
and 85% aligned coverage of the shorter sequence. Results are in
`all_assembled_selected_contigs.tot_95-85.clstr` and its BLAST/ANI tables.

Every longest-first cluster representative is kept for annotation. There is
no final top-100 selection, no additional mapping, and no summing of member
RPKMs. Existing `all_assembled_top_contigs` output names are retained, but now
contain all representatives. Their abundance columns remain the representative
contig's measurements from the original full-catalogue mapping, as recorded in
`abundance_scope`. `selected_contigs_in_cluster` counts the selected dereplicated
sequences; `cluster_size` includes their known MMseqs members. The separate final
membership table links each representative to these renamed original members
and their dereplicated representatives. Unselected relatives are not assigned
to these ANI clusters. The main filtered viral vOTU catalogue is unchanged.

geNomad evidence is collected by the current contig IDs from the existing
per-sample DNA results in `04_VIRAL_ID/{sample}_geNomad_tot` and, when RNA
enrichment is enabled, each `03_CONTIGS/RNA/{sample}/{assembler}_genomad` result.
No combined `all_assembled_geNomad_tot` run is requested. The selected-set
metadata retains virus/plasmid calls, scores, taxonomy/topology, provirus-region
and hallmark evidence, including `not reported` where evidence is absent.
Previously generated combined geNomad results are not used or deleted.

With `sourmash_contig_catalogue=True`, selected representatives at least
`sourmash_contig_catalogue_min_length` (default 5,000 bp) also receive the existing
GTDB Sourmash screen. `sourmash_` evidence is joined by sequence ID into
`03_CONTIGS/ALL_ASSEMBLED/all_assembled_top_contigs_metadata_tot.tsv`.
Shorter/no-match rows remain visible; ranking, cluster membership and abundance
values are unchanged. See [representative screening](sourmash_contig_catalogue.md).

### Annotation reuse audit

| Analysis | Current decision | Reason |
| --- | --- | --- |
| Combined full-assembly geNomad | Removed from the workflow dependencies | DNA and individual RNA assembly results already provide the needed evidence. |
| Isolate and all-assembled concatenation/MMseqs dereplication | Shared when RNA enrichment is off | One `all_assembled` DNA catalogue is generated. Isolate FASTA/provenance and representative/membership filenames are symlinks to those outputs. With RNA enrichment on, both catalogues are generated separately. |
| Selected-cluster VIBRANT and VirSorter2 | Keep | The main viral-catalogue analyses do not cover every selected original assembly contig; the RNA candidate run also uses a different group selection. |
| Selected-cluster CheckV | Keep | The main CheckV inputs are viral candidates, sometimes extracted provirus regions, rather than every selected whole contig. The isolate catalogue already reuses per-sample isolate CheckV outputs. |
| Selected-cluster RefSeq/METAVR BLAST and Pharokka | Keep | The existing main-catalogue outputs do not cover every selected representative. |
| Selected-cluster Cenote | Optional, unchanged | Requested only with both `map_to_all_assembled=True` and `run_cenote=True`; it adds separate discovery and annotation evidence. |
| Main-catalogue `genomad_vOTUs` | Keep, possible simplification to evaluate separately | Its conservative, default and relaxed outputs are all consumed by vOTU filtering. Removing a preset would change that filtering. |

Overlapping contig IDs alone do not establish reusable results: the sequence
may have been extracted or trimmed, and tool parameters may differ. Reusing
additional annotations would require checking identical sequences and
parameters, then running each tool only on missing sequences. That broader
change is not part of removing the redundant combined geNomad run.

Repeated occurrences of the same input path do not cause repeated jobs:
Snakemake schedules one producer job for that output. With RNA enrichment off,
`reuse_dna_contigs_for_isolates` and `reuse_dna_dereplication_for_isolates`
preserve the isolate filenames without repeating concatenation or MMseqs.
The symlink targets are relative filenames in the same directory, so moving
the project directory does not break these links.
Requesting isolate outputs alone also uses this shared DNA catalogue, but does
not request all-assembled mapping or selected-set annotation. Subsequent ANI
clustering stays separate: the isolate catalogue clusters all retained DNA
contigs, while the abundance branch clusters its selected subset. Existing
study outputs are not migrated or deleted automatically by this code change.

Full-catalogue files use the `all_assembled` prefix, including
`all_assembled_RPKM_raw_tot.txt` and `all_assembled_mapping_summary_tot.tsv`.
Full-catalogue normalization outputs are unchanged. This branch no longer
creates an `all_assembled_top` mapping index, BAMs, coverage tables or mapping
summary. Old remapping outputs from earlier runs are not used or automatically
deleted.

The flag also requests `NOTEBOOKS/07_mapping_statistics_tot.ipynb`, even without
`assembly_stats=True`. Its HTML table includes the full-catalogue mapping
percentage, and a per-sample barplot is displayed in the notebook and saved as
PNG/SVG in `FIGURES_AND_TABLES/07_mapping_statistics_all_assembled_tot.*`.
There is no selected-set mapping percentage because the selection is not
remapped. The report retains its existing mapping comparisons, so their mapping
inputs are also requested if not already available.

An empty abundance selection produces empty clustered FASTAs and header-only
selection/metadata tables without running BLAST. This branch does not require
`run_cenote=True`; that flag only enables the additional Cenote annotation.

## Rules and environment

`downloadCenoteDB` downloads the core HMM, hallmark-taxonomy, RefSeq-taxonomy,
MMseqs CDD and viral-domain databases into `cenote_db`, by default
`db/cenote-taker3`. This directory is ignored by Git. `download_date.txt` records
the successful installation day in UTC. The download is approximately 3 GB
decompressed; optional HHsuite databases are not downloaded.

`annotate_cenote` runs discovery plus annotation on
`03_CONTIGS/ALL_ASSEMBLED/all_assembled_top_contigs_tot.fasta`. It uses virion and
RdRP hallmark evidence, not forced annotation-only mode. The generic rule also
uses the existing `annotation_fastas` mechanism, but the workflow automatically
requests it only for the top selection.

Both rules use `envs/cenote.yaml`, pinning Cenote-Taker3 3.4.4 and Python 3.11.
The current Pharokka environment cannot share this version without changing its
existing Mash/GSL pins. The separate environment leaves those pins untouched.

Configuration:

```yaml
run_cenote: False
cenote_min_contig_length: 1000
cenote_threads: 8
cenote_mem_mb: 16000
cenote_prune_prophage: False
cenote_db: "db/cenote-taker3"
```

The memory value is a Snakemake scheduling reservation, not a hard memory limit
inside Cenote. Its circularity check requires at least 1,000 nt; the linear
length cutoff follows `cenote_min_contig_length`.

## Results

Results are in `07_ANNOTATION/Cenote_all_assembled_top_contigs.tot/`, including:

- `mosaic_ct3_virus_summary.tsv`: reported viruses, hallmark counts, taxonomy and
  annotation information, including the original identifier in `input_name`.
- `final_genes_to_contigs_annotation_summary.tsv`: gene annotations using Cenote's
  internal identifiers.
- `contig_name_map.tsv`: internal identifiers mapped to input FASTA headers.
- `sequin_and_genome_maps/`: GenBank files and related annotation outputs, when
  produced by Cenote.
- `mosaic_ct3_prune_summary.tsv`: region coordinates, when pruning produces regions.
- Run arguments, filtered sequences and the tool's own log, when produced.

The outer log is `07_ANNOTATION/Cenote_all_assembled_top_contigs.tot.log`.
Both rules also have their own benchmarks.

`03_CONTIGS/ALL_ASSEMBLED/all_assembled_top_contigs_metadata_tot.tsv` gains
`Cenote_` columns when this branch is enabled. Results are joined using
`input_name`, not Cenote's renamed contig IDs. Every selected contig remains in
the table. Short contigs are marked `below length cutoff` in `Cenote_assessed`;
unreported classifications and annotations say `not reported`, not “non-viral”.
Multiple reported regions are retained as ` | `-separated values, not silently
reduced to the first hit.

Enabling Cenote does not change the assembly catalogue, top FASTA, raw RPKM
ranking, mapping or normalization. Rotation is disabled, but Cenote can still trim terminal
repeats in its own output sequences. Metadata distinguishes whole-contig,
processed-contig and region evidence. Do not transfer gene coordinates from a
processed sequence directly to the original assembly without accounting for
that change.

Cenote's sequence-submission exports use its default `DNA` molecule-type label.
This is a supplied export label, not an inferred genome type. In mixed RNA/DNA
analyses, inspect that field before using the GenBank files for submission or
genome-type summaries; RNA assembly origin alone does not establish an RNA
virus genome type. The hallmark/taxonomy results are collected independently.

Temporary working files are cleaned on success and failure. The identifier map
is retained, while `ct_processing` is removed. An empty selection, a selection
below the length cutoff, or a completed search without qualifying hallmark hits
produces header-only annotation tables. Missing outputs from an unsuccessful
tool invocation are not silently treated as a negative result.

Upstream installation, modes and output documentation:
[Cenote-Taker3](https://github.com/mtisza1/Cenote-Taker3).
