# Sourmash screening of catalogue representatives

For viromes, add the screen to `runWorkflow`:

```bash
snakemake --use-conda -p runWorkflow --config \
  input_dir=/path/to/00_RAW_DATA \
  sourmash_contig_catalogue=True sourmash_contig_catalogue_min_length=5000 \
  map_to_all_assembled=True \
  -j 32 --rerun-incomplete -n
```

Remove `-n` to execute. The final filtered vOTUs are screened; with
`map_to_all_assembled=True`, the abundance-selected cluster representatives
are screened too. Omit the mapping flag if that branch is not needed.

From `mosaic/`, add these options to the usual `phage_isolates` command:

```bash
snakemake --use-conda -p phage_isolates --config \
  input_dir=/path/to/00_RAW_DATA isolates=True metagenome=False \
  sourmash_contig_catalogue=True sourmash_contig_catalogue_min_length=5000 \
  -j 32 --rerun-incomplete
```

This flag is independent of `sourmash_clean_reads`, the final-clean-read report.
It defaults to `False`. The cutoff is inclusive and applies only to Sourmash
queries, not to assembly saving, viral acceptance, abundance ranking or
retained/excluded decisions. The Python wrapper also accepts
`--sourmash-contig-catalogue --sourmash-contig-catalogue-min-length 5000`.

## Catalogues and metadata

| Context | FASTA screened | Metadata receiving the evidence |
| --- | --- | --- |
| Virome | `05_vOTUs/filtered_95-85_positive_viral_contigs.tot.fasta` | `FIGURES_AND_TABLES/10_abundance_lifestyle_summary.tot/vOTU_metadata_polished.csv` |
| Selected all-assembled branch | `03_CONTIGS/ALL_ASSEMBLED/all_assembled_top_contigs_tot.fasta` | `03_CONTIGS/ALL_ASSEMBLED/all_assembled_top_contigs_metadata_tot.tsv` |
| Phage isolates | `03_CONTIGS/ALL_ASSEMBLED/phage_isolates_cluster_representatives.tot.fasta` | Existing all-contig and cluster metadata tables |

The top-contig FASTA contains all selected cluster representatives, not a fixed
100 sequences. Selection and the existing independent 95/85 clustering are
unchanged. No new clustering or remapping is introduced. The existing combined
microbial catalogue keeps its original unfiltered Sourmash behavior.

## Reused rules

The global isolate catalogue includes every saved SPAdes contig >=1 kb,
including both retained and excluded contigs. Exact MMseqs dereplication is followed by the original CheckV
independent >=95% ANI and >=85% target-coverage clustering criterion.

`vOUTclustering_get_new_references` also extracts the global cluster centroids
directly from the existing `.clstr` file. It writes
`03_CONTIGS/ALL_ASSEMBLED/phage_isolates_cluster_representatives.tot.fasta` and
its lengths table. This is not the viral-only representative selection: host,
plasmid, unresolved and viral clusters are all retained.

The existing `single_fasta_microbial`, `sourmash_sketch_microbial`,
`sourmash_gather_microbial` and `sourmash_tax_microbial` rules now accept this
catalogue, filtered vOTUs, selected all-assembled representatives and the existing
combined microbial catalogue. They use the
same environments, k=31/scaled=1000 sketches and configured `sourmash_rocksdb`
and `sourmash_tax` files. Their existing download rules provide GTDB RS226 by
default. No new database downloader or environment is required. The shared
`sourmash_contig_metadata` rule reuses the isolate notebook's evidence parser
for the two virome catalogues. Its executed parsing notebook is a benchmark/log
artifact, not another user-facing analysis notebook.

Only representatives whose actual FASTA length is **>=5,000 bp** enter this
screen, unless the configurable cutoff is changed. The cutoff does not remove
shorter sequences from the representative FASTA or either metadata table.

## Results and interpretation

The polished virome table uses `Sourmash_` labels in its existing metadata topic
order, before abundance. Selected-contig and isolate tables keep `sourmash_`
labels. Fields are joined by representative ID; shorter/no-match rows remain
visible with explicit screening statuses. Missing evidence is `not reported`,
not an invented zero.

The catalogue notebook adds `sourmash_` columns to both existing isolate outputs; the
isolate-summary notebook displays that evidence alongside retention decisions:

- `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/cluster_metadata.tsv`
- `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/all_contig_metadata.tsv`

Columns include representative length, sketch hashes, screening status, GTDB
lineage and rank, taxonomy status and assigned fraction, best reference match,
total matched-query fraction, shared hashes, estimated shared bases, number of
reference matches and the database path. `sourmash_taxonomy_estimated_ani` is
Sourmash's sketch-based estimate, not the BLAST ANI used for MOSAIC clustering.

The evidence scope is always `cluster_representative`. Every original member
in the isolate catalogue inherits its representative's annotation; members were not searched
separately. Their own geNomad, CheckV and host BLAST evidence stays unchanged.

Screening statuses distinguish `below_length_cutoff`, `no_usable_hashes`,
`no_reported_match` and `screened` (a gather match was reported). The separate
taxonomy status comes directly from Sourmash; `sourmash_gtdb_match=True` means
taxonomy status `match`, not confirmed host contamination. Short/no-match
clusters remain in the catalogue and are not labelled host-free.

The existing gather threshold of 0 and Sourmash's default taxonomy thresholds
are retained; taxonomy is not forced to strain rank. A gather match can have
very few supporting hashes, especially on short fragments. Review the fractions
and hash support alongside host BLAST and geNomad before concluding host origin:
viral or plasmid sequences can also share sequence with microbial genomes.
The screen does not change retained/excluded decisions or discard contigs.
No six-category biological classification is imposed in the isolate workflow.

Raw taxonomy and per-query sketch measurements remain in `07_ANNOTATION/`:

- `sourmash_phage_isolates_cluster_representatives_tot.classifications.csv`
- `phage_isolates_cluster_representatives_tot_sourmash_queries.tsv`

Virome outputs use catalogue names `filtered_votus_tot` and
`all_assembled_top_contigs_tot` with the same naming convention. Collected
evidence is saved as `sourmash_<catalogue>.metadata.tsv`, including the screening
cutoff, and joined into final metadata. Catalogue IDs remain independent;
evidence is not silently transferred between main and selected catalogues.

Single-contig split FASTAs, sketch ZIPs and gather intermediates retain the
existing temporary-file behavior. Empty selections and no-match searches
produce header-only results. With the flag disabled, the existing metadata
columns and other workflow modes remain unchanged.

See the [isolate workflow guide](isolate_purity.md) for inputs, read accounting,
the output index and the step-by-step workflow diagrams.
