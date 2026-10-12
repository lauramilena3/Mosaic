# RefSeq Viral read mapping

Enable in `config.yaml`:

```yaml
map_to_RefSeq: True
RefSeqViral_db: "db/RefSeqViral/RefSeq_viral.fasta"
```

The flag defaults to `False`. When enabled, `all`, `runWorkflow`,
`runAbundance`, and `phage_isolates` request RefSeq mapping and the detection
report. It is independent of `RNA_enriched`,
`microdiversity`, MetaVR BLAST and `additional_reference_contigs`. It does not
replace or change the shared filtered vOTU catalogue or its abundance outputs.

## Dated RefSeq download

The default `RefSeqViral_db: "db/RefSeqViral/RefSeq_viral.fasta"` schedules
`downloadRefSeqViral` automatically if the FASTA is missing. Override this
setting with an existing local FASTA to skip the download. The rule downloads
NCBI's viral RefSeq genomic FASTA into a
UTC-dated directory under `db/RefSeqViral/`, builds its nucleotide BLAST
index there, and creates undated symlinks for the FASTA and index files.
The dated directory also contains the original archive, HTTP response headers,
and `download_info.tsv` with the source URL, UTC timestamps and SHA-256 hashes.
`db/RefSeqViral/RefSeq_viral.download_info.tsv` links to the active manifest.
The download date records when this copy was retrieved, not the NCBI release date.
The source URL can be overridden with `refseq_viral_download_url`.

## Standalone execution

From the workflow directory with your usual configuration and mapping environments:

```bash
snakemake --use-conda --cores 16 map_to_RefSeq
```

The `map_to_RefSeq` target requires `map_to_RefSeq: True`. It requests the
cleaned reads and post-QC read-length table, building missing QC prerequisites,
but does not require assembly or vOTU discovery.

The wrapper also supports `mosaic run refseq_mapping --raw /path/00_RAW_DATA
--results /path/results`, or `--map-to-refseq` / `--no-map-to-refseq` with an
end-to-end or abundance run. A flag does not expand unrelated targets such as
QC-only or assembly-only modes. Set `RefSeqViral_db` in your configuration, or
pass `--config RefSeqViral_db=/path/to/RefSeq_viral.fasta` to the wrapper.

## Mapping and normalization

Every discovered short-read sample uses its full cleaned paired-end `tot` files.
In virome mode this branch does not map raw/unpaired/long reads or create a `sub` RefSeq analysis,
even when other branches use subsampling. Isolate mode additionally maps every
QC-passed orphan read and includes primary single-mate/discordant alignments in
the no-`XS:i:` unique subset. RNA and DNA libraries use the same
mapping procedure. The configured database must be a nucleotide FASTA.

The rules in `07_abundance.smk` follow `mapReads_reference`: Bowtie2 `--fast`,
CoverM filtering at 95% read identity and 85% aligned read length, alignment
flag statistics, and coverage statistics for filtered alignments and the existing
proper-pair/no-`XS:i:` subset. The latter is the workflow's operational
"unique"-mapping heuristic, not proof that a read has only one possible match.

The six-file Bowtie2 large index is stored under the results directory. The
shared source database is not modified and `additional_reference_contigs` is
not repointed. This branch uses the existing `env1_mapping.yaml` environment.

`normalise_reads_RefSeq` reuses `07_Normalise.py.ipynb`, including its existing
ambiguity adjustment and 200 uniquely covered-base cutoff. Normalized RPKM uses
the notebook's adjusted mapped-read denominator, not total sequenced library
reads. `raw` in output filenames means before this normalization, **not raw
input reads**. Counts are read-alignment counts, not deduplicated molecule counts.

The existing abundance-table filters are also preserved: >5,000 covered bases
for references of length >=6,667 nt, or >75% breadth for shorter references;
the separate `filtered_75` tables use >75% breadth for all lengths. Unfiltered
tables remain available. These reporting filters do not remove FASTA records.

## Outputs

Under `results_dir/06_MAPPING/REFSEQ_VIRAL/`:

- `RefSeqViral_RPKM_raw_tot.txt` and `RefSeqViral_RPKM_normalised_tot.txt`.
- `RefSeqViral_counts_raw_tot.txt` and `RefSeqViral_counts_normalised_tot.txt`.
- `RefSeqViral_breadth_coverage_percent_tot.txt` and
  `RefSeqViral_breadth_coverage_bases_tot.txt`.
- `RefSeqViral_mean_depth_tot.txt`: CoverM mean depth, not count divided by length.
- Corresponding `filtered_` count/RPKM and `filtered_75_` RPKM tables.
- Per-sample BAMs, Bowtie2 logs, flag statistics, full/unique coverage statistics,
  and optional nonzero-position base coverage files. Intermediate BAMs are temporary.

`per_base_coverage` defaults to `False`. Enable it with
`--config per_base_coverage=True` to export
`bowtie2_RefSeqViral_<sample>_tot_basecov.txt.gz` and
`bowtie2_RefSeqViral_<sample>_tot_unique_basecov.txt.gz`. These are gzipped,
tab-separated, zero-based positions with positive depth. CoverM summaries,
normalisation and detection tables are produced with either setting. The same
flag controls the optional exports for the general `map_to_reference` target.

Merged `.txt` matrices are comma-separated, matching the existing normalization
outputs. Rows are labeled `RefSeq_accession`, columns are configured samples;
zero-hit samples remain as zero columns. Per-sample coverage-statistics files
are tab-separated. No taxonomy/species aggregation or segmented-genome grouping
is performed: each FASTA sequence is reported separately.

The executed notebook is `results_dir/NOTEBOOKS/07_Normalise_RefSeqViral.tot.ipynb`.
The shared notebook now uses declared input files/sample order and handles
zero-hit samples and small plotting layouts without changing its positive-count
normalization formula. Reference matches should not be treated as independent,
confirmed virus detections; inspect breadth and ambiguity alongside counts.

## Detection report

Set `negative_control` to an exact paired-end sample ID in the run, or leave it
empty. No metadata table is accepted. The notebook downloads accession records
from NCBI Virus, molecule/strand information from NCBI Nucleotide, and ranked
lineage from NCBI Taxonomy in batches using Python. Set `NCBI_API_KEY` and
optionally `NCBI_EMAIL` in the environment if available; neither is required.
Network access is required on the first run. Downloads are resumed from
`06_MAPPING/REFSEQ_VIRAL/RefSeqViral_NCBI_metadata_cache.sqlite`; set
`refseq_metadata_refresh: True` to re-download. The RefSeq FASTA supplies local
descriptions and fallback names. Missing NCBI annotations remain blank rather
than being guessed. Genome-type plots are made only if RNA/DNA data are returned.

`07_RefSeq_detection.py.ipynb` reads the two per-sample CoverM files, optional
gzipped base-coverage files and the existing raw RPKM, count, breadth and
mean-depth matrices. Without per-base exports, positional diagnostic columns
are empty rather than zero; detection calls still use the CoverM summaries. It
creates one row for every accession/sample, including zero hits. A high-confidence
call passes the existing >5,000 covered bases for genomes >=6,667 bp or >75%
breadth for shorter genomes, plus at least five unique mappings. Candidates need
at least ten unique mappings and 25% breadth. Signals with at least five mapped
reads and a unique/all ratio <=0.2 are ambiguous. With a negative control,
signal with >=5 control reads and <=2-fold sample/control RPKM enrichment is
also ambiguous. The control library's high-confidence/candidate evidence is
reported separately as potential contamination. All
thresholds and the enrichment pseudocount have named `refseq_*` config keys.

The five `RefSeqViral_detection_*` / `RefSeqViral_detected_*` TSV outputs
are in `06_MAPPING/REFSEQ_VIRAL`. PNG and SVG figures are in
`FIGURES_AND_TABLES/07_RefSeq_detection`; the executed notebook is in
`NOTEBOOKS/07_RefSeq_detection.py.ipynb`. Scatterplots include mapped-read
evidence; heatmaps show high-confidence/candidate accessions, limited to the
configured top accessions. The TSVs are not truncated.

## Tests

```bash
MPLBACKEND=Agg MOSAIC_SNAKEMAKE=/path/to/snakemake \
MOSAIC_MAPPING_BIN=/path/to/mapping/environment/bin \
    python -m unittest discover -s mosaic/tests -p 'test_refseq*.py' -v
```

Tests exercise automatic flag gating, standalone DAG/CLI wiring, known
normalization values, zero-hit samples, existing-reference compatibility,
detection classes/control handling, metadata/plots, and actual mapping against
a tiny synthetic FASTA. Omitting the environment
variables skips the corresponding external-tool tests.
