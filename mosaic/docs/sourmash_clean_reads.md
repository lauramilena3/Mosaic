# Clean-read Sourmash profiles

Enable the optional report in a virome run:

```bash
snakemake --use-conda -p runWorkflow \
  --config input_dir=/path/to/00_RAW_DATA sourmash_clean_reads=True \
  -j 32 --rerun-incomplete -n
```

Remove `-n` to execute. For isolates, use:

```bash
snakemake --use-conda -p phage_isolates \
  --config input_dir=/path/to/00_RAW_DATA isolates=True metagenome=False sourmash_clean_reads=True \
  -j 32 --rerun-incomplete -n
```

The corresponding defaults in `config.yaml` are:

```yaml
sourmash_clean_reads: False
sourmash_clean_min_shared_bp: 50000
sourmash_clean_min_reference_fraction: 0.10
```

`sourmash_clean_reads` is the single read-profiling flag. It replaces the older
trimmed-pair `sourmash` setting and its `--sourmash` wrapper option. Use
`--sourmash-clean-reads` with the Python wrapper. QC, `runWorkflow`,
`phage_isolates` and microbial targets that request QC include the shared report
when enabled. The `microbial_metagenome` wrapper enables it by default; other
presets default to off. The assembled-cluster `sourmash_contig_catalogue` screen
remains separate and unchanged.

To request only this report and any missing cleaned-read prerequisites:

```bash
snakemake --use-conda -p runSourmashCleanReads \
  --config input_dir=/path/to/00_RAW_DATA sourmash_clean_reads=True \
  -j 32 --rerun-incomplete -n
```

## Reads and pooling

Every sample contributes all three final-clean files:

```text
02_CLEAN_DATA/{sample}_forward_paired_clean.tot.fastq.gz
02_CLEAN_DATA/{sample}_reverse_paired_clean.tot.fastq.gz
02_CLEAN_DATA/{sample}_unpaired_clean.tot.fastq.gz
```

These are full libraries, not the 2M subset or BBNorm-normalised reads.
Abundance-tracking DNA sketches use `k=31`, `scaled=1000`.
`sourmash_pool_clean_reads` merges every sample sketch, including any negative
control. No additional pooled FASTQ is created. The pool is weighted by actual
sequencing contribution, not an equal-weight average of sample profiles.

The existing `sourmash_sketch_clean_reads`, `sourmash_pool_clean_reads`,
`sourmash_gather` and `sourmash_tax` rules, environment, GTDB RS226 RocksDB
and matching taxonomy downloaders are reused. All read profiles use standard
`sourmash gather`, retaining abundance information. The duplicated trimmed-pair
sketch and `fastmultigather` branch are no longer scheduled. Each profiling
rule and the shared report have their own benchmark. Old results are not deleted.

## Evidence and database scope

Gather returns reference genomes with at least the configured shared-base
overlap. A returned match is labelled `supported` when both apply:

- `unique_intersect_bp >= sourmash_clean_min_shared_bp`;
- `f_match_orig >= sourmash_clean_min_reference_fraction`.

Other returned matches are `weak`, not discarded. Matches below the gather
threshold are not searched/reported. This is supporting evidence, not a
confirmed contamination or species-identification call.

Reference evidence includes accession/name, lineage, family/genus/species,
estimated shared bases/hashes, reference k-mer fraction/percent,
abundance-weighted query fraction and mean/median matching k-mer abundance.
Shared bases are estimated from sampled hashes. Reference k-mer containment
is not read-mapping breadth, and k-mer abundance is not sequencing depth.
Weighted query fractions are not exact percentages of sequencing reads.

The configured GTDB database screens bacteria and archaea, not all possible
contaminants. Unmatched signal may come from phages, eukaryotes, PhiX or taxa
absent from the database. In isolates, biological filtering is bypassed, so
these are all fastp-passed pairs and orphans. In viromes, the report describes
what remains after the configured eukaryotic/contaminant filtering.

No reads are removed and no contig decision is changed. A specified negative
control has a red sample label, but is not subtracted or tested for enrichment.
It remains in both the individual profiles and the pool.

## Outputs

Raw sketches, gather results and taxonomy outputs stay under
`07_ANNOTATION/SOURMASH_CLEAN/`:

```text
SAMPLES/{sample}_sourmash.sig.zip
SAMPLES/{sample}_gather_sourmash.csv
SAMPLES/{sample}_gather_sourmash.with-lineages.csv
SAMPLES/{sample}_sourmash.summarized.csv
SAMPLES/{sample}_sourmash.kreport.txt
POOLED/pooled_sourmash.sig.zip
POOLED/pooled_gather_sourmash.csv
POOLED/pooled_gather_sourmash.with-lineages.csv
POOLED/pooled_sourmash.summarized.csv
POOLED/pooled_sourmash.kreport.txt
```

The executed notebook is `NOTEBOOKS/01_sourmash_clean_reads.tot.ipynb`.
Its HTML, tables and figures are under `FIGURES_AND_TABLES/`:

```text
01_sourmash_clean_reads.tot.html
01_sourmash_clean_reads.tot.tsv
01_sourmash_clean_taxonomy.tot.tsv
01_sourmash_clean_profile_summary.tot.tsv
01_sourmash_clean_genus.tot.tsv
01_sourmash_clean_pooled.tot.png / .svg
01_sourmash_clean_genus.tot.png / .svg
```

The summary retains every sample and the pool: returned/supported/weak
reference counts, distinct assigned species counts and matched weighted
query fraction. Species counts are database labels, not confirmed species.
The genus matrix contains full untransformed weighted query percentages.
All evidence stays in the TSV files; notebook previews show the first 30 rows.

The notebook has explicit parsing and plotting code, a pooled genus barplot
and per-sample genus heatmap. Plots show up to 30 genera and at most 50 samples
per row, keeping the same sample-cell width in the final row. A heatmap star
means at least one reference assigned to that genus passed both support
filters. Blank stars do not mean absence.

The isolate notebook/HTML reads the shared tables and reuses the figures,
without repeating profiling, parsing or plotting. It keeps the existing copies:

```text
08_phage_isolates_sourmash_clean_reads.tot.tsv
08_phage_isolates_sourmash_clean_taxonomy.tot.tsv
08_phage_isolates_sourmash_clean_pooled.tot.png / .svg
08_phage_isolates_sourmash_clean_genus.tot.png / .svg
```

Zero-hit samples remain visible. Empty/no-match profiles produce header-only
evidence tables, zero-match summary rows and explanatory placeholder plots,
not a claim that the sample is contamination-free.

See the [isolate guide](isolate_purity.md) and
[viral metagenome manual](viral_metagenome.md) for the wider workflows.
