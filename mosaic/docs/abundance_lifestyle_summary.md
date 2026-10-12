# Final abundance/lifestyle and microdiversity notebooks

`runWorkflow` now includes `10_abundance_lifestyle_summary.tot.html`. The report
uses the existing final `filtered_...tot.fasta`, normalized abundance matrices,
coverage tables and annotations. It does not change catalogue filtering or
normalization. You can also request `runAbundanceSummary` directly.

The source notebook is `notebooks/10_abundance_lifestyle_summary.py.ipynb`.
Snakemake writes the executed notebook to `NOTEBOOKS/` and the HTML report to
`FIGURES_AND_TABLES/`. Tables and PNG/SVG figures are under
`FIGURES_AND_TABLES/10_abundance_lifestyle_summary.tot/`.

## Metadata and interpretation

`vOTU_metadata_polished.csv` uses the DANIEL notebook's capitalization and topic
ordering: viral identification, CheckV, seqtk composition, BACPHLIP, viral taxonomy,
ViralRefSeq/METAVR similarity, resolved hosts, iPHoP and CRISPR evidence, then
abundance. Additional workflow fields (including RNA acceptance routes and
TaxMyPhage taxonomy) precede the final abundance columns. Original vOTU IDs and
all catalogue rows are retained, including undetected and unknown-quality vOTUs.

All workflow samples remain in matrices and prevalence denominators. No control
samples are removed by name. Presence means nonzero **filtered normalized RPKM**.
The existing length-dependent coverage filter is reused, not replaced by a new
75%-only or RPKM cutoff. Relative abundance is each sample's filtered RPKM divided
by its catalogue RPKM sum. It is not a percentage of all sequenced reads.
All-zero libraries have NA relative abundance; mean relative abundance averages
only libraries with a defined denominator. Their counts/presence remain zero.

The report produces sample summaries, category counts and abundance proportions,
genome summaries, prevalence/presence matrices, cumulative abundance, heatmaps,
host/taxonomy panels and sample-accumulation curves. Accumulation uses 100 random
sample orders with seed 42 and percentile bands, not read-depth rarefaction or
confidence intervals. No automatic significance tests are performed.

BACPHLIP scores and VIBRANT lifecycle evidence remain separate. Valid BACPHLIP
scores use the DANIEL 0.5 virulent-score cutoff; missing/invalid scores are Unknown,
not automatically Temperate. Recognized RNA/NCLDV/ssDNA/lavidaviridae evidence is
marked Not applicable for the dsDNA-phage lifestyle summary. Source predictions
remain in the metadata for inspection.

iPHoP candidates require confidence >=90, preferring BLAST-supported candidates
before the highest confidence within that pool, as in DANIEL. Original GTDB names
and lineage ranks are retained without implicit NCBI downloads, suffix stripping,
or project-specific overrides. MetaVR fields describe reference-hit metadata;
they are not promoted to confirmed hosts of the query vOTUs.

## PhaBOX 2 taxonomy evidence

The existing `PhaGCN_family` and `PhaGCN_score` labels are preserved when the
taxonomy rule runs through PhaBOX 2. Full lineage, individual ranks, raw scores,
genus grouping and classifier/database provenance are added as `PhaBOX_` fields
immediately after those existing labels. A family confidence is retained only
when it matches the original scored family; network-propagated taxonomy stays
visible without an invented score. Older PhaGCN CSVs remain readable. See
[PhaBOX taxonomy](phabox_taxonomy.md).

## Optional representative Sourmash evidence

`sourmash_contig_catalogue=True` screens the actual final-catalogue representatives
against GTDB, using the inclusive `sourmash_contig_catalogue_min_length=5000`
cutoff. `Sourmash_` columns are added to `vOTU_metadata_polished.csv`, after
host/CRISPR evidence and before abundance. They retain screening status, query
length/hashes, lineage/rank, taxonomy status, assigned fractions, reference
matches and the database. Short/no-match vOTUs remain in the table. These are
microbial-similarity annotations, not confirmed host assignments or a reason
to remove a vOTU. The same screen annotates selected all-assembled representatives
in their own metadata table. See [catalogue screening](sourmash_contig_catalogue.md).

The read-profile report uses the separate `sourmash_clean_reads` flag;
neither flag enables the other.

## Optional sample groups and CRISPR host metadata

```yaml
abundance_sample_metadata: "/path/to/sample_metadata.csv"
abundance_group_column: "treatment"
abundance_common_prevalence: 0.2
abundance_core_prevalence: 0.5
```

The sample metadata CSV must have unique `sample` IDs matching the read libraries.
Other columns are preserved; the selected grouping column drives group summaries
and an UpSet-style plot of vOTUs detected in every sample of each group. Libraries
without a group remain in the overall analysis. Common/core categories are only
reporting labels and never filter the final catalogue.

When `microbial_spacers` enables SpacePHARER, raw match counts are reported.
To additionally resolve hosts, set `abundance_crispr_host_metadata` to a CSV with
`microbial_contig`, `microbial_contig_length_bp`, and `genus`; optional higher-rank
columns are `domain`, `phylum`, `class`, `order`, `family`. IDs must match spacer
accessions before `_CRISPR_`. Informative matches on microbial contigs >5,000 bp
are ranked by alignment length, adjusted p-value, contig length and accession.
CRISPR takes precedence over iPHoP; genus conflicts are recorded, and missing
CRISPR ranks are not filled from a conflicting iPHoP lineage. No host files from
DANIEL or PhylloVir are implicitly loaded.

## Separate microdiversity report

`notebooks/11_microdiversity_summary.py.ipynb` consumes only the final catalogue
and per-sample `MICRODIVERSITY/.../summary.tsv` files. It has no dependency on
the abundance/lifestyle notebook, annotation tools or RNA enrichment.

`microdiversity: True` schedules its HTML report through the existing automatic
microdiversity branch. Explicit targets `runMicrodiversity` and
`runMicrodiversitySummary` request it even when the automatic flag is off.
No additional reporting flag is needed.

Outputs are `FIGURES_AND_TABLES/11_microdiversity_summary.tot.html`, its matching
table/figure directory and an executed notebook in `NOTEBOOKS/`. Tables preserve
every sample/vOTU pair, analysis thresholds, callability and depth-capped sites.
Pi/entropy are NA when no sites are callable. Sample-level diversity is weighted
by callable sites, including invariant sites; mean depth is weighted by length.
Heatmaps show pi alongside callable fraction. These are observed read-level
estimates, not error-corrected population diversity, validated variants or
evidence of viral activity.

## Checks

`python -m unittest discover -s mosaic/tests -p test_summary_notebooks.py -v`
executes the actual notebook cells on synthetic data, including optional evidence,
column ordering, zero-detection samples, empty catalogues, conflicting hosts,
callable-site weighting and missing/mismatched inputs. Run in the workflow
environment with matplotlib, seaborn, pandas, NumPy, Biopython and IPython.
