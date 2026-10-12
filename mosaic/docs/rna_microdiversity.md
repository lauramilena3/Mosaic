# Optional RNA assembly and microdiversity

These are independent, opt-in flags in `mosaic/config.yaml`:

```yaml
RNA_enriched: False
microdiversity: True
```

| RNA_enriched | microdiversity | Additional work in `all` / `runWorkflow` |
| --- | --- | --- |
| False | False | None; existing workflow |
| True | False | RNA assembly and candidates added to the existing vOTU clustering/filtering |
| False | True | Per-sample variation against the existing final filtered vOTU catalogue |
| True | True | RNA-enriched shared catalogue, then per-sample variation against it |

Both branches use every paired-end sample discovered in `input_dir` (normally
`00_RAW_DATA`), following the existing filename/tag conventions. There is no
`rna_samples` list. They use full **cleaned paired reads**, before subsampling or
BBnorm, not raw reads or orphan reads. Neither branch currently profiles long reads.

## Running

From the repository root, with input and output paths adjusted:

```bash
# Existing assembly plus microdiversity, without the RNA branch:
python mosaic/mosaic.py run end_to_end --raw /data/00_RAW_DATA \
    --results /data/results --no-rna-enriched --microdiversity

# Enable both for every paired-end library:
python mosaic/mosaic.py run end_to_end --raw /data/00_RAW_DATA \
    --results /data/results_rna --rna-enriched --microdiversity

# Just microdiversity and its missing prerequisites:
python mosaic/mosaic.py run microdiversity --raw /data/00_RAW_DATA \
    --results /data/results --no-rna-enriched
```

The wrapper enables Conda by default. Add `--dry-run` to inspect scheduled jobs.
Direct Snakemake users can set `--config RNA_enriched=True microdiversity=True`
when requesting `runWorkflow`, or use the dedicated targets `runRNAAssembly`,
`runRNA` and `runMicrodiversity`. Explicitly requesting `runMicrodiversity` runs it
even when the automatic `microdiversity` flag is off. `runQC` keeps its existing
scope. With `RNA_enriched=True`, `runAssembly` also requests the RNA assemblies
and their assembly report. Any target that needs to build the final catalogue
also includes its RNA prerequisites when `RNA_enriched=True`.

There is only one final catalogue path, with or without RNA enrichment:
`05_vOTUs/filtered_95-85_positive_viral_contigs.tot.fasta` by default
(`filtered_<REPRESENTATIVE_CONTIGS_BASE>.tot.fasta` for customized catalogue names).
Changing the RNA flag can rebuild that file and its downstream results. Use separate
results directories if you want to preserve both versions for comparison.
The former `microdiversity_reference` override is no longer supported.

## RNA branch

`rna_assemblers` defaults to `"rnaviralspades megahit trinity"`. All three run
independently on each sample. SPAdes uses **`--rnaviral`**, not ordinary `--rna`.
These assemblies are additional to the standard MOSAIC assembly.
Each assembler has its own explicit rule in `rules/12_rna_microdiversity.smk`,
with its own command, log and benchmark. DNA SPAdes, RNAviralSPAdes, hybrid SPAdes
and SPAdes assembly-depth tests share `env3.yaml` (SPAdes 4.3.0 and seqtk 1.5).
MEGAHIT and Trinity retain their dedicated environments. Assembly reconciliation
and microdiversity calculations are also written directly in their rule blocks.

The companion `03_assembly_short_RNA.py.ipynb` runs after the selected assemblies
are combined for all samples. Its executed notebook is saved in `NOTEBOOKS/`;
the summary and provenance CSVs, plus four PNG/SVG plots, are saved in
`FIGURES_AND_TABLES/` with the `03_assembly_short_RNA_` prefix. It reports
per-assembler and combined contig counts, assembled bases, length distributions,
N50, and the counts retained or collapsed as exact duplicates. All saved assembler
FASTAs already meet `rna_min_contig_length`; the legacy below-minimum provenance
column is normally zero. Empty assemblies remain visible. The `combined_derreplicated` FASTA is
not a vOTU catalogue and has not yet been screened by VirSorter2.

With `RNA_enriched=True`, `assemblyStats_RNA` also runs reference-free QUAST on
each selected RNA assembler and each sample's combined FASTA, reusing `env3.yaml`.
The report is in `03_CONTIGS/RNA/assembly_quast_report.tot.txt`, and the full
QUAST output, including `transposed_report.tsv`, is in
`03_CONTIGS/RNA/statistics_quast_tot/`. Unique `sample_assembler` labels distinguish
the sample-prefixed assembly filenames. All saved contigs are
included (`--min-contig 1`); empty assemblies are skipped and remain visible in
the notebook's existing FASTA summary. If all RNA assemblies are empty, a
header-only QUAST table is produced. The QUAST table is displayed in the RNA
assembly notebook. The DNA QUAST report and notebook are unchanged.

Every saved assembler FASTA is filtered at `rna_min_contig_length` (default 500 nt,
inclusive) before assigning contig identifiers, running geNomad, testing terminal
repeats or combining assemblies. Sequences are then pooled per sample, and exact
duplicates/reverse complements are collapsed. VirSorter2
selects RNA-virus candidates using `rna_viral_groups: "RNA"`. Identification
is a candidate screen, not confirmation that every retained sequence is viral.

RNA candidates enter the existing `tot` catalogue workflow **before clustering**,
alongside the standard viral contigs. The existing MMseqs exact dereplication,
BLAST, `anicalc_checkv.py`, CheckV's original `aniclust`, and representative selection are
reused. Clustering arguments stay `--min_ani 95 --min_tcov 85 --min_qcov 0`.
ANI must independently reach 95% and aligned coverage of the shorter sequence
must reach 85%. The modified product-based helper is not used. There is no CD-HIT catalogue
branch and no additional RNA-specific final vOTU FASTA.

RNA candidates receive CheckV assessments, merged into the existing quality summary
for representative selection. DNA and RNA assemblies use the shared
`sample_assembler_number_len_length` identifiers. `<sample>_assembly_provenance.tsv`
records the renamed member and representative IDs. Original assembler headers
remain in the optional-lookup `.ids.tsv` sidecars; downstream rules do not read
them. Assembly FASTAs, `.ids.tsv` sidecars and combination provenance are directly
under `03_CONTIGS/RNA/`, using sample-prefixed filenames such as
`<sample>_trinity.fasta` and `<sample>_combined_derreplicated.fasta`. Work directories, logs,
per-assembler geNomad results, VirSorter2 output and CheckV assessments stay
separate under `03_CONTIGS/RNA/<sample>/`. VirSorter2-positive sequences are in
`<sample>/virsorter/final-viral-combined.fa`, not `<sample>_combined_derreplicated.fasta`.
The helper rule `rna_viral_positive_fasta` exposes that result as a relative symlink
`03_CONTIGS/RNA/<sample>_virsorter2_RNA_viral_positive.fasta`, without duplicating
sequence data or rerunning VirSorter2 just to create the link. CheckV and the shared
vOTU combination use this symlink.
The RNA branch also writes `final-viral-score.tsv` there, but no separate ID-only
positive list.
After a successful RNA VirSorter2 job, only the final FASTA, score table and
boundary table (when produced) remain in its working directory. Intermediate
files, internal logs and `<sample>/virsorter.log` are removed. Failed jobs retain
their diagnostics. CheckV and the viral-positive symlink still use the same
final outputs.
The existing clustering tables track membership in the shared catalogue.

Existing FASTAs in the former `RNA/<sample>/<assembler>.fasta` layout are not
moved or filtered automatically. They need an explicit migration to reuse them
under the new filenames; otherwise the missing flat outputs can trigger assembly
jobs again. geNomad outputs also need refreshing for the new input basename.

### Final acceptance routes

Standard acceptance is unchanged: all three identification tools agree **or**
conservative geNomad support is present, together with the existing length/quality
gate (default length >=5,000 nt or the existing high-quality condition).

When RNA enrichment is enabled, a second route accepts a selected representative
if it is in the final representative-level VirSorter2 positive list, its
`max_score_group` is `RNA`, and its actual FASTA length is at least
`rna_min_contig_length` (default 500 nt). It does not require three-tool agreement
or a high CheckV quality estimate. Evidence on a different cluster member, or an
RNA assembler name alone, is insufficient. A standard-assembly contig can qualify
if the representative itself has this RNA evidence.

**Both routes still require the same mapping support:** filtered read count >1
in at least one of the configured key samples. These are the existing
unfiltered-catalogue mapping statistics, not the stricter microdiversity callability
thresholds. The RNA exception does not apply to subsampled catalogues.

The existing `vOTU_clustering_summary.tot.csv` now records `sequence_length`,
`VirSorter2_RNA_positive`, `mapping_support`, and `acceptance_route`
(`standard`, `RNA`, or `standard+RNA`) for retained representatives.

Abundance, annotation and optional microdiversity all use the same final
`filtered_...tot.fasta`. Phage-specific annotations are not thereby made applicable
to every RNA virus; their biological interpretation still needs care.

RNAviralSPAdes and MEGAHIT use `rna_assembly_threads` (16 by default), with
separate tool allowances of 32 and 16 GiB and reservations of 40 and 20 GiB.
Trinity uses `rna_trinity_threads` (8) and `rna_trinity_bfly_heap_gb` (20), passed
as `--CPU 8` and `--bflyHeapSpaceMax 20G`. Its main `rna_trinity_mem_gb` allowance
is 64 GiB, independent of the heap and total reservation. The default reservation
is 280 GiB (`286720` MiB), including concurrent heaps and overhead; at 16 allocated
workers it becomes 480 GiB. `--max_memory` is not an aggregate Java heap limit.
Use `--resources mem_mb=<available RAM in MiB>` to limit concurrent declared
reservations. Tool limits and reservations can be overridden independently;
raising only the reservation does not change the tool allowance.
Trinity runs without
read normalization and uses `--no_salmon` to skip its final Salmon expression filter.

Keep 8 Trinity workers as the starting point; the other RNA assemblers use their
separate 16-thread default. See [memory allowances and scheduling](memory_resources.md)
for `complexity=low|medium|high|extreme`, individual memory overrides and command
examples. Profiles scale DNA SPAdes/BBtools, while RNA allowances have their own
defaults. Python abundance normalisation is unchanged.

The assembled transcripts continue through the existing viral-candidate and vOTU filters.
Trinity uses `--full_cleanup` to remove its working directory
(`03_CONTIGS/RNA/<sample>/trinity_out`) after successful assembly. The rule copies
Trinity's sibling output
`trinity_out.Trinity.fasta` to a temporary FASTA, then filters and names it as
`03_CONTIGS/RNA/<sample>_trinity.fasta`; it does not
look for `Trinity.fasta` inside the working directory. After checking that the
exported FASTA is non-empty, it removes any remaining working directory, including
those retained by Trinity during resumed runs. Failed assemblies retain their
working files for resumption. Existing completed working directories are not
removed automatically: check their final FASTAs and active processes before
cleaning them manually, without rerunning assembly. SPAdes and MEGAHIT use
fresh working directories that are automatically removed when their jobs exit.
Leave `rna_trinity_strandedness` empty for mixed library types;
only set `RF` or `FR` if that orientation is appropriate for every input library.
Review `remove_euk` and contaminant filtering for your experiment. The optional
wrapper `rna_viruses` preset sets `remove_euk=False`; adding `--rna-enriched` to
another preset does not.

## Microdiversity outputs and interpretation

The branch also executes the separate `11_microdiversity_summary.py.ipynb`
notebook, producing `FIGURES_AND_TABLES/11_microdiversity_summary.tot.html` and
combined tables/figures. It does not depend on the abundance/lifestyle notebook.
See [the report documentation](abundance_lifestyle_summary.md) for targets and outputs.

Every sample is mapped competitively against the shared filtered catalogue with Bowtie2
end-to-end `--very-sensitive`. Results are under
`06_MAPPING/MICRODIVERSITY/<sample>/`:

- `aligned.sorted.bam`, its index, mapping log and `flagstat.txt`.
- `positions.tsv.gz`: one-based positions, quality-filtered A/C/G/T counts, reverse
  strand counts, callability, observed nucleotide diversity and Shannon entropy.
  This is sparse: positions absent from the pileup have no row.
- `snv_candidates.tsv`: non-reference alleles passing depth/count/frequency filters.
  Fixed or near-fixed reference differences are distinguished from polymorphic sites.
- `summary.tsv`: per-vOTU coverage, callable fraction, polymorphic-site count,
  mean callable-site diversity (`pi_callable`) and entropy, with analysis thresholds.

Defaults are depth >=100, base quality >=30, mapping quality >=30, alternate
count >=5 and frequency >=0.03. All are configurable. Unmapped, secondary,
supplementary, QC-failed and duplicate-flagged alignments are excluded. This does
not perform duplicate marking. Pysam overlap handling avoids counting overlapping
mates twice at the default quality threshold; improper/orphan paired alignments
are excluded. No both-strand requirement is imposed on RNA libraries; strand
counts are provided for inspection. BAQ correction is disabled.

Sites reaching `microdiversity_max_depth` (default 100,000) are conservatively
excluded from diversity/candidate calling and counted separately. Mean coverage
at capped sites can be underestimated. The maximum must exceed the minimum depth.
Only callable A/C/G/T reference sites contribute to mean diversity, **including
invariant sites**. No callable sites produces `NA`, not a misleading zero.

Per-site diversity is `pi = (n*n - sum(count_base**2)) / (n*(n-1))`.
Allele count/frequency thresholds govern candidate and polymorphic-site reporting;
they do not remove low-frequency bases from the pi calculation. These are
**observed read-level estimates and SNV candidates**, not error-corrected population
diversity or statistically validated variant calls. Sequencing/RT/PCR errors,
reference bias and cross-mapping between related viruses can affect results.
Compare callable coverage as well as pi. This branch does not implement haplotype
reconstruction, indel calling, consensus refinement or the complete nf-core workflow.

## Validation

Synthetic tests cover nucleotide diversity, invariant/uncallable sites, quality and
alignment-flag filtering, depth caps, sequence deduplication, and standard/RNA
acceptance and rejection cases. Two-sample Snakemake dry-runs check independent
flag combinations and RNA assembly through the existing CheckV clustering and
final filtering to microdiversity.

```bash
MOSAIC_SNAKEMAKE=/path/to/snakemake python -m unittest discover -s mosaic/tests -v
```

The test Python environment needs `pysam`, `pandas`, and `numpy`; CLI checks need `click` and
`ruamel.yaml` (or set `MOSAIC_CLI_PYTHON` to an interpreter with those dependencies).
The new assembler environments and a
full real-data run still need validation on the deployment system.

Migration from the earlier experimental implementation: the CD-HIT catalogue
script/environment have been removed. Any previously generated `05_vOTUs/RNA/`
results are not used or automatically deleted.

Tool references: [SPAdes modes](https://ablab.github.io/spades/running.html),
[MEGAHIT](https://github.com/voutcn/megahit),
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity),
[VirSorter2](https://github.com/jiarong/VirSorter2),
[pysam pileup API](https://pysam.readthedocs.io/en/latest/api.html).
