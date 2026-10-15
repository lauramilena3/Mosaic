# Assembly identifiers

One final FASTA is retained per assembly. RNA FASTAs are directly under `03_CONTIGS/RNA/`; DNA locations are unchanged. There
is no `RENAMED` directory and no retained original-ID copy of that FASTA.

The assembler produces a temporary `*.unrenamed.fasta`. The shared
`rename_assembled_contigs` rule in `02_assembly_short.smk` writes the final
FASTA and a neighboring `*.ids.tsv` table. Snakemake removes the temporary
input after successful processing. The helper uses Biopython already available
in the workflow environment; no new environment is needed.

MOSAIC writes saved assemblies and its own extracted/reference FASTAs with one
header line and one complete sequence line per record (`fasta-2line`). FASTA
splits use `seqkit -w 0`, and geNomad candidate exports are unwrapped before
use. This changes only line wrapping, not sequence content, IDs or descriptions.
Downloaded databases and third-party tools' internal working files can retain
their native formatting.

```text
03_CONTIGS/{sample}_spades_filtered_scaffolds.tot.fasta
03_CONTIGS/{sample}_spades_filtered_scaffolds.tot.ids.tsv
03_CONTIGS/RNA/{sample}_rnaviralspades.fasta
03_CONTIGS/RNA/{sample}_rnaviralspades.ids.tsv
03_CONTIGS/RNA/{sample}_megahit.fasta
03_CONTIGS/RNA/{sample}_megahit.ids.tsv
03_CONTIGS/RNA/{sample}_trinity.fasta
03_CONTIGS/RNA/{sample}_trinity.ids.tsv
```

Identifiers follow `{sample}_{assembler}_{number:0{width}d}_len_{length}`, for
example `sample01_metaspades_00001_len_1000`. The counter width is
`max(5, len(str(total_contigs)))`, calculated separately for each assembly.
An assembly of 200,000 contigs therefore starts at `000001` and ends at
`200000`; an assembly of 2,000,000 starts at `0000001`. An initial counting
pass determines the width without loading all sequences into memory. For RNA,
only sequences of at least `rna_min_contig_length` are counted, numbered and saved;
shorter sequences are omitted before renaming or downstream analysis. Counters
follow FASTA order and restart at 1 per assembly. Length is calculated from the
actual sequence, including Ns. Renaming preserves sequence content and case. Existing DNA length/coverage
filtering happens first and is unchanged.

SPAdes assemblies use `metaspades` when the actual command uses `--meta`,
otherwise `spades`. Subsampled assemblies add `_sub` to the sample prefix;
assembly-depth tests add `_{percentage}pct`. Long-read assembly copies use the
actual assembler, and polished stages include their stage in the prefix.

ID tables record `contig_id`, `sample`, `assembler`, `original_id`, `length_bp`
and `original_description`. Assembly graphs keep their original identifiers;
the ID table provides the link to the old FASTA headers. No separate retained
FASTA is required to recover those headers.

The `.ids.tsv` files are lookup-only sidecars. No downstream rule reads or
requires them. Combination, classification, clustering, mapping and metadata
joins use the renamed FASTA identifiers directly. Original assembler headers
are kept only in the sidecars, not copied into the combined provenance, top
metadata or cluster membership tables.

RNA combination and all-assembled concatenation retain the named contig IDs.
The per-sample RNA pool is `03_CONTIGS/RNA/{sample}_combined_derreplicated.fasta`.
Existing RNA exact-sequence/reverse-complement deduplication keeps the first
encountered named ID. Its provenance table records both member and representative
IDs (`contig_id` and `representative`). All-assembled provenance records
`contig_id`, `sample`, `assembler` and `length_bp`, derived from the named FASTAs
and their source paths. MMseqs/CheckV clustering and the final filtered vOTU catalogue filename
remain unchanged.
If a tool extracts or trims a region, the retained identifier refers to its
parent assembly; use the measured sequence length in that tool's output for
the processed sequence.

The isolate-summary notebook reuses this existing `contig_id`/`sample`
provenance lookup before loading contig evidence. Provirus regions inherit
their parent contig's sample. Standardized identifiers and legacy `_NODE_`
identifiers are also recognized when a contig is not in that lookup; no
`.ids.tsv` dependency is added.

## Full-assembly geNomad evidence

geNomad runs on each sample's named DNA FASTA and, when RNA enrichment is
enabled, each named RNA assembler FASTA. These results cover the assemblies
before abundance selection. `map_to_all_assembled=True` reuses these existing
per-assembly results; it does not request a second run on the combined FASTA.

Selected cluster metadata collects geNomad results and terminal-repeat reports
directly by `contig_id`. The recorded assembly sample and assembler identify
the original result directory; original headers and `.ids.tsv` lookups are
not needed. Missing classifications are recorded as `not reported`. Existing
`04_VIRAL_ID/all_assembled_geNomad_tot/` outputs are no longer used and are not
automatically deleted. VIBRANT, VirSorter, CheckV, Pharokka and optional Cenote
retain their selected-cluster-representative scope. geNomad uses `--restart`
to regenerate per-assembly results when an input assembly changes.

## Existing projects

Existing assembly FASTAs require a one-time conversion with neighboring ID
tables. Convert only approved projects, while no workflow is using those
files. Conversion uses a same-directory temporary file followed by atomic
replacement, verifies sequence content/order/length, and retains no backup
FASTA. The saved original descriptions allow original headers to be restored
if required.

After conversion, rerun downstream classification, combined catalogues,
mapping, normalization and reports. Previous results contain old IDs and
must not be mixed with renamed assemblies. Assemblies themselves need not
be rebuilt solely for conversion.
Clear the old assembler's Snakemake metadata for each converted final FASTA
with `snakemake --cleanup-metadata <converted FASTA paths>`; otherwise its
old producer/code record can cause an unnecessary assembly rerun. Leave
preprocessing and downstream metadata intact.
