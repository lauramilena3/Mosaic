# Memory allowances and scheduling

Choose memory from the expected community complexity, then check the benchmarks
after the first run. These are starting recommendations, not guarantees based
on the sample type. Sequencing depth, strain variation and coassembly can change
the requirement considerably.

Set `complexity=low|medium|high|extreme` to select the DNA SPAdes and BBtools
memory profile. The default is `medium`; isolate runs should explicitly select
`low`. Profiles change memory settings, not biological workflow modes. Each tool
allowance and scheduling reservation can also be overridden independently.
Python abundance normalisation is unchanged.

## Tool allowance versus job reservation

The tool allowance controls settings such as SPAdes `--memory`, MEGAHIT `-m` or
Java `-Xmx`. The job reservation includes that allowance plus overhead and is
used by Snakemake to decide which jobs can run together. A reservation alone
does not increase a tool's memory limit or enforce an operating-system limit.
See the [Snakemake resource documentation](https://snakemake.readthedocs.io/en/v7.18.2/snakefiles/rules.html#resources).

Use GiB for human-readable allowances and convert
reservations consistently to MiB: `mem_mb = GiB * 1024`. For example, an 80 GiB
reservation is `mem_mb=81920`. Convert the tool allowance separately to the
units expected by each tool; do not pass the larger reservation as its heap or
working-memory allowance. MEGAHIT receives bytes (`GiB * 1024 ** 3`), while
BBtools receives a Java heap in MiB. The rule commands remain explicit.

## Default allowances and reservations

| Tool / step | Tool allowance | Total job reservation | `mem_mb` |
| --- | --- | --- | --- |
| DNA SPAdes and SPAdes assembly-depth tests | 64 GiB | 80 GiB | 81920 |
| RNAviralSPAdes | 32 GiB | 40 GiB | 40960 |
| MEGAHIT | 16 GiB | 20 GiB | 20480 |
| BBNorm read normalisation | 6 GiB Java heap | 8 GiB | 8192 |
| BBCountUnique / k-mer rarefaction | 6 GiB Java heap | 8 GiB | 8192 |
| Trinity, 8 workers | 64 GiB main allowance; 20 GiB per Butterfly process | 280 GiB | 286720 |
| Standard DNA/RNA QUAST | No tool-enforced limit | 8 GiB | 8192 |
| Assembly-depth QUAST | No tool-enforced limit | 12 GiB | 12288 |
| Hybrid SPAdes | 350 GiB, pending suitable benchmarks for a reduction | 440 GiB | 450560 |

Most reservations add approximately 25% overhead, rounding small reservations
upwards. These values are based on sampled peak RSS from existing runs, not
virtual-memory size. Moderate virome runs peaked around 37 GiB for DNA SPAdes,
11 GiB for RNAviralSPAdes, 4 GiB for MEGAHIT and 6.2 GiB for BBNorm.
The DNA and BBtools entries are the `medium` defaults; the other profile values
are below. RNA defaults do not scale automatically with complexity: raise their
individual allowances when the RNA assemblies require more RAM. Python
abundance-normalisation notebooks have no new tool limit or reservation, and
existing behaviour stays intact.

BBNorm-normalised paired and unpaired FASTQs are temporary assembly inputs.
Snakemake removes them after their scheduled consumers, including read counting,
finish. Raw and clean reads, read-count tables, histograms and peak reports remain
available. This does not change abundance normalisation.

### Trinity workers

Use **8 workers per sample** as the starting point. This is already the supplied
`rna_trinity_threads` default; RNAviralSPAdes and MEGAHIT keep their separate
`rna_assembly_threads` default of 16. The run-wide `-j` value controls concurrent
allocated cores, not the worker count of each Trinity job.

Keep `rna_trinity_bfly_heap_gb=20`. The earlier Butterfly failures were Java heap
failures; increasing only the main memory allowance does not increase that heap.
Trinity's `--max_memory` applies where memory limiting is supported and is not a
whole-job ceiling. See the [Trinity documentation](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity).

The Trinity reservation is deliberately conservative:

```text
(main allowance + workers * Butterfly heap) * 1.25
8 workers:  (64 + 8 * 20) * 1.25 = 280 GiB
16 workers: (64 + 16 * 20) * 1.25 = 480 GiB
```

These are parallel-heap allowances, not observed requirements. Existing RNA runs
peaked around 21 GiB overall, but future partitions may use more heap at once.

## Suggested analysis profiles

Select one of these profiles explicitly. MOSAIC does not detect complexity or
infer a profile from the sample type or project name. Omitting `complexity`
uses `medium`, including for isolates.

| Profile / analysis | Main DNA SPAdes allowance | Main DNA job reservation | BBNorm / k-mer heap and reservation |
| --- | --- | --- | --- |
| `low`: phage isolates, one expected dominant genome with possible host contamination | 32 GiB | 40 GiB (`40960`) | 6 GiB heap / 8 GiB reserved |
| `medium`: RNA-enriched or other moderately complex viromes | 64 GiB | 80 GiB (`81920`) | 6 GiB heap / 8 GiB reserved |
| `high`: complex viromes or microbial communities, including leaf-associated microbiomes | 128 GiB | 160 GiB (`163840`) | 16 GiB heap / 24 GiB reserved |
| `extreme`: very complex communities, such as soil or highly diverse metagenomes | 450 GiB | 550 GiB (`563200`) | 16 GiB heap / 24 GiB reserved |

With separate reservations, a 16 GiB BBNorm/k-mer heap gets 24 GiB
(`24576`) per job. Isolate recommendations are provisional: many existing isolate
benchmarks lack usable memory measurements. Complex runs reached about 90 GiB
for DNA SPAdes; very complex runs reached about 410 GiB, and assembly-depth tests
reached about 374 GiB. Do not apply the 64 GiB baseline to those large datasets.
Length or contig count alone is not a reliable predictor of assembly memory.

## Running with profiles

Run from `mosaic/`, replacing the input path. The examples use a run-wide budget
of `900000` MiB (about 879 GiB) for a high-memory server. Choose a budget that
fits the server's available RAM, leaving headroom for the operating system,
unrestricted notebooks and other users. It is not a universal budget for every
server. Use `free -m` to inspect total and currently available memory.

Each example starts as a dry-run; remove `-n` to execute. Add optional biological
analysis flags only when required for the experiment.

### Phage isolates

```bash
snakemake --use-conda -p phage_isolates \
  --config input_dir=/path/to/project/00_RAW_DATA isolates=True metagenome=False complexity=low \
  --resources mem_mb=900000 \
  -j 128 -k --rerun-incomplete -n
```

The profile sets `ecc_memory` to 6144 MiB and reserves 8192 MiB for BBNorm and
k-mer rarefaction. Isolate mode uses all QC-passed reads for mapping and does not
enable the RNA assembly branch.

### RNA-enriched viromes

```bash
snakemake --use-conda -p runWorkflow \
  --config input_dir=/path/to/project/00_RAW_DATA RNA_enriched=True complexity=medium \
  --resources mem_mb=900000 \
  -j 128 -k --rerun-incomplete -n
```

This uses 16 threads for RNAviralSPAdes/MEGAHIT and 8 Trinity workers. RNAviralSPAdes
gets 32 GiB with a 40 GiB reservation; MEGAHIT gets 16 GiB with a 20 GiB reservation.
Trinity gets `--max_memory 64G` and `--bflyHeapSpaceMax 20G`, with a 280 GiB
reservation. At 16 allocated Trinity workers the reservation becomes 480 GiB.
Changing only `mem_mb` does not change any of these tool allowances.

### Complex viromes and microbial communities

```bash
snakemake --use-conda -p runWorkflow \
  --config input_dir=/path/to/project/00_RAW_DATA complexity=high \
  --resources mem_mb=900000 \
  -j 128 -k --rerun-incomplete -n
```

For microbial rather than viral analysis, replace `runWorkflow` with
`microbial_metagenome` and add `microbial=True sourmash_clean_reads=True` to `--config`.
The memory profile does not change the biological workflow mode. Add
`RNA_enriched=True` only when the additional RNA assemblies are wanted.

### Very complex communities

```bash
snakemake --use-conda -p runWorkflow \
  --config input_dir=/path/to/project/00_RAW_DATA complexity=extreme \
  --resources mem_mb=900000 \
  -j 128 -k --rerun-incomplete -n
```

At this budget, two main DNA jobs reserving 563200 each cannot overlap. Other
jobs may still run if their declared reservations and core requirements fit.
For soil microbial profiling, use the microbial target/config adjustment above.

### Wrapper example for a 1 TB server

Use `extreme` for very complex samples, with a starting scheduling budget of
`900000` MiB (about 879 GiB). Check currently available RAM and other server users
before choosing that budget. With 28 available cores, for example:

```bash
python mosaic.py run viral_metagenome \
  --raw /path/to/project/00_RAW_DATA \
  --complexity extreme \
  --subassembly \
  --assembly-stats \
  --subsampling \
  --vcontact \
  --dram \
  --metavr-blast \
  --map-to-refseq \
  --config min_votu_length=10000 \
  --config microbial_spacers=/path/to/spacers.fa \
  -j 28 --dry-run \
  -- --resources mem_mb=900000
```

Remove `--dry-run` to execute. These optional analysis flags reproduce the
extended virome analysis; omit any that are not required. The spacer FASTA must
already exist. Add `--kraken-db /path/to/database` only to select a different
database from the configured default. Omitting `--ecc-memory` lets the extreme
profile supply its 16 GiB BBtools heap. Updating the checkout is required before
using `--complexity`; no obsolete classifier-selection flag is needed.

## Individual overrides

`null` in `config.yaml` means use the profile or the tool's default. An explicit
allowance takes precedence; if its reservation remains `null`, the reservation
scales using the selected profile/tool's headroom ratio, rounded up to GiB.
For example, a 96 GiB DNA allowance under `medium` reserves 120 GiB.

| Config parameter | Meaning |
| --- | --- |
| `assembly_mem_gb`, `assembly_mem_mb` | DNA SPAdes allowance and reservation; also used by assembly-depth tests and coassembly. |
| `ecc_memory`, `bbtools_mem_mb` | BBNorm/BBCountUnique Java heap in MiB and separate total job reservation. |
| `rna_spades_mem_gb`, `rna_spades_mem_mb` | RNAviralSPAdes allowance and reservation. |
| `rna_megahit_mem_gb`, `rna_megahit_mem_mb` | MEGAHIT allowance and reservation. |
| `rna_trinity_mem_gb`, `rna_trinity_mem_mb` | Trinity main allowance and total reservation including concurrent Butterfly heaps. |
| `rna_trinity_bfly_heap_gb`, `rna_trinity_threads` | Java heap per Butterfly process and requested worker count. |
| `hybrid_assembly_mem_gb`, `hybrid_assembly_mem_mb` | Both individual and pooled hybrid SPAdes allowances/reservations. |
| `quast_mem_mb`, `quast_depth_mem_mb` | Standard DNA/RNA and assembly-depth QUAST reservations; no new tool-enforced limit. |

The legacy `rna_assembly_mem_mb` is an optional shared RNA tool-allowance fallback,
now consistently interpreted as MiB, rounded up to whole GiB for tool arguments.
An explicit per-tool `rna_*_mem_gb` takes precedence. Remove old shared overrides
to use the separate 32/16/64 GiB defaults. `ecc_memory=6000` remains a valid
explicit heap override; otherwise profiles supply 6144 or 16384 MiB.

For example, add `assembly_mem_gb=96` to `--config` to increase the DNA allowance.
Add `assembly_mem_mb=143360` to reserve 140 GiB without changing that allowance.
Alternatively, `--set-resources shortReadAsemblySpadesPE:mem_mb=143360` changes
the reservation of that rule only, not assembly-depth tests. The old
`shortReadAsemblySpadesPE:mem_gb=...` resource override no longer controls SPAdes:
use the config allowance instead. Do not lower a reservation below its allowance.

The wrapper accepts `--complexity low|medium|high|extreme`, or repeated
`--config KEY=VALUE` overrides. Pass the global budget after `--`, for example:

```bash
python mosaic.py run phage_isolates --raw /path/to/project/00_RAW_DATA \
  --complexity low -j 128 --dry-run -- --resources mem_mb=900000
```

The run-wide budget constrains declared resources only. Other tools, unreserved
notebooks and other users can consume additional RAM; this is not a whole-server
enforcement mechanism. Keep the budget at least as large as the largest scheduled
job reservation: Snakemake 7.18 can cap a reservation at the global budget without
reducing the tool allowance. For a smaller server, lower the tool allowance and
its reservation rather than relying on the budget to reduce a single job's RAM.
Do not add a blanket `--default-resources mem_mb=...` to
these examples, since it would assign new reservations to otherwise unrestricted
notebook rules.

## Check and adjust

Use the [benchmark summary](benchmarking.md) after the run. Check valid sampled
`max_rss` values, not `max_vms`, and do not interpret missing or implausibly low
values as evidence that a job needs almost no RAM. Keep headroom above observed
peaks; benchmark sampling can miss short-lived peaks and completed runs do not
capture every failed job's requirement.

Raise the tool allowance and its reservation together when a tool reports memory
exhaustion. Reduce concurrent jobs or the available scheduling budget when total
server RAM is under pressure. A larger `-j` does not grant an individual tool more
memory. Changing a tool allowance changes rule parameters/commands and can
legitimately trigger a rerun; inspect `-n` before changing profiles on existing
results. This guide does not touch completed outputs or bypass rerun checks.

## Storage after assembly

Trimmed paired FASTQs and the merged orphan FASTQ remain persistent for later
reuse. The two separate fastp orphan files are temporary because their reads are
preserved in the merged file. BBNorm-normalised FASTQs are temporary assembly
inputs. Successful MMseqs clustering keeps representative sequences and membership
tables, but removes redundant `*_all_seqs.fasta` exports (the DRAM export is
already declared temporary).

Optional per-base mapping exports default to off (`per_base_coverage: False`).
Add `--config per_base_coverage=True` to save gzipped nonzero-position coverage
for assemblies, final vOTUs, all-assembled, RefSeq, user references, contaminants,
microbial contigs and host mappings. Counts, RPKM, breadth and mean-depth summaries
are always saved. Isolate activity/discovery reports have a separate
`host_prophage_activity` flag and retain their required gzipped host-depth input
while enabled. See [RefSeq mapping](refseq_mapping.md) and
[isolate activity](isolate_purity.md#host-prophage-activity).

These settings affect future successful jobs. Editing them does not clean
previously completed project files by itself; old plain-text coverage files
and redundant exports are not migrated automatically.
