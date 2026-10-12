**MOSAIC** is a modular Snakemake workflow for viral and microbial genomics.
It is designed mainly for viral metagenomics, but also includes modes for
phage isolates, microbial metagenomes, bacterial isolates, and fungal isolates.

## Requirements

- **Conda**: A package manager for managing environments and dependencies
- **Mamba**: A fast, drop-in replacement for Conda, also used for managing environments and dependencies. It provides faster package resolution and environment creation compared to Conda.

You can learn how to install it here:

[Installing Mamba/Anaconda](https://www.notion.so/Installing-Mamba-Anaconda-15f1465d41b780c7be0bfa9b13fbf605?pvs=21)

## Installation

### 1. Clone the Repository

Create a dedicated folder for your bioinformatics applications to keep them organized. Navigate to this folder before cloning the repository.

```bash
mkdir ~/apps
cd ~/apps
```

Clone the MOSAIC repository to your local machine:

```bash
git clone https://github.com/lauramilena3/MOSAIC
cd MOSAIC/mosaic
```

### 2. Create Conda Environment

Create a conda environment with the required dependencies by running:

```bash
mamba env create -n Mosaic_workflow -f Mosaic.yaml
```

Activate the environment:

```bash
conda activate Mosaic_workflow
```

## Running MOSAIC

MOSAIC now includes a Python wrapper around Snakemake. The wrapper keeps the
command explicit, but users choose a biological workflow mode instead of
remembering long Snakemake target/config combinations.

From a source checkout, run the wrapper from the `mosaic/` directory:

```bash
python mosaic.py run MODE --raw /path/to/00_RAW_DATA -j 144
```

Mode names can be written with underscores or hyphens for convenience, for
example `phage_isolates` and `phage-isolates`.

Every run requires a raw FASTQ folder via `--raw` or `--input-dir`. If
`--results` is not provided, MOSAIC infers the results directory from the parent
of the raw data folder.

The wrapper always prints the full Snakemake command before executing it, so it
is clear which target and config values are being used. To print the command
without running Snakemake, use `--print-only`. To build the DAG without running
jobs, use `--dry-run`.

If Snakemake reports that the workflow directory is locked and no other MOSAIC
run is active, clear the stale lock with:

```bash
python mosaic.py unlock
```

```bash
python mosaic.py run phage_isolates \
  --raw /home/lmf/LEISE/ERWINIA_PHAGES_INSHERA/00_RAW_DATA/ \
  --kraken-db /home/lmf/db/KRAKEN/kraken/ \
  --ecc-memory 16000 \
  -j 144 \
  --dry-run
```

When MOSAIC is installed as a Python package, the same interface will become:

```bash
mosaic run phage_isolates --raw /path/to/00_RAW_DATA -j 144
```

## Main Workflow Modes

Use `python mosaic.py modes` to list all available wrapper modes. The main
biological modes are:

| Mode | Snakemake target | Use case | Default analysis flags |
| --- | --- | --- | --- |
| `viral_metagenome` | `runWorkflow` | End-to-end viral metagenome workflow from Illumina reads | `metagenome=True`, `isolates=False`, `microbial=False`, `sourmash_clean_reads=False`, `remove_euk=True` |
| `end-to-end` | `runWorkflow` | Hyphenated alias for `viral_metagenome` | Same as `viral_metagenome` |
| `phage_isolates` | `phage_isolates` | Phage isolate workflow from Illumina reads | `metagenome=False`, `isolates=True`, `microbial=False`, `sourmash_clean_reads=False`, `remove_euk=False` |
| `microbial_metagenome` | `microbial_metagenome` | Microbial metagenome workflow, including bacteria/fungi contexts | `metagenome=True`, `isolates=False`, `microbial=True`, `sourmash_clean_reads=True`, `remove_euk=True` |
| `short_read_bacteria` | `assembly_short_bacteria` | Short-read bacterial isolate assembly and analysis | `metagenome=False`, `isolates=True`, `microbial=True`, `sourmash_clean_reads=False`, `remove_euk=False` |
| `short_read_fungi` | `assembly_short_fungi` | Short-read fungal isolate assembly and analysis | `metagenome=False`, `isolates=True`, `microbial=True`, `sourmash_clean_reads=False`, `remove_euk=False` |

The old target names `assembly_phage` and `microbial` have been replaced by
`phage_isolates` and `microbial_metagenome` for consistency.

Initial viral identification, including subassemblies, uses geNomad. VirSorter2
still classifies clustered vOTU representatives later in `runWorkflow`; no
enabling flag is needed. The obsolete `--virsorter` / `--no-virsorter` wrapper
options and `VirSorter` config selector have been removed. Remove them from old
commands; subassembly FASTAs and CheckV inputs consistently use geNomad outputs.

The [memory guide](mosaic/docs/memory_resources.md) gives suggested profiles and
current Snakemake memory flags for isolates, RNA-enriched viromes, complex
microbiomes and very complex communities. Use `complexity=low|medium|high|extreme`
or the wrapper's `--complexity` flag; Python abundance normalisation is unchanged.

The [isolate recovery and purity workflow](mosaic/docs/isolate_purity.md)
documents inputs, outputs and decision thresholds, with a colour-coded overview
and seven step diagrams. It uses all QC-passed reads, preserves host-associated
exclusions and unexplained reads, and treats geNomad as evidence rather than a
viral selection gate.

The [viral metagenome manual](mosaic/docs/viral_metagenome.md) follows the same
step-by-step format, with a colour-coded overview and seven step diagrams. It
covers strict viral selection, the shared DNA/RNA catalogue, read subsets,
abundance and annotation outputs, optional branches, databases, memory and storage.

## Long-Read and Hybrid Modes

These modes are also available through the wrapper:

| Mode | Snakemake target | Required reads | Preset |
| --- | --- | --- | --- |
| `nanopore_assembly` | `assembly_nanopore` | Nanopore | `viral_metagenome` |
| `long_read_phage` | `assembly_long_only_phage` | Nanopore | `phage_isolates` |
| `nanopore_bacteria` | `assembly_nanopore_only_bacteria` | Nanopore | `bacteria_isolate` |
| `hybrid_nanopore_bacteria` | `assembly_nanopore_hybrid_bacteria` | Illumina plus Nanopore | `bacteria_isolate` |
| `pacbio_bacteria` | `assembly_pacbio_only_bacteria` | PacBio | `bacteria_isolate` |
| `hybrid_pacbio_bacteria` | `assembly_pacbio_hybrid_bacteria` | Illumina plus PacBio | `bacteria_isolate` |
| `pacbio_fungi` | `assembly_pacbio_only_fungi` | PacBio | `fungi_isolate` |
| `hybrid_pacbio_fungi` | `assembly_pacbio_hybrid_fungi` | Illumina plus PacBio | `fungi_isolate` |

## Stage Modes

Stage modes run a section of the viral metagenome workflow by default, but can
use another preset with `--preset` when appropriate.

| Mode | Snakemake target | Description |
| --- | --- | --- |
| `qc` or `only_qc` | `runQC` | Quality control and read-cleaning reports |
| `assembly` | `runAssembly` | QC plus short-read assembly |
| `viral_id` | `runViralID` | QC, assembly, and viral identification |
| `votu` | `runvOTUClustering` | Viral identification and vOTU clustering |
| `abundance` | `runAbundance` | Read mapping and abundance tables |
| `annotation` | `runAnnotation` | Viral annotation, taxonomy, and host prediction |

Example: QC for phage isolates, keeping the phage isolate defaults:

```bash
python mosaic.py run qc \
  --preset phage_isolates \
  --raw /path/to/00_RAW_DATA \
  --kraken-db /path/to/kraken_db \
  -j 144
```

## Read Requirements

Raw reads are mandatory for all wrapper modes. MOSAIC recognizes these default
file patterns:

| Read type | Expected file pattern |
| --- | --- |
| Illumina paired-end | `sample_R1.fastq.gz` and `sample_R2.fastq.gz` |
| Nanopore | `sample_nanopore.fastq.gz` |
| PacBio | `sample_pacbio.fastq.gz` |

The tags can still be overridden through Snakemake config values, for example:

```bash
python mosaic.py run viral_metagenome \
  --raw /path/to/00_RAW_DATA \
  --config forward_tag=1 \
  --config reverse_tag=2
```

## Important Flags

Most optional analyses default to `False` unless a preset enables them. The
wrapper passes these values explicitly to Snakemake so the run is reproducible.

- `remove_euk`: controls whether Kraken-classified eukaryotic reads are removed
  before later read-cleaning steps.
- `sourmash_clean_reads`: enables the shared per-sample and pooled final-clean-read
  GTDB report. Use `--sourmash-clean-reads`; it defaults to `True` only for the
  `microbial_metagenome` wrapper. It replaces the older `sourmash` option.
- `microbial`: enables microbial/fungal analysis behavior in the Snakefile.
- `visualization_tool`: chooses the small filtered-vOTU genome visualization
  tool. Defaults to `lovis4u`; use `clinker` to generate clinker HTML instead.
- `visualization_max_contigs`: renders the filtered-vOTU visualization only
  when the FASTA has fewer than this many contigs. Defaults to `50`.
- `metagenome` and `isolates`: internal workflow context flags. Users normally
  choose a wrapper mode instead of setting these directly.

## Passing Extra Snakemake Options

Any extra Snakemake options can be passed after `--`:

```bash
python mosaic.py run phage_isolates \
  --raw /path/to/00_RAW_DATA \
  --dry-run \
  -- --quiet rules
```
