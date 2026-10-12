# Rule/environment audit

Checked on 2026-10-05. This audit compares commands in the current rules with
their assigned YAML files and, where available, matching installed Snakemake
Conda environments. Dependencies supplied indirectly by another package count
as available; absence from a YAML's direct package list alone is not a failure.

## Fixed in this change

All five SPAdes rules now use `env3.yaml`: `shortReadAsemblySpadesPE`,
`hybridAsemblySpades`, `hybridAsemblySpadesPooled`, `metaspadesPE_test_depth`, and
`rna_assemble_spades`. `env2.yaml` is an iPHoP environment and does not provide
SPAdes. The shared environment pins SPAdes 4.3.0 and now includes seqtk 1.5 for
assembly filtering. The redundant `rna_spades.yaml` definition (SPAdes 4.0.0)
was removed; existing installed environments were not deleted.

The installed SPAdes 4.3.0 executable reports both `--meta` and `--rnaviral` in
its help. A fresh Linux-64 Conda **dry-run solve** of the updated `env3.yaml`
succeeded with flexible channel priority. No environment was installed or
modified by that check; the pip section and a full assembly run were not tested.
With `--use-conda`, Snakemake will create a new environment for the updated YAML.

## Requested dependency fixes

The following omissions found in the original audit are now fixed in the YAML
definitions. Existing installed environments were not modified; Snakemake will
create new hashed environments when these definitions are used with `--use-conda`.

| Environment | Added packages | Rules that need them |
| --- | --- | --- |
| `env6.yaml` | BLAST 2.17.0 | `vOUTclustering`, `combine_with_taxmyphage`, `clustering_isolates` |
| `env1_mapping.yaml` | bcftools 1.24 | `call_SNPs_sub` |
| `env4.yaml` | seqtk 1.4, git 2.38.1 | `merge_microbial`, `getKrakenTools`; the legacy PhaGCN downloader has since been replaced |
| `satellite_finder.yaml` | seqtk 1.5 | `satellite_finder_get_fasta` |
| `viga.yaml` | seqtk 1.3 | `fasta_to_a2m` |
| `vcontact.yaml` | dos2unix 7.5.2, OpenJDK 11 | `clusterTaxonomy`, including ClusterONE |
| `bacterial.yaml` | OpenJDK 11 | `predict_spacers`, using the downloaded MinCED JAR |
| `wtp.yaml` | wget 1.21.4 | `downloadGtdbtk_db`, `get_vcontact2` |

`env5.yaml` already contains BLAST 2.17.0 and did not need an edit. BLAST was
added to `env6.yaml` following confirmation of the environment-name correction.

All eight amended environments, plus the existing `env3.yaml` and `env5.yaml`,
passed fresh-prefix Linux-64 Conda dry-run solves with flexible channel priority.
Every declared Conda package constraint was checked against the resolved
package records using Conda's `MatchSpec`. No existing dependency declarations
were removed or changed. Only package indexes were downloaded into `/tmp`;
these checks did not install environments or change global Conda settings.

Compatibility choices preserve the existing analysis tools and Python versions:

- `env4.yaml`: seqtk 1.5 requires newer zlib than legacy PyTables permits;
  seqtk 1.4 resolves with the existing Python 3.7/NumPy/PyTables pins. Git 2.38.1
  resolves with the Perl 5.32 runtime required by the existing MMseqs2 package.
- `viga.yaml`: seqtk 1.5 conflicts with the pinned GCC runtime. seqtk 1.3
  resolves while retaining Python 2.7, libgcc-ng 11.2.0, zlib 1.2.13 and OpenSSL
  1.1.1q. The required `seqtk subseq` command is available in that version.
- `wtp.yaml`: wget 1.25.0 requires zlib 1.3, while the Mash 2.3 builds compatible
  with the pinned GSL 2.6 need zlib below 1.3. wget 1.21.4 resolves without
  changing the Mash/GSL pins.
- bcftools 1.24 matches the existing samtools 1.24 pin. OpenJDK 11 is used in
  both Java-dependent environments; the existing MinCED and ClusterONE JARs
  were also successfully loaded with an installed Java 11 runtime.

## PhaBOX taxonomy migration (2026-10-12)

`PhaGCNTaxonomy` now uses `phabox.yaml` and the PhaBOX 2.2 taxonomy-only task.
`downloadPhaBOXDB` reuses `env5.yaml` for wget and Python ZIP extraction;
`env4.yaml` is unchanged. A minimal Conda dry-run containing PhaBOX 2.2 and
BACPHLIP 0.9.6 failed: their scikit-learn/Python dependencies cannot coexist.
The dedicated environment keeps the existing BACPHLIP and MMseqs rules separate.
See [PhaBOX taxonomy](phabox_taxonomy.md) for versions, paths and report fields.

## Assembly-depth QUAST

`assemblyStatsILLUMINA_test_depth` and `viralStatsILLUMINA_test_depth` now use
`env3.yaml` and invoke `quast.py` from that environment. They no longer take
`config["quast_dir"]` as an input or depend on `tools/quast-5.0.2`.

The existing report directories and text-report filenames are preserved. Each
rule additionally declares its `transposed_report.tsv`, writes a separate log,
uses its allocated threads, skips empty FASTAs and labels the remaining
assemblies by filename. `--min-contig 1` reports every saved contig, matching the
RNA QUAST rule; the upstream assembly length filtering is unchanged. If all
inputs are empty, the rule writes an empty table and a short explanatory report
without running QUAST. Commands remain explicit inside each rule block. The
normal DNA and RNA QUAST rules are unchanged.

Twelve focused tests passed, covering the package declarations, actual QUAST
5.3.0 execution for both depth rules, all-empty inputs, real Snakemake execution
on existing depth FASTAs, unchanged DNA/RNA report behavior and benchmarks.

## Host tools and externally installed programs

These are distinct from the missing analysis packages above:

- Source-build rules including `get_minced`, `get_mmseqs`, `get_VIGA` and
  `get_ALE` rely on host build tools such as make and, where required by the
  source, compilers/JDKs. Their assigned environments do not guarantee the
  complete build toolchain.
- `get_ALE`, `get_weeSAM`, `getClusterONE` and `downloadCanu` do not declare a
  Conda environment; their shell commands rely on the launch environment's
  git/download/build tools as applicable.
- The bacterial environment's existing pip section downloads VAMB from a Git
  URL. Its fresh pip installation still needs access to Git and the network in
  the launch environment; those pip declarations were not changed.

## Scope and limits

The audit covered direct tool invocations in rule blocks, including optional
short-read, RNA, hybrid/long-read and assembly-depth branches. Repository tools
and tools executed inside containers were distinguished from commands expected
on the Conda PATH. Commented commands were excluded: for example, the commented
git clone in `downloadVirSorterDB` is not a git requirement of that rule.

Conda solves do not install or validate the pip sections. The existing bacterial
cache passed `pip check`; the existing VIGA cache reports missing requirements
for extra `gravity`/`comparem` packages that are not declared in `viga.yaml`.
Those unrelated cached packages were not changed. Full fresh pip installations,
Python imports in every notebook, hidden tool subprocesses beyond the Java
cases above, database downloads and complete biological runs were not
exhaustively executed. QUAST itself does not accept paths containing spaces.

Changing an environment definition can cause Snakemake's normal provenance
checks to schedule existing jobs again, including download rules that share
that environment. No existing project outputs or Snakemake metadata were
changed or cleaned up by this implementation.
