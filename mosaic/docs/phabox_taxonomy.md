# PhaBOX 2 taxonomy

MOSAIC now runs the PhaGCN taxonomy module through PhaBOX 2. The old
`PhaGCN_newICTV` repository is unavailable; fresh installations no longer need
to clone it. This is a taxonomy update, not another viral-identification step.

## Environment and database

`envs/phabox.yaml` pins PhaBOX 2.2, Python 3.10, NumPy 1.26.4 and a CPU-only
PyTorch 2.4 build. Snakemake creates the environment with `--use-conda`.
`env4.yaml` remains unchanged for BACPHLIP and the existing MMseqs rules.
PhaBOX and BACPHLIP require different scikit-learn versions and cannot share
their packaged environment.

`downloadPhaBOXDB` uses the existing `env5.yaml` download environment:

```yaml
phabox_db: "db/PHABOX/phabox_db_v2_2"
phabox_db_url: "https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2_2.zip"
phabox_db_sha256: "4bdcd1954fff40a84784218828c475f3470cabb831fe65ec74e557d3525f9bc6"
phabox_min_length: 2000
```

The ZIP is approximately 705 MB (672 MiB). Its official SHA-256 is checked
before extraction. The rule finds the database directory inside the archive,
checks the required taxonomy files and writes `download_info.tsv` with the
database version, UTC download time, URL and actual checksum. The ZIP and
extraction workspace are temporary. The extracted database is kept and reused
by all samples and projects sharing this workflow installation.

Database paths are relative to the workflow working directory, as for other
MOSAIC databases. Use an existing compatible v2.2 database through `phabox_db`
if it is already installed elsewhere. If overriding the download URL, also
supply that archive's checksum; an empty `phabox_db_sha256` disables the hash
comparison, not ZIP integrity or required-file checks.

See the [official PhaBOX release instructions](https://github.com/KennthShang/PhaBOX#-the-most-recent-update-logs).

## Classification and outputs

The existing `PhaGCNTaxonomy` rule name is retained. It now runs only
`phabox2 --task phagcn`, with the allocated `{threads}` and the existing
2,000 bp minimum classification length. Changing `phabox_min_length` changes
which records receive this annotation, not which records exist in the catalogue.
No PhaMer filtering, CHERRY host prediction, PhaTYP lifestyle prediction or
PhaBOX vOTU clustering is added. MOSAIC's existing clustering remains unchanged.

Each job starts in a new temporary output directory so PhaBOX cannot reuse
filtered contigs or alignments from an older FASTA. Working files are removed
on success and failure. MOSAIC retains:

| Output under `07_ANNOTATION/` | Contents |
| --- | --- |
| `PhaGCN_taxonomy_report_<sequence>.csv` | Existing filename and `contig_name`, `prediction`, `score` columns, plus PhaBOX rank/source fields. |
| `PhaBOX_<sequence>/phagcn_prediction.tsv` | Unmodified PhaBOX prediction table, including its full lineage, raw scores and genus grouping. |
| `PhaBOX_<sequence>.log` | Tool output or a short note for an empty catalogue. |

The compatible CSV's `prediction` is the explicitly reported family, not a
genus substituted for a missing family. `score` is the family-level score from
the original per-contig lineage, only when that family still matches the final
assignment. PhaBOX can propagate taxonomy within its network; if the final
family differs from the original scored family, MOSAIC keeps the final taxonomy
but does not invent a family confidence. Original/final lineage and score strings
remain available for review. Missing evidence is `not reported`.

`vOTU_metadata_polished.csv` keeps `PhaGCN_family` and `PhaGCN_score` in their
existing topic order, followed by the new `PhaBOX_` fields. Realm through species,
genus cluster, classification status, source versions, database path and cutoff
are retained where reported. Short/unclassified records remain in the metadata.
Legacy PhaGCN CSVs can still be read; new PhaBOX fields are not fabricated for
those older results. Results from the old 2022 classifier and the new database
are not assumed to be identical.

## Running and updating

No new enable flag is needed: targets that already required PhaGCN taxonomy
now use PhaBOX. The short-read `phage_isolates` target does not gain this step.
To refresh an existing final catalogue explicitly, from `mosaic/`:

```bash
snakemake --use-conda -p runAbundanceSummary --config \
  input_dir=/path/to/00_RAW_DATA \
  --forcerun PhaGCNTaxonomy --cores 32 --rerun-incomplete -n
```

Keep the same biological/optional flags as the original run. Remove `-n` to
execute. This target regenerates taxonomy and the integrated abundance report,
plus any missing prerequisites. Do not touch failed outputs to bypass the
missing-repository error. No existing cached PhaGCN checkout or project results
are deleted by this migration.
