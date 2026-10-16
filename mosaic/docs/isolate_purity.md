# Phage isolate recovery and purity

The isolate workflow assesses genome recovery, purity and relatedness across a collection of sequenced phage isolates. The expected result is one dominant viral genome per sample, but fragmentation, host DNA and possible mixtures remain visible for review. A warning does not automatically invalidate a recovered genome.

Unlike the virome workflow, it does not require a geNomad viral-positive call. All saved SPAdes contigs enter the catalogue; retained/excluded is a coverage-and-assigned-host decision. Novel or unclassified sequences are not removed simply because a viral classifier did not recognise them.

## Quick start

Run from `mosaic/`:

```bash
snakemake --use-conda -p phage_isolates \
  --config \
    input_dir=/path/to/project/00_RAW_DATA \
    isolates=True metagenome=False \
    sourmash_clean_reads=True sourmash_contig_catalogue=True \
  -j 128 -k --rerun-incomplete -n
```

Remove `-n` to execute. This example enables both optional Sourmash screens, with the default 5,000 bp cutoff for cluster representatives. Omit those two flags for the core isolate analysis. `isolates=True` enforces DNA assembly mode and bypasses biological read removal, regardless of the virome defaults for `metagenome`, `RNA_enriched` or `remove_euk`.

See [memory allowances and scheduling](memory_resources.md#phage-isolates) for
`complexity=low` isolate memory profile and a command with separate DNA SPAdes
allowance and scheduling reservation. Without a selection, `complexity=medium`
is used, including for isolates. Python abundance normalisation is unchanged.

Add `map_to_RefSeq=True` for independent RefSeq read mapping/detection, `metavr_blast=True` for METAVR similarity searches, or `host_identification_test=True` for the all-host comparison. These are not prerequisites for the core recovery/purity assessment. RefSeq BLAST of cluster representatives is already part of the core report; it is separate from the optional read-mapping branch.

## Inputs and project layout

```text
project/
├── 00_RAW_DATA/
│   ├── IsolateA_R1.fastq.gz
│   ├── IsolateA_R2.fastq.gz
│   ├── IsolateB_R1.fastq.gz
│   └── IsolateB_R2.fastq.gz
├── HOST/
│   ├── HostA.fasta
│   └── HostB.fasta
└── host_mapping_file.tsv
```

`host_mapping_file.tsv` is tab-separated:

```tsv
sample	host
IsolateA	HostA
IsolateB	HostB
```

| Input | Requirement |
| --- | --- |
| Paired Illumina FASTQs | Required. Default names are `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`; change `forward_tag`/`reverse_tag` if needed. Sample IDs are discovered from the reverse-read filenames. |
| `HOST/{host}.fasta` | Supply the available bacterial host assemblies directly under `HOST/`. The host ID is the filename without `.fasta`. |
| `host_mapping_file.tsv` | Required whenever host FASTAs exist. Every analysed sample needs one assignment to an available host. Several samples can use the same host. The file can instead be under `HOST/`; if both copies exist they must agree. |
| `config.yaml` and Conda | Use the repository configuration and `--use-conda`. The command-line options override the YAML defaults. |
| Databases and external tools | Missing default databases/tools are provided by their existing download rules. First execution needs network access and sufficient disk space; existing installations are reused. Custom database paths must refer to usable installations or outputs supported by a download rule. |

By default, results are written beside `00_RAW_DATA`, in its parent project directory. `results_dir=/path/to/project` can override that location; place `HOST/` and the host assignment table under the selected results directory. Runs without host FASTAs are allowed, but host origin is explicitly unassessed and the isolate decision requires review. There is no inferred host assignment from a sample name.

### Databases and tools

Existing download rules provide these shared dependencies when needed. Sourmash downloads are only required when either optional Sourmash screen is enabled.

| Download rule | Used for |
| --- | --- |
| `downloadKrakenDB`, `getKrakenTools` | Read-level contamination reporting. |
| `downloadGenomadDB` | Assembly and host classification evidence. |
| `downloadCheckvDB` | Isolate-contig and host viral-region quality assessment. |
| `downloadCheckMDB` | Host assembly completeness/contamination. |
| `downloadRefSeqViral` | Core cluster-representative BLAST; also reused by optional RefSeq read mapping. The downloader keeps a dated snapshot and undated links. |
| `downloadTaxmyphageDB`, `downloadPholdDB` | Taxonomy and Phold protein annotation of exported enriched host loci when host activity is enabled. The Phold downloader installs both the reference database and prediction model. |
| `downloadSourmashRocksDB`, `downloadSourmashTaxonomy` | GTDB RS226 references and matching taxonomy, shared by both optional Sourmash screens. |

## Workflow overview

![Phage isolate workflow: reads, assembled contigs and assigned hosts](figures/phage_isolates_overview.png)

[Zoomable overview (SVG)](figures/phage_isolates_overview.svg)

The diagrams show the biological steps, not individual Snakemake jobs. Blue is reads, teal is contigs, purple is hosts, amber is exclusions/review, and grey is unexplained reads/reports. Each major step has its own diagram below. Thresholds shown are the repository defaults; a run with config overrides can use different values. The saved SPAdes filter is **1,000 bp (1 kb)**, not 1,000 kb.

| Step | What happens |
| --- | --- |
| [1. Read QC](#1-read-qc) | fastp trims adapters/poly-G and low-quality sequence. FastQC/MultiQC, read counts, SuperDeduper statistics and Kraken provide QC/contamination evidence. Biological read removal is bypassed. |
| [2. Assembly and contig evidence](#2-assembly-and-contig-evidence) | BBnorm-normalised QC-passed reads enter SPAdes. Save contigs >=1 kb, assign persistent IDs, then collect geNomad, CheckV and terminal-repeat evidence. Assembly reads are not the purity-mapping denominator. |
| [3. Host characterisation](#3-host-characterisation) | CheckM, geNomad, CheckV and eligible BACPHLIP predictions describe assigned hosts and viral regions. Optional activity testing measures coverage enrichment. |
| [4. Retained and excluded contigs](#4-retained-and-excluded-contigs) | Full-read own-assembly mapping supplies depth. Apply the length/depth rule and standard assigned-host BLAST screen; record plasmid similarity as a report-only flag. Preserve retained, combined excluded and split excluded FASTAs. |
| [5. Dereplication and clustering](#5-dereplication-and-clustering) | Concatenate all saved contigs, exactly dereplicate with MMseqs, then cluster at independent >=95% ANI and >=85% target coverage. Extract representatives and search RefSeq; optional Sourmash/METAVR screens add evidence. |
| [6. Priority read accounting](#6-priority-read-accounting) | Assign each QC-passed read once through eight stages. Host-excluded contigs precede threshold-only excluded contigs, both for own and other-sample references. Preserve alternative/single-mate evidence and the final unexplained reads. |
| [7. Reports and review](#7-reports-and-review) | Sample decisions, contig/cluster metadata, host summaries and figures. Optional full-read Sourmash, RefSeq mapping and host identification remain separate comparisons. |

To regenerate the compact diagrams from the repository root:

```bash
conda activate Mosaic_workflow
python mosaic/docs/generate_isolate_diagrams.py
```

This documentation-only generator uses matplotlib and PyYAML from the workflow environment. It reads `config.yaml` for configurable thresholds and saves PNG plus editable-text SVG files under `mosaic/docs/figures/`. It does not call Snakemake, execute analysis, download databases or change project results. Layouts and shared alignment/clustering policies are explicit in the generator, so update those labels if the workflow policy changes.

See [clean-read Sourmash](sourmash_clean_reads.md), [cluster-representative Sourmash](sourmash_contig_catalogue.md), [contig identifiers](contig_identifiers.md) and [optional RefSeq mapping](refseq_mapping.md) for the branch-specific details.

## 1. Read QC

![Read QC: fastp, raw-read duplicate reporting and all QC-passed reads](figures/phage_isolates_01_read_qc.png)

[Zoomable read-QC diagram (SVG)](figures/phage_isolates_01_read_qc.svg)

With `isolates=True`, biological read removal is bypassed, including configured-contaminant and eukaryotic removal. Kraken classifies the QC-passed reads once for reporting; no post-decontamination Kraken job is required. All fastp-passed paired reads and orphans enter mapping, never the 2M subset. BBnorm normalisation is for assembly, not purity mapping. SuperDeduper measures PCR duplicates on the raw pairs for QC reporting; the current assembly path does not consume its deduplicated output, and neither fastp deduplication nor a separate PCR-duplicate filter is enabled.

Fewer than 20,000 surviving fastp pairs produces a report-only QC warning (`isolate_min_fastp_pairs`). The trimmed R1 count counts pairs once; R2 and orphan reads are not added to it. Exactly 20,000 pairs passes this check. Missing measurements and non-paired libraries remain unassessed. This warning does not change contig retention or the recovery/purity decision.

## 2. Assembly and contig evidence

![Assembly: BBnorm, SPAdes, saved contigs and annotation evidence](figures/phage_isolates_02_assembly_evidence.png)

[Zoomable assembly diagram (SVG)](figures/phage_isolates_02_assembly_evidence.svg)

Existing saved SPAdes assemblies retain their filenames and >=1 kb filtering. Persistent renamed IDs identify contigs throughout the analysis; neighbouring `.ids.tsv` sidecars preserve old assembler names for provenance only. geNomad, CheckV and terminal-repeat checks run on the saved assembly. Their results are evidence, not mandatory viral selection gates. Six biological categories are not imposed. RNA assembly, virome-positive filtering, VIBRANT, VirSorter, Pharokka, Phynteny and legacy AAI analysis are not required by the isolate target.

New DNA assembly runs preserve the SPAdes log at `03_CONTIGS/LOGS/{sample}_spades.{sampling}.log` before removing the working directory. The log includes available insert-size estimates and assembly diagnostics; failed runs retain the captured console log at the same path. Previously deleted logs cannot be recovered, and this does not retroactively populate logs for existing assemblies. Snakemake may schedule existing assemblies again because their rule code changed; inspect a dry-run before executing.

These saved sequences come from SPAdes `scaffolds.fasta`: one scaffold can contain unresolved joins, so one retained sequence is not necessarily one gap-free genome. The final contig metadata records `assembly_gapped`, `assembly_n_bases`, `assembly_gap_count`, `assembly_longest_gap_bp` and `assembly_gap_coordinates` (1-based inclusive N-run intervals). Ns are unknown sequence, not a measured gap length; SPAdes can use 100 Ns as a repeat-resolution placeholder. Gap flags are report-only and never change retention, clustering or PASS/REVIEW/FAIL. A sample can PASS while retaining a gapped scaffold. Sequences are not edited or filled from another sample.

MOSAIC-written FASTAs use one header line and one complete sequence line per record. This includes saved assemblies, retained/excluded sets, pooled mapping references, host viral/enriched-region exports and network/VIRIDIC inputs. Sequence content and identifiers are unchanged; third-party working files can retain their native format.

The DNA assembly log-scale plots show ordinary numerical ticks (`1`, `10`, `100`, `1,000`, `10,000`, …), with a lower limit of 1 and an upper power-of-ten limit large enough for the plotted values. Zero bars remain zero in the tables; all-zero panels use the existing linear-scale display.

## 3. Host characterisation

![Assigned hosts: classification, viral regions and optional activity testing](figures/phage_isolates_03_hosts.png)

[Zoomable host diagram (SVG)](figures/phage_isolates_03_hosts.svg)

Host FASTAs require explicit sample-to-host assignments. CheckM describes host assembly quality; geNomad supplies predicted embedded prophages and whole-contig viral candidates, kept as separate evidence scopes. CheckV and terminal-repeat evidence describe those viral regions. BACPHLIP is limited to CheckV-complete Caudoviricetes genomes.

The assigned host supplies the masked chromosome and viral-region references used below. The optional all-host comparison is report-only. [Host activity testing](#host-prophage-activity) uses the assigned unmasked host mapping, compares each embedded prophage with the pooled non-viral, non-plasmid background of that host assembly, and does not change retention. It is enabled by default for isolate runs with hosts.

## 4. Retained and excluded contigs

![Selection: coverage support and assigned-host matches, with retained and excluded FASTAs preserved](figures/phage_isolates_04_selection.png)

[Zoomable selection diagram (SVG)](figures/phage_isolates_04_selection.svg)

Full-read own-assembly mapping supplies contig depth. Retention requires:

```text
depth >= isolate_min_depth
AND (length >= isolate_min_length_bp OR depth >= isolate_short_min_depth)
AND no qualifying standard assigned-host or predicted viral-region match
```

Defaults are 5x, 4,000 bp and 10x. The standard host BLAST screen is unchanged: >=90% identity AND >=90% query coverage, controlled by `isolate_host_min_identity` and `isolate_host_min_query_coverage`. Only the sample's assigned host drives exclusion. Supported host-prophage-associated genomes are excluded from the retained set but preserved and flagged. Other-host matches remain available in metadata and the all-host BLAST table/heatmap; they never automatically exclude or relabel a contig.

### Report-only host-plasmid flag

An isolate plasmid can be longer than its related contig in the host assembly, so it can fail the standard query-coverage requirement even when most of the host plasmid matches exactly. The following comparison records that association without adding an exclusion route or lowering the general 90/90 screen:

```text
isolate contig AND matching assigned-host contig have saved geNomad plasmid-positive calls
AND identity >= isolate_host_plasmid_min_identity                  # default 95%
AND isolate-contig coverage >= isolate_host_plasmid_min_query_coverage   # default 75%
AND host-contig coverage >= isolate_host_plasmid_min_subject_coverage    # default 90%
```

All three alignment thresholds must pass for the same isolate/host-contig pair to establish the flag. Query coverage uses the original isolate-contig length; subject coverage uses the matched host-contig length. Inclusive coordinate unions prevent overlapping HSPs from increasing coverage by simple double-counting. Saved plasmid-summary membership supplies the calls; no extra geNomad score cutoff is invented. A small matching region in a much longer isolate contig does not establish the flag because the 75% query requirement still applies.

For example, a 2,537-bp isolate plasmid matching a 2,019-bp host plasmid at 100% identity, 78.4% query coverage and 94.7% subject coverage fails the standard host screen but receives `host_plasmid_match=True`. It remains retained if its length/depth support passes and it has no qualifying standard host or predicted viral-region match. The plasmid flag does not create an exclusion reason. `host_match` retains its original 90/90 definition. Matches below the flag thresholds remain visible as evidence; a negative flag does not establish absence from the host.

Selection reuses each sample's existing geNomad directory and only its assigned host's directory, alongside the existing BLAST results. No new classifier run, BLAST format, database or environment is introduced. Missing subject coordinates in historical eleven-column BLAST files cannot establish this flag; regenerate those BLAST outputs when migrating. geNomad viral positivity is still not required to retain a contig, and the standard host screen does not require a plasmid call. This is host-association evidence, not proof of contamination or proof that two assemblies are identical.

The plasmid flag does not choose a FASTA or mapping stage. Flagged contigs that pass the normal selection enter the retained FASTA and retained-reference mappings; flagged contigs failing length/depth alone enter `excluded_threshold`. A plasmid that also passes the standard 90/90 host screen remains host-excluded under that existing rule. The eight-stage order, full-read denominator, colours, combined excluded totals and recovery/purity warning thresholds are unchanged. Rerunning the isolate target refreshes selection, affected reference pools, mappings and reports; existing project results are not edited in place by a source-code update.

The isolate summary shows retained, host-excluded and threshold-excluded contigs in its overview. An additional selection-breakdown figure preserves the overlap: retained, excluded by length/depth thresholds only, excluded by assigned-host matches only, and excluded by both. These four groups count every saved original contig once. Host exclusion uses the recorded standard and prophage/viral-region matches, not the report-only plasmid flag; the plot reads those decisions without applying another filter. The notebook and HTML display the figure with up to 50 samples per panel and consistent bar widths. `FIGURES_AND_TABLES/08_phage_isolates_contig_selection.tot.tsv` saves all per-sample counts plus excluded and assembled totals; matching PNG/SVG files contain the figure.

Isolate host and RefSeq/METAVR BLAST outputs append `sstart send bitscore btop` to the existing eleven columns. Query/reference breadth uses merged inclusive coordinate intervals, without double-counting overlapping alignments. Host identity is calculated from the same selected query portions, using the highest-scoring alignment where HSPs overlap and counting mismatches/gaps from BTOP. Each query/reference pair is assessed separately. Historical eleven-column files remain readable with identity explicitly labelled `legacy HSP identity estimate` and unavailable reference breadth left missing; regenerate BLAST for exact traceback-based identity. Other annotation datasets retain their existing BLAST output format.

## 5. Dereplication and clustering

![All-contig relatedness: exact representatives, independent 95/85 clustering and provenance](figures/phage_isolates_05_clustering.png)

[Zoomable clustering diagram (SVG)](figures/phage_isolates_05_clustering.svg)

All saved contigs retain the hierarchy `original_contig -> exact_rep -> mosaic_cluster -> cluster_rep`. Existing concatenation, exact MMseqs dereplication, provenance, per-assembly geNomad/CheckV and CheckV ANI calculation are reused, including DNA-only sharing links. Both retained and excluded sequences enter this global catalogue, without a top-N restriction. MMseqs uses 100% identity and 100% target coverage.

Every clustering command uses the original installed CheckV tool:

```bash
aniclust --fna ... --ani ... --out ... \
  --min_ani 95 --min_tcov 85 --min_qcov 0
```

ANI and target coverage pass independently. The local modified product-based `scripts/aniclust_checkv.py` is not invoked or edited. Optional Sourmash screens representatives >=5 kb but never removes shorter clusters from the catalogue.

### Fragment-candidate flags

The catalogue also screens retained contigs for containment in longer retained contigs, including members from other samples and other clusters. It reuses the existing ANI output; no new BLAST, clustering or mapping run is needed.

```yaml
isolate_fragment_min_identity: 95
isolate_fragment_min_short_coverage: 95
isolate_fragment_max_length_ratio: 0.80
```

`fragment_candidate=True` means that a qualifying match has at least the configured identity and coverage of the shorter sequence, and the shorter/longer length ratio is at most the configured maximum. There is no absolute length cutoff. These thresholds are separate from the unchanged 95/85 clustering policy. Clustering is representative-based: a fragment can match a cluster member without passing the identity threshold against that cluster's representative.

The best match is recorded in `all_contig_metadata.tsv` as `contained_in_contig`, `contained_in_sample`, `contained_in_exact_rep`, `contained_in_cluster`, `contained_in_length_bp`, `containment_identity_percent`, `containment_short_coverage_percent`, `containment_length_ratio` and `containment_evidence_scope`. Matches are ranked by identity, shorter-sequence coverage, longer target length and finally target ID. Every qualifying original-contig pair is preserved in `fragment_candidate_matches.tsv`, with its `match_rank`; exact duplicates remain attributable to their original samples. `cluster_metadata.tsv` adds `number_fragment_candidates`, counting flagged retained original contigs rather than independent genomes.

ANI percentages are expanded only to verified full-length identical copies of exact representatives, including reverse-complement copies. Existing MMseqs groups can also contain shorter perfect matches: these originals are checked directly for exact sequence containment within longer retained members of their group, on either strand, and labelled with that evidence scope. This avoids assigning a longer representative's alignment percentages to a shorter member. Excluded contigs are not screened; an unflagged retained sequence means no qualifying evidence in these existing comparisons, not a complete genome.

This is a report-only candidate flag, not proof of incompleteness. Complete plasmids or other mobile elements can also occur within longer sequences, so CheckV, geNomad plasmid and terminal-repeat evidence stays visible for review. Retention, cluster membership, read-accounting fractions and PASS/REVIEW/FAIL decisions are unchanged. The flag does not remove any catalogue sequence. The optional within-cluster VIRIDIC report omits flagged fragments from its own comparison inputs only. The isolate summary displays the candidates and supporting matches, and the main sample table adds `retained_fragment_candidates` and `dominant_contig_fragment_candidate`.

## 6. Priority read accounting

![Read accounting: eight priority mapping stages and final unexplained reads](figures/phage_isolates_06_read_accounting.png)

[Zoomable read-accounting diagram (SVG)](figures/phage_isolates_06_read_accounting.svg)

```text
own retained -> masked host chromosome -> host viral regions
             -> own host-excluded -> own threshold-excluded
             -> other samples' retained -> other samples' host-excluded
             -> other samples' threshold-excluded -> unexplained
```

Only reads without a passing alignment advance. Every stage uses >=95% identity and >=85% aligned read length with deterministic Bowtie2 mapping. The first stage reuses the own-assembly filtered BAM, which is retained in isolate mode so a revised host/selection does not require repeating that alignment.

The shared `mapping_max_insert_bp` setting defaults to 1,000 bp and supplies Bowtie2 `-X` to every paired-end mapping rule in MOSAIC, including assembly, host, catalogue, RefSeq, microbial, polishing and microdiversity mappings. It describes the maximum concordant paired-end fragment span, not read length or a hard filter rejecting every longer pair: discordant/single-mate alignments remain possible. Override it with `--config mapping_max_insert_bp=1500`. Single-end and long-read mappings are unchanged. Changing it requires regenerating affected alignments; touching old BAMs would retain the previous mapping behavior.

The `06_other_retained` alignment remains one stage. A separate lightweight reporting rule reads its saved accepted-reference lists and the existing MOSAIC catalogue, and partitions every assigned read into exactly one group:

| Reporting group | Accepted-reference evidence |
| --- | --- |
| `06_other_retained_same_cluster` | All accepted references belong to clusters containing this sample's retained contigs. |
| `06_other_retained_mixed_clusters` | Accepted references include both those clusters and other clusters. |
| `06_other_retained_different_cluster` | No accepted reference belongs to this sample's retained-contig clusters. |

All accepted references are used, not just the mapper's selected placement. Samples with no retained contigs have no same-cluster group. This adds no alignments and does not change downstream residual reads, reference selection or coverage. Same-cluster recovery can reflect scaffold gaps or related sequence variants; it is not evidence of another phage. Mixed/different-cluster signal is not automatically contamination. The three groups sum to the preserved `06_other_retained` total; the `05_other_retained` compatibility alias remains unchanged. Tables and plots show the split instead of stacking both it and its aggregate. Their entering-read denominator is the original stage-06 input, not a new sequential flow.

Read-level classifications are saved in `06_MAPPING/ISOLATES/other_retained_cluster_reads.tot.tsv.gz`; the adjacent `other_retained_cluster_summary.tot.tsv` and `other_retained_cluster_reference_reads.tot.tsv` supply both notebooks. These reporting files can be regenerated from existing mappings without rerunning Bowtie2.

Excluded contigs are partitioned into two non-overlapping FASTAs before mapping. `excluded_host` contains qualifying standard or predicted viral-region matches to the source sample's assigned host, including contigs that also fail length/depth. `excluded_threshold` contains excluded contigs with no qualifying standard or viral-region match. The report-only plasmid flag does not choose either group. A separate `SELECTION/{sample}.tot.exclusion_groups.tsv` adds `exclusion_group` and `primary_exclusion_reason` while preserving every selection flag and reason. The catalogue reads this enriched selection table; the original selection TSV and combined excluded FASTA remain available.

The primary reason is host-prophage association first, then the standard host match, minimum-depth failure or short-contig depth failure. The standard host screen does not require a geNomad call. Plasmid similarity is preserved in metadata rather than used as an exclusion reason. Contigs failing both standard host and length/depth checks go into the host-excluded mapping set; the overlap remains visible in metadata and the selection-breakdown figure.

The priority-assignment plots and diagram use matching colour families: own retained contigs use blue `#3A72ED`, and other-retained same/mixed/different-cluster evidence uses light blue `#93B4F7`, blue-grey `#9297CE` and dark blue `#244B91`. Host chromosomes/viral regions retain purple `#7C3AED` / lavender `#B69AEF`. Own/other host-excluded contigs use dark/pale purple `#4C1D95` / `#DDD1F3`, threshold-excluded contigs use dark/light orange `#B45309` / `#E6B07C`, and unexplained reads remain grey `#8A8A8A`. These colours describe reference groups, not PASS/REVIEW/FAIL status; the conceptual diagram keeps stage 06 together, with its report-only split identified below.

This is genuinely sequential mapping, not a retrospective label: a read matching both excluded sets is assigned to host-excluded first, and only its residual reads enter threshold-excluded. The split can change selected references, pairing/coverage details and downstream assignments compared with the previous combined-reference mapping. Empty reference sets pass all entering reads to the next stage without invoking Bowtie2.

Reads, not alignment records or pairs, are counted. Broken pairs become orphans with original mate identity recorded. Single-mate, discordant and non-unique evidence remains explicit. Passing alternative references are retained, but one deterministic placement per read contributes to stage coverage. Both full-QC and entering-stage percentages are reported.

Host FASTAs are discovered under `HOST/`. When host FASTAs exist and `isolates=True`, `host_mapping_file.tsv` is mandatory at the project root or under `HOST/`, with `sample` and `host` columns. Every analysed sample must have one unambiguous assignment to an available host FASTA. Missing assignments, unavailable assigned references and conflicting rows/files stop DAG construction before jobs run. If both metadata locations exist, their assignments must agree. Sample names are not used to guess the host.

Residual host stages always use the assigned host; there is no all-host fallback. Runs without any host FASTAs remain possible and record host origin as not assessed, not as evidence of purity. Only reported viral coordinates are masked in isolate mode, without extra flanks. Whole-contig viral candidates are moved out of the chromosome reference into the host-viral reference.

The startup summary is printed only by the main Snakemake process. Discovered Illumina, Nanopore, PacBio and host lists use natural numeric order (`S1`, `S2`, ..., `S10`). One grouped summary shows paths, sample/host counts, wrapped names and effective enabled options; key samples do not repeat the full Illumina list when identical. Contig extraction, reference pooling, residual-read extraction and BACPHLIP input selection use explicit inline Python in `shell:` blocks, avoiding additional Snakemake workers for these lightweight steps. They reuse `env5.yaml` for pandas/Biopython and `env1_mapping.yaml` for residual-read extraction; no new environment is required. Other existing `run:` or shadow rules elsewhere can still start normal Snakemake workers, without repeating the startup summary or restarting the full workflow.

Cross-sample references reuse pooled exact representatives and membership tables. Self-only representatives are ineligible; representatives shared with other samples are eligible. Exclusion groups follow each original member's source sample, not the query sample's host or the cluster representative's classification. Membership tables record `source_expected_host`. The same exact representative can occur in both excluded pools if its original members have different selection outcomes; host-excluded mapping takes priority. A cross-sample match is sequence sharing, not proof of its source or contamination.

The main read-accounting TSV has eight mapping stages plus `09_unexplained`, not aggregate rows. The separate `07_mapping_statistics_isolates_aggregates_tot.tsv` preserves combined own/other excluded totals and the legacy other-retained alias. Reads and percentages are added across the corresponding stages; reference counts use the union of distinct positively assigned reference IDs. These overlapping compatibility totals must not be added to the nine partition fractions. Both mapping reports display the detailed groups and retain the combined totals in tables.

`host_identification_test=True` separately enables full-read mappings against every masked/unmasked host. Breadth/depth/count comparisons are not summed into the purity fractions and never automatically relabel samples. Minimum evidence defaults to 100 reads and 1% masked-host breadth; equal best-host measurements remain ambiguous.

## 7. Reports and review

![Reports: PASS, REVIEW, FAIL and separate QC warnings](figures/phage_isolates_07_reports.png)

[Zoomable reports diagram (SVG)](figures/phage_isolates_07_reports.svg)

Start with the sample summary, then use contig/cluster tables and host reports to review the evidence. `REVIEW` flags a concern without automatically rejecting a recovered genome. Read-QC warnings remain separate. The [full decision definitions](#recoverypurity-decisions-and-separate-qc-warnings) distinguish review triggers from unusable data or invalid accounting.

| Location | Outputs |
| --- | --- |
| `03_CONTIGS/ISOLATES/` | Retained/excluded FASTAs and shared reference pools |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/` | `all_contig_metadata.tsv`, `cluster_metadata.tsv`, `sample_host_assignments.tsv`, `fragment_candidate_matches.tsv` |
| `06_MAPPING/ISOLATES/{sample}/` | Stage summaries, read-level evidence, counts, coverage/BAMs and unexplained FASTQs |
| `NOTEBOOKS/` | Executed notebooks and matching HTML exports of public reports; the internal catalogue notebook remains available as `.ipynb` only |
| `FIGURES_AND_TABLES/` | Tables/HTML and purity, contig-count, relatedness, host-match and depth plots |

### Which outputs to open first

All paths below are relative to the project results directory. `tot` identifies the full-data analysis, not a 2M mapping subset.

Public reports also have a matching HTML beside the executed notebook, for example `NOTEBOOKS/08_phage_isolates_summary.tot.html`. These exports retain code, tables and displayed plots with embedded images; conversion does not rerun the analysis. Existing HTML summaries in `FIGURES_AND_TABLES/` are unchanged. Activity/discovery exports are requested only when those analyses are enabled; individual selection notebooks and the internal catalogue parser are not exported.

For an older completed project, backfill the existing public report notebooks from `mosaic/`:

```bash
snakemake -p runNotebookHTML \
  --config input_dir=/path/to/project/00_RAW_DATA metagenome=False isolates=True \
  --allowed-rules runNotebookHTML notebook_to_html -j 4
```

The allowed-rule list prevents upstream analyses from running. Add `map_to_RefSeq=True` if the project has a RefSeq detection notebook to export. Jupyter from the existing `Mosaic_workflow` environment handles conversion.

| File or directory | What to use it for |
| --- | --- |
| `FIGURES_AND_TABLES/08_phage_isolates_summary.tot.html` | Start here for the combined figures, sample decisions, evidence and review notes. |
| `FIGURES_AND_TABLES/08_phage_isolates_summary.tot.csv` | Main machine-readable table: one row per sample, with recovery/purity decision, QC, mapping fractions and dominant-contig evidence. |
| `NOTEBOOKS/08_phage_isolates_summary.tot.ipynb` | Executed report with explicit analysis/plotting code and displayed results. |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/all_contig_metadata.tsv` | One row per saved original contig, including retained/excluded status, host evidence, annotation and cluster membership. |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/cluster_metadata.tsv` | One row per global cluster: member/sample counts, annotations, retention counts and optional Sourmash evidence. |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates.tot/fragment_candidate_matches.tsv` | All qualifying retained shorter/longer original-contig pairs, alignment or exact-containment evidence, original cluster IDs and match ranks. The best match is also in `all_contig_metadata.tsv`. |
| `03_CONTIGS/ALL_ASSEMBLED/phage_isolates_cluster_representatives.tot.fasta` | Representative sequence from every global all-contig cluster, not only viral-positive or retained clusters. |
| `03_CONTIGS/ISOLATES/{sample}_{retained,excluded}.tot.fasta` | Review the sequences kept for own-retained mapping versus those excluded by the length/depth/assigned-host policy. Excluded does not mean deleted. |
| `03_CONTIGS/ISOLATES/{sample}_{excluded_host,excluded_threshold}.tot.fasta` | Non-overlapping host-first exclusion sets used by the new sequential mappings. |
| `03_CONTIGS/ISOLATES/SELECTION/{sample}.tot.tsv` | Per-sample selection decision and exclusion reason before global reporting. |
| `03_CONTIGS/ISOLATES/SELECTION/{sample}.tot.exclusion_groups.tsv` | The same original-contig evidence, plus the mapping group and host-first primary exclusion reason. |
| `FIGURES_AND_TABLES/08_phage_isolates_read_accounting.tot.tsv` | Sample/stage read totals, fractions and mapped-reference counts; includes the terminal unexplained group. |
| `FIGURES_AND_TABLES/07_mapping_statistics_isolates_aggregates_tot.tsv` | Combined own/other excluded totals and the legacy other-retained alias; separate from the mutually exclusive read partition. |
| `06_MAPPING/ISOLATES/{sample}/` | Detailed stage summaries, reference-level counts, read evidence, BAM/coverage files and final unexplained FASTQs. Some intermediates are temporary and may already be cleaned. |
| `FIGURES_AND_TABLES/08_phage_isolates_sample_clusters.tot.tsv` | Sample-by-cluster membership, retained/excluded counts, representative breadth and annotation scope. |
| `FIGURES_AND_TABLES/08_phage_isolates_cluster_annotations.tot.tsv` | Representative CheckV completeness and non-exclusive member-evidence counts, in the same row order as the hierarchically clustered coverage heatmap. |
| `FIGURES_AND_TABLES/08_hosts_summary.tot.html` and `.tsv` | Host CheckM quality and numbers/types of predicted viral regions. The executed notebook is `NOTEBOOKS/08_hosts_summary.tot.ipynb`. |
| `FIGURES_AND_TABLES/08_host_viral_regions.tot.tsv` | Host-region coordinates, embedded versus whole-contig scope, CheckV, terminal-repeat and eligible BACPHLIP evidence. |
| `FIGURES_AND_TABLES/08_host_prophage_activity.tot.tsv` | Optional assigned-host sample/region coverage-enrichment calls. The matching notebook and coverage profiles explain the measurements. |
| `NOTEBOOKS/08_host_enriched_regions.tot.ipynb` | Additional exploratory coverage profiles: geNomad regions in orange, unpredicted enriched host-DNA blocks in red, with 200 kb context on each side. Coordinates and per-sample summaries are in `FIGURES_AND_TABLES/08_host_enriched_regions.tot.tsv` and `08_host_enriched_regions_summary.tot.tsv`. |
| `HOST/enriched_regions/host_enriched_regions.tot.fasta` and `.tsv` | Unique enriched host loci, merging overlapping calls across samples. The table lists union coordinates, original intervals and supporting samples; plotting flanks are not exported. |
| `HOST/enriched_regions/host_enriched_regions_metadata.tot.tsv` | One row per exported locus, with provenance, prophage proximity, Phold gene-function evidence, geNomad, CheckV, TaxMyPhage and RefSeq BLAST; optional METAVR best-hit metadata. |
| `FIGURES_AND_TABLES/08_phage_isolates_sourmash_clean_reads.tot.tsv` and `08_phage_isolates_sourmash_clean_taxonomy.tot.tsv` | Optional per-sample and pooled read-profile evidence and taxonomy. Raw Sourmash outputs are under `07_ANNOTATION/SOURMASH_CLEAN/`. |
| `NOTEBOOKS/01_sourmash_clean_reads.tot.ipynb` and `FIGURES_AND_TABLES/01_sourmash_clean_reads.tot.html` | Shared clean-read GTDB report also available to virome/QC runs. The isolate summary reuses its tables and figures while keeping the existing isolate output names. |
| `06_MAPPING/REFSEQ_VIRAL/` | Independent RefSeq count/RPKM/breadth/depth matrices and detection tables, only with `map_to_RefSeq=True`. |

### Main sample-table columns

The sample CSV places the ten detailed mapping percentages together in priority order,
followed by `unexplained_percent`, then matching read-count and mapped-reference-count
blocks. Compatibility totals are grouped after these blocks, not mixed into the detailed
partition. Reporting TSVs use `qc_passed_reads`, `percent_qc_passed` and
`remaining_percent_qc_passed`; existing per-stage caches remain readable without remapping.

| Column(s) | Meaning |
| --- | --- |
| `sample`, `expected_host`, `host_assessment` | Library ID, explicit host assignment and whether an assigned host reference was assessed. |
| `decision`, `warnings`, `failure_reasons` | Recovery/purity verdict and measured reasons. `FAIL` takes precedence over `REVIEW`. |
| `recovery_state` | `single_retained_contig`, `fragmented_or_mixed`, `no_retained_genome`, or `supported_host_associated_genome_excluded`. These describe recovery, not a definitive biological diagnosis. |
| `qc_status`, `qc_warnings` | Independent read-QC assessment. A PCR-duplicate warning does not automatically change the recovery/purity decision. |
| `fastp_pairs`, `fastp_pairs_status`, `fastp_pairs_minimum` | Surviving paired reads counted once using trimmed R1, the report-only PASS/WARN/INFO assessment and its configured threshold. |
| `qc_passed_reads` | All QC-passed paired and orphan reads entering accounting, not raw/pre-QC reads. Each mate is one read, not one pair. |
| `assembled_contigs`, `retained_contigs`, `excluded_contigs` | Saved >=1 kb contig counts before and after the coverage/assigned-host decision. |
| `excluded_host_contigs`, `excluded_threshold_contigs` | Host-first exclusion counts; their sum is `excluded_contigs`. A contig failing both checks belongs to the host group. |
| `host_plasmid_associated_contigs` | Report-only count of original contigs with a plasmid-positive assigned-host match at 95% identity / 75% query / 90% subject coverage; includes retained and excluded contigs. |
| `retained_contig_lengths_bp` | Semicolon-separated lengths of original retained contigs, largest first; duplicate lengths are kept and no retained contigs is `none`. |
| `assembly_gap_flag`, `retained_gapped_scaffolds`, `retained_assembly_n_bases`, `retained_assembly_gap_count` | Report-only evidence of unresolved sequence in retained scaffolds. Gap flags do not change PASS/REVIEW/FAIL. Individual coordinates and maximum N-run length are in the contig metadata. |
| `all_contig_clusters`, `retained_contig_clusters` | Distinct global cluster representatives represented by all saved contigs or retained contigs in this sample. |
| `{stage}_reads`, `{stage}_percent` | Reads assigned at one of the eight stages `01_own_retained` through `08_other_threshold_excluded`, and percent of `qc_passed_reads`. Stage 06 is split into three reporting groups; the ten detailed fractions plus `unexplained_percent` sum to 100%. |
| `04_own_excluded_{reads,percent,mapped_references}` | Compatibility total across `04_own_host_excluded` and `05_own_threshold_excluded`; not an additional mapping stage. |
| `05_other_retained_{reads,percent,mapped_references}` | Compatibility alias of the new `06_other_retained` stage. |
| `06_other_retained_{same_cluster,mixed_clusters,different_cluster}_{reads,percent,mapped_references}` | Mutually exclusive reporting groups within stage 06, using all accepted references. Their read counts and percentages sum to the original other-retained total; distinct-reference counts are not additive. |
| `06_other_excluded_{reads,percent,mapped_references}` | Compatibility total across `07_other_host_excluded` and `08_other_threshold_excluded`; reference counts are deduplicated across the two stages. |
| `{stage}_mapped_references` | Distinct reference sequences receiving assigned reads at that stage; the heatmap cell's `n`. Not an additive genome count. |
| `retained_mapping_percent` | Alias of `01_own_retained_percent`, used for the recovery warning. |
| `unexplained_reads`, `unexplained_percent` | Reads with no accepted match after all eight stages. These are not automatically contaminants. |
| `dominant_contig`, `dominant_contig_reads`, `dominant_contig_length_bp` | Most read-supported retained contig, its stage-01 assigned reads and sequence length. |
| `dominant_percent_full_qc`, `dominant_percent_retained` | Dominant-contig fraction of all QC-passed reads versus only own-retained assigned reads. The denominators differ. |
| `dominant_contig_genomad_*`, `dominant_contig_checkv_*` | Evidence for that original retained contig, not annotations borrowed from its cluster representative. |
| `retained_genomad_virus_calls`, `retained_genomad_plasmid_calls` | Non-exclusive counts of retained original contigs with each geNomad call. |
| `supported_host_prophage_contigs` | Coverage-supported original contigs excluded because they match their assigned host's predicted viral regions. |
| `host_prophage_activity_assessment`, `active_host_prophage_count`, `active_host_prophage_regions`, `active_host_prophage_details` | Assigned-host activity evidence, including region IDs/depth/breadth. Unassessed counts are missing, not zero. |
| `host_identification_status`, `host_identification_best_host` | Optional all-host comparison; never an automatic reassignment. The best-host column may be absent when no comparison supplies a best hit. |
| `read_accounting_conserved`, `read_accounting_valid` | Total-count conservation versus the stronger check of valid counts, consistent denominators and stage transitions. |

In the per-contig table, `original_id` is the saved, renamed assembly ID before dereplication, not the original SPAdes header. `exact_rep`, `mosaic_cluster` and `cluster_rep` record the clustering hierarchy. `coverage_supported`, `host_match`, `host_prophage_match`, `retained`, `contig_set` and `exclusion_reason` explain selection. `host_plasmid_match` is a separate report-only association flag. `best_host_blast_contig` and `host_blast_subject_coverage_percent` retain the standard best-match subject evidence. The separate plasmid fields record `host_plasmid_assessment`, `host_plasmid_query_call`, both geNomad scores (`host_plasmid_query_score` / `host_plasmid_host_score`), `best_host_plasmid_reference`, `host_plasmid_contig`, and identity/query/subject percentages under `host_plasmid_*_percent`. The best plasmid candidate need not be the standard best BLAST hit. Neither this flag nor `other_host_*` evidence drives exclusion. `genomad_*`, `CheckV_*`, terminal-repeat evidence and optional `sourmash_*` fields describe their recorded evidence scopes. The old assembler header can be looked up in the neighbouring `.ids.tsv` sidecar; downstream tools do not need that table.

In the cluster table, `number_exact_representatives`, `number_original_contigs` and `number_samples` use different units. `number_retained`/`number_excluded` count original contigs. `number_genomad_virus`, `number_genomad_plasmid` and their fractions are non-exclusive original-member evidence. `best_CheckV_quality`/`maximum_CheckV_completeness` are the best member measurements, not proof that every member or the representative has that quality. `number_host_matching` retains the standard host-match count; `number_host_plasmid_matching` counts the report-only plasmid-association flags, irrespective of retention. These counts can overlap and must not be added as disjoint exclusion totals. `number_phage_plasmid_candidate` retains the underlying dual-biology evidence count. Use the per-contig table to identify the contributing members.

Text annotations without evidence use `not reported`; missing numeric values remain missing. Neither means that the sequence was proven non-viral or host-free.

### Figures and reporting dependencies

Plots are explicit notebook cells using existing MOSAIC conventions. The cluster-coverage heatmap includes all clusters. Per-contig raw RPKM and depth measurements remain in the metadata, but the cluster RPKM heatmap and dominant-contig depth heatmap are not generated or embedded in the report.

Catalogue generation uses `notebooks/08_isolate_contig_catalogue.py.ipynb`; the final report uses `notebooks/08_phage_isolates_summary.py.ipynb`. Existing executed notebook paths and main table names are unchanged. Report/plot edits therefore invalidate only the report, not catalogue selection, retained/excluded FASTAs or the eight-stage read accounting. Changes to the shared catalogue/selection notebook source or actual selection inputs/parameters still invalidate the relevant biological results. Global annotation tables are no longer mapping dependencies.

Bowtie indexes, non-filtered assembly-mapping intermediates, residual FASTQs and Sourmash sketches/gather intermediates remain temporary. The isolate own-assembly filtered BAM and stage read-evidence tables remain available for reuse. Missing temporary files alone do not require rerunning an up-to-date consumer. Residual FASTQs are reconstructed from the full QC-passed reads and cumulative assigned-read evidence by `extract_isolate_remaining_reads`; their deletion does not require repeating earlier alignments. Report-only refreshes must not force the catalogue, disable provenance checking globally or touch biological outputs.

## Incremental assembly and host updates

`select_isolate_contigs` runs the catalogue notebook internally in sample-selection mode, using only that sample's assembly, own-assembly coverage and geNomad directory, plus its assigned host's BLAST, viral-region BLAST and geNomad directory. No pooled sample-annotation table is a selection dependency. It writes `03_CONTIGS/ISOLATES/SELECTION/{sample}.tot.tsv`, without saving individual notebooks under `NOTEBOOKS/ISOLATES/`. The global catalogue collects these decisions alongside clustering and annotation results; it no longer determines the per-sample mapping references. The same selection and BLAST parsers are reused in both notebook modes. Use the combined catalogue and final summary notebooks for interpretation; the per-sample TSVs retain the complete selection evidence.

Host-stage references and indexes are shared only among samples with the same assigned host, under `03_CONTIGS/ISOLATES/REFERENCES/HOST/{host}/`. Runs without hosts use empty `REFERENCES/UNASSIGNED/` references. Host assignments are tracked as per-sample rule parameters, not as an input dependency on a globally regenerated catalogue table. Per-sample host BLAST reuses `run_BLASTn_host` with `-subject`, avoiding concurrent writes to the same host BLAST database; existing all-host BLAST output paths remain available for reporting.

An improved assembly reruns its own coverage/selection and early accounting stages. A changed host reruns host-specific analysis and early accounting for its assigned samples. Changes to retained/excluded sets or exact dereplication can still require cross-sample stages 06/07/08 for every sample, because those shared references really change. Their pools depend on the selected FASTAs and existing exact-cluster membership, not global geNomad/CheckV/Sourmash reporting metadata. Recreating deleted residual FASTQs at this boundary uses saved read evidence, not earlier alignments. Global clustering, annotation and combined reports retain their genuine shared dependencies.

Existing projects need a one-time refresh to create the split exclusion FASTAs and reference pools and execute the new sequential stages. Snakemake may also refresh earlier accounting jobs because the shared rule implementation changed; these reuse the saved own-assembly BAM where applicable. This can schedule a large initial run; the narrower invalidation applies once outputs are current. Do not use `--touch` or disable code/input/parameter tracking to hide this migration. Assembly, classification and clustering settings remain unchanged; plasmid similarity is report-only as described above. Old six-stage BAMs, tables and pooled references are not deleted automatically, but the new reports read the eight-stage outputs rather than obsolete stage files. Combined excluded FASTAs, main report paths and legacy summary-total columns remain available.

## Reading reports and figures

`01_QC.tot.ipynb` displays one k-mer sampling endpoint per sample and writes `FIGURES_AND_TABLES/01_kmer_rarefaction_summary.tot.csv`. It includes the maximum reported read-pair count, the corresponding individual-read count (twice the pairs; no orphans), and known/new first-position k-mer percentages in that last block. The table covers all samples even when plots are limited. The isolate summary reuses it in the notebook, HTML and main sample CSV. Known means previously observed in this file, not a database match or genome completeness. BBCountUnique defaults omit the final incomplete block; header-only histograms have missing measurements, not zero counts. A single endpoint does not establish a plateau and does not change warnings or PASS/REVIEW/FAIL. Per-block percentages can fluctuate with read order, library composition and sequencing quality; the reported histogram and its generation parameters are unchanged.

The main sample-level table is `FIGURES_AND_TABLES/08_phage_isolates_summary.tot.csv`: `decision`, `recovery_state`, `warnings`, and `failure_reasons` explain recovery/purity; `qc_status` and `qc_warnings` are independent read-QC assessments. The mapping heatmap `08_phage_isolates_mapping_heatmap.tot.png/.svg` shows all eleven detailed reporting fractions with PASS/REVIEW/FAIL at the right. A thin column-colour strip identifies the same reference groups as the stacked plots; the cells retain the continuous `mako` scale. `REVIEW` labels also show the corresponding recovery/purity reasons from `warnings`. A square-root colour scale makes small fractions easier to distinguish; the numerical percentages retain their original denominator, with colourbar ticks at 0, 1, 5, 10, 25, 50 and 100%. Cells show percent of all QC-passed reads and `n`, the number of distinct reference sequences receiving assigned reads at that stage, not all available contigs or alternative placements. Host columns count scaffolds/viral regions, and other-sample columns count pooled exact representatives. Unexplained reads have no reference count. Counts are saved as `{stage}_mapped_references` in the main CSV and `mapped_references` in the read-accounting TSV; they are not additive biological genome counts. Stacked-bar figures remain available. QC, isolate, host-summary and optional RefSeq stacked bars have no white borders between segments. The retained/excluded contig-count overview now separates host-excluded from threshold-excluded originals; the additional selection figure still shows failures of both tests. The heatmap is displayed in the notebook and included in the HTML report.

Review notes give the measured percentages and configured limits, including the chromosome/host-viral split and retained-contig/cluster counts. The main CSV also reports the dominant **retained** contig's length, geNomad classification/evidence scope and CheckV quality/completeness. These are its own annotations, never inherited from the cluster representative. If there is no mapped retained contig, text fields say `not reported` and numeric annotations are empty. `retained_genomad_virus_calls` and `retained_genomad_plasmid_calls` count original retained contigs; the counts can overlap and are not additional filtering criteria.

When host activity is enabled, the summary reuses its existing TSV for the assigned host and adds `host_prophage_activity_assessment`, `active_host_prophage_count`, `active_host_prophage_regions` and `active_host_prophage_details` (region ID with mean depth and breadth). Only embedded regions called `active` are counted; whole-contig viral candidates are not treated as embedded prophages. Disabled, unavailable and partially assessed activity is labelled explicitly; an unassessed count is empty, not zero. A host with no predicted embedded prophages has a zero count and the status `no_predicted_embedded_prophages`. Activity is supporting evidence, not confirmed induction, and never changes retention or the sample decision. These fields and the decision notes appear in both the notebook and HTML report without new mapping or activity calculations.

`08_phage_isolates_sample_clusters.tot.tsv` reports each sample/cluster combination, member counts split into retained/excluded, representative origin/length/classification/scores/CheckV quality, own-member quality, host-prophage evidence and assigned retained-read percentages. Representative quality is never inherited as the quality of every member. The cluster-summary CSV additionally flags `mixed_genomad_calls`; this means a cluster has both virus and plasmid calls, not that every member is a phage-plasmid. The existing individual `phage_plasmid_candidate` flag remains separate.

The `08_phage_isolates_cluster_coverage.tot` heatmap shows representative versus sample, with original-contig counts in the cells. Coverage is the union of the existing clustering BLAST subject intervals, using exact-representative membership to recover duplicate origins. It is not summed member length, a completeness estimate or read-mapping breadth. Five annotation blocks show the representative's own CheckV completeness (purple), existing geNomad virus calls (blue), plasmid calls (green), assigned-host-prophage matches (red), and fragment candidates among retained members (amber). Counts are `n/N` original cluster members, except fragments use retained members only; excluded contigs are not screened for fragments. Gray `NA` means no reported CheckV estimate or no assessed/eligible members, not a negative call. `?` marks partial assessment. The new `08_phage_isolates_cluster_annotations.tot.tsv` saves these counts, assessment denominators and representative annotations. No new score threshold or selection criterion is applied.

All clusters appear in one coverage heatmap, with samples on the x-axis and representatives on the y-axis. Average-linkage hierarchical clustering orders both axes using Euclidean distance between aligned-breadth profiles; counts and annotation blocks follow the same ordering. This is a reporting order, not a change to biological cluster membership. Dark representative labels are bold only for the representative's own CheckV Complete assessment; `†` marks a representative that itself is flagged as a fragment candidate. Member evidence is non-exclusive and is never inherited by the representative or every member. Sample labels use dark green (`PASS`), amber (`WARNING`) or red (`FAIL`); `WARNING` displays the existing `REVIEW` decision without changing tables, cutoffs or the separate read-QC status. The figure scales with the full sample and cluster counts, without pagination. Saved matrix tables retain their cluster-by-sample orientation and use the plotted row/column order. Individual Complete-representative alignment panels are no longer generated. The separate within-cluster VIRIDIC comparisons described below use actual VIRIDIC similarity scores, not these alignment-coverage values.

## Within-cluster VIRIDIC comparisons

`isolate_viridic=True` is enabled by default only when `isolates=True`. It adds one annotated VIRIDIC heatmap per eligible existing MOSAIC cluster to `08_phage_isolates_summary.tot.ipynb` and its HTML report. There is no new user-facing notebook and no whole-catalogue all-versus-all run.

The selection uses retained contigs only, excludes existing `fragment_candidate=True` records from these comparisons, and requires at least two distinct remaining sequences per cluster. Full sequences and their reverse complements are hashed; sharing an MMseqs representative alone is not treated as exact sequence identity. Identical copies within a cluster are compared once, while `members.tsv` preserves every included original contig, sample, comparison ID and its own annotation evidence. There is no new minimum-length cutoff or viral-positive filter. A genuine small genome without a fragment flag remains eligible.

By default, every eligible retained cluster is included, regardless of CheckV quality. Set `isolate_viridic_complete_only=True` to additionally require the cluster representative itself to be CheckV Complete, or `isolate_viridic=False` to disable this branch. `isolate_viridic_threads` defaults to 8 and is passed to VIRIDIC as `ncor`. The launch rule also caps the container's R future workers at the allocated thread count, because the bundled similarity script otherwise uses all visible server cores despite `ncor`.

The checkpoint `select_isolate_viridic_clusters` reuses `env5.yaml` and writes `03_CONTIGS/ALL_ASSEMBLED/phage_isolates_viridic.tot/`: `clusters.tsv` records eligibility/skipping reasons, `members.tsv` records selected originals and collapsed comparison IDs, and `INPUTS/{cluster}.fasta` contains the distinct sequences submitted for each eligible cluster. `viridic_isolate_cluster` inherits the existing `viridic_relatives_phages` command. It runs VIRIDIC with `steps=sim_clust`; the stock heatmap is replaced by explicit plotting code in the existing isolate notebook.

Tool results are under `07_ANNOTATION/VIRIDIC_CLUSTERS/tot/{cluster}/`, including `04_VIRIDIC_out/sim_MA_genCol.csv` and VIRIDIC's own `clusters.csv`. Logs are under `07_ANNOTATION/VIRIDIC_CLUSTERS/tot/LOGS/`, with one benchmark per cluster. Plots are saved as PNG/SVG under `FIGURES_AND_TABLES/08_phage_isolates_viridic.tot/`, alongside `plots.tsv` linking each plot to its matrix and sequence/sample counts. Singletons and exact-duplicate-only groups have no plot but remain in the selection table.

Cells show VIRIDIC intergenomic similarity percentages, with darker colors for higher similarity. Both axes share one average-linkage hierarchical tree based on pairwise distance `100 − VIRIDIC similarity`, grouping related sequences together. Labels show the compared original contig's length and its own CheckV quality/completeness; Complete contigs are bold. Green labels and `[plasmid]` tags reuse that contig's existing geNomad plasmid call, without a new threshold. Unknown completeness says `not reported`. A representative's quality is never inherited by other genomes or by identical copies with different annotation evidence. VIRIDIC's internal species/genus grouping is saved but does not replace the existing MOSAIC 95/85 clusters. Retention, fragment flags, mapping, purity decisions and the optional Mash network are unchanged. These within-cluster plots do not add external references.

No new Conda environment is needed: selection and the launch rule reuse `env5.yaml`, while VIRIDIC itself runs in the existing Singularity bundle configured by `viridic_folder` (default `tools/viridic_v1.1`). Install Singularity >=3.5 on the server and place the VIRIDIC bundle there, or set the path to an existing installation; it is not automatically downloaded by MOSAIC. The launcher selects the first compatible Singularity on PATH, skipping legacy 2.x binaries that can shadow a working installation. It does not modify Conda or the server installation. See the [official VIRIDIC installation instructions](https://github.com/CristinaMoraru/VIRIDIC#readme). The legacy relatives rule and its original output location remain available.

## Optional interactive genome network

Set `isolate_network=True` to add `FIGURES_AND_TABLES/08_isolate_network.tot.html` to the isolate target. The separate `isolate_network` target builds just this branch from the existing isolate results:

```bash
snakemake --use-conda -p isolate_network \
  --config input_dir=/path/to/project/00_RAW_DATA isolates=True isolate_network=True \
  -j 32 --rerun-incomplete
```

No extra environment or Singularity installation is needed: reference extraction reuses `env5.yaml`, and Mash 2.3 is already in `env7.yaml`. This branch does not run VIRIDIC. The separate legacy VIRIDIC rule and its output paths are unchanged.

The network includes every retained original contig, without a geNomad-positive filter. Identical sequences or reverse complements are compared once, while all original IDs and metadata remain available. This exact-sequence check is independent of MMseqs membership, which can also include shorter contained sequences. RefSeq relatives are selected using the existing overlap-aware BLAST parser; by default up to three hits per cluster-representative query pass >=70% identity and >=50% query coverage. Those cutoffs select reference candidates, not final graph connections. A hit to the existing cluster representative is recorded as discovery evidence, not an annotation inherited by its retained members.

`isolate_network_references="RefSeq METAVR"` adds METAVR relatives and requests the existing METAVR BLAST rule. This explicit network setting works independently of `metavr_blast`, which controls the broader report. Set `isolate_network_references=""` for own contigs only. RefSeq sequences are extracted from the configured FASTA; METAVR sequences are retrieved from the configured BLAST database. Missing sequences are listed rather than represented by invented nodes or scores.

Mash compares all distinct selected sequences, including reference-to-reference pairs. Defaults are `isolate_network_kmer_size=15`, `isolate_network_sketch_size=25000`, `isolate_network_seed=42` and `isolate_network_threads=8`. The HTML initially displays distances <=0.15 (`isolate_network_default_distance`), adjustable without recomputation. Lower distance means closer sketches. Edge details include the measured distance, matching hashes and Mash p-value. Comparisons sharing at least one hash and with distance <1 are retained in `08_isolate_network.tot/edges.tsv`.

Mash distance is approximate: it is not a VIRIDIC similarity percentage, an aligned-genome fraction or a new ANI/coverage clustering criterion. Short fragments and unequal genome sizes can affect sketch comparisons, and zero distance does not prove exact sequence identity. The identical-sequence view instead uses full-sequence hashes. Existing 95/85 clusters remain unchanged. Connected components do not imply that every pair is within the selected distance, and layout distances have no calibrated evolutionary meaning.

The report opens offline without a server. The pinned Cytoscape.js library is downloaded once, checksum-checked and embedded into the HTML. Controls select sample, assigned host, purity decision, CheckV quality, minimum length, fragment/plasmid evidence and reference source. Labels and colours can be changed, dots dragged, and visible node/edge tables or a PNG exported. Source shapes distinguish own contigs, RefSeq and METAVR; node size follows genome length. Dark borders indicate CheckV Complete, green borders plasmid evidence and dashed borders fragment candidates. Unconnected own contigs remain visible. Reference-only clouds unrelated to the selected own contigs are omitted while browsing own sequences.

Views show individual retained contigs, identical-sequence groups or one retained member per existing MOSAIC cluster. Collapsed views choose the longest retained member among the currently selected samples, with ID breaking ties. Their quality fields and connections belong to that displayed member, not the best score across all members. Member IDs, samples and hosts remain listed for review. No retention flags, cluster memberships, mapping results or PASS/REVIEW/FAIL decisions are changed.

Supporting files are in `FIGURES_AND_TABLES/08_isolate_network.tot/`: `nodes.tsv`, `edges.tsv`, `reference_hits.tsv`, `missing_references.tsv` and `provenance.json`. The provenance records input hashes, selection/sketch settings, Mash version and sequence counts. Prepared sequences, the `.mash.msh` sketch, version file and full `.mash_distances.tsv` table are under `07_ANNOTATION/isolate_network.tot*`. Reference selection/preparation, Mash, library download and report generation each have a benchmark.

Methods: [Mash distances](https://mash.readthedocs.io/en/latest/distances.html); the view is inspired by [PhageClouds](https://doi.org/10.1089/phage.2021.0008). PhageClouds used Dashing sketches and Mash-distance estimates, not Mash's MinHash implementation, so numerical results are not expected to be identical.

## Host prophage activity

`host_prophage_activity=True` (default) adds `NOTEBOOKS/08_host_prophage_activity.tot.ipynb` to isolate runs with host references. Set it to `False` to skip this optional report. It reuses full-read `map_to_host` for each sample's assigned host only; it does not enable the all-host identification test. The existing rule produces both masked and unmasked host mappings in the same job; only the unmasked BAM contributes to activity.

The activity depth rule excludes unmapped, secondary and supplementary alignment records, then writes compressed per-base coverage under `06_MAPPING/HOST/ACTIVITY/`. All QC-passed paired and orphan reads are used. Host contig identifiers are read as text in both the region and depth tables, including numeric names such as `5` or `001`. Both host coverage notebooks pool all usable bases from the sample's assigned host assembly. Predicted viral intervals, whole-contig viral candidates, geNomad plasmid-positive contigs and scaffold-end masks are excluded. Clean portions of an embedded prophage's parent contig still contribute. Uncovered positions are restored as zero; the mean and variance are calculated across bases, not by averaging contig means. Masked positions do not contribute to means or breadth denominators; assessed lengths are included. Saved geNomad plasmid-summary membership supplies the exclusions, without a new score threshold or classifier run. Display binning does not affect statistics.

`per_base_coverage=False` (default) skips optional own-assembly, masked/unmasked
host and sequential-stage per-base exports, but keeps every CoverM summary.
The dedicated `HOST/ACTIVITY/*.basecov.tsv.gz` files remain required when
`host_prophage_activity=True`, independently of the optional export flag.
Set `host_prophage_activity=False` as well to skip those files and both the
activity and enriched-region reports. Enable `per_base_coverage=True` to save
the optional exports as gzip files; existing BAMs and summary tables are unchanged.

Coverage profiles are grouped by host, then prophage region, with visible host and region headings. Hosts are alphabetical; regions follow parent-contig and genomic-position order; samples within each region follow the full sample-ID order, keeping the same region's profiles together across samples. They show the prophage plus 50,000 bp to the left and right by default (`host_prophage_activity_plot_flank_bp`), clipped at contig boundaries. Titles and `coverage_profiles.tsv` include the full host-contig length, predicted prophage-region length and Cohen's d; titles also show fold change. The displayed dashed baseline and activity statistics use the pooled assigned-host background, not just the local plotting flanks. Profile coordinates are recorded as 1-based inclusive `plot_start`/`plot_end`.

The predicted interval is shaded orange for `active`, pink for `ambiguous`, and grey for `dormant` or `not_assessed_*` profiles. The legend and title retain the exact status; grey does not necessarily mean biological inactivity. This display change does not alter activity statistics or thresholds.

The per-base coverage-check table and profile titles show the existing CheckV quality category and estimated completeness (%). Missing assessments are `not reported`; the profile index also records `checkv_quality` and `checkv_completeness`. These sequence-quality assessments do not change coverage thresholds or activity calls. No additional CheckV run is required.

Fold change is the mean prophage depth divided by the pooled assigned-host background mean defined above. Both means include uncovered positions as zero. Zero host-background depth is not assessed; no pseudocount or infinite ratio is used. Changing the plotted flanks does not change these measurements. MOSAIC retains PropagAtE's activity thresholds and unweighted combination of foreground/background variances for Cohen's d, but uses this assembly-wide background rather than PropagAtE's parent-scaffold comparison. It is not an identical stock PropagAtE run.

Each host has a separate enrichment heatmap containing only its own embedded prophages and assigned samples. Hosts and regions without any assessable coverage are omitted from the plots, not the results table; each omitted host has a short explanatory notebook entry. All hosts share a data-derived log2 fold-change colour scale spanning every displayed value, without a fixed upper cap. Assessed cells display the log2 value to one decimal place, with `*` for active and `?` for ambiguous. A dash (`–`) marks no prophage coverage; `NA` marks activity that cannot be assessed. These cells remain uncoloured, rather than showing a misleading log2 value of zero. Hosts and samples are ordered alphabetically. Individual PNG/SVG figures and `heatmaps.tsv` are written under `08_host_prophage_activity.tot/heatmaps/`; the existing overview PNG/SVG contains the same host-specific panels. The host summary TSV, notebook and HTML report also include separate `complete_viral_regions` and `high_quality_viral_regions` counts across embedded and whole-contig predictions; unreported CheckV qualities do not count as either, and Complete is not counted again as High-quality. The host prediction bar plot is ordered alphabetically and retains hosts with zero predictions.

Named activity thresholds use [PropagAtE v1.1 defaults](https://github.com/AnantharamanLab/PropagAtE#flag-descriptions): ratio >=2, Cohen's d >=0.70, prophage mean depth >=1x, breadth >=50% at >=1x, scaffold-end mask 150 bp and minimum prophage length 1,000 bp. The effect uses the unweighted pooled variance in the original implementation. MOSAIC retains its 95% identity/85% aligned-read-length policy, rather than the tool's 97% identity default: this is PropagAtE-style analysis, not an identical stock run.

The TSV `08_host_prophage_activity.tot.tsv` contains all assigned-host sample/region combinations, depth, breadth, background support, ratios, effects and individual threshold checks. Embedded-prophage calls are `active`, `dormant`, `ambiguous` and `not_present`. Zero host background, short regions and fully masked regions receive explicit `not_assessed_*` calls. `active` means coverage enrichment compatible with replication, not confirmed induction; `dormant` does not prove inactivity. The main heatmap and per-region coverage profiles are saved and displayed inside the notebook. This evidence never changes contig retention or purity accounting. Coverage analysis itself uses the existing environment and no new external tool.

### Whole-contig viral candidates

An additional section in the same activity notebook assesses geNomad predictions spanning an entire host-assembly contig. These are not embedded prophages. Each candidate uses the same pooled non-viral, non-plasmid assigned-host background as embedded prophages; the whole candidate contig is excluded from that background. Candidate scaffold ends are masked too. No reference from a different host is used.

The same ratio, effect, depth, breadth and length thresholds apply. Calls are `enriched`, `not_enriched`, `ambiguous`, `not_present` or explicit `not_assessed_*` states. No usable background, or background with zero mean depth, is `not_assessed_host_background`. These enrichment calls do not establish prophage induction or physical packaging, and do not contribute to `active_host_prophage_count`.

The section displays its coverage-check table, per-host annotated heatmaps and coverage profiles, including existing CheckV quality/completeness. Heatmaps retain unassessable candidates as `NA`; `–` means no coverage, `*` enriched and `?` ambiguous. Profiles are orange for enriched, pink for ambiguous and grey otherwise. PNG/SVG files, `heatmaps/whole_heatmaps.tsv` and `coverage_profiles.tsv` are under `FIGURES_AND_TABLES/08_host_prophage_activity.tot/whole_contig_candidates/`. Both evidence scopes remain in the existing activity TSV, with `background_scope` explaining which comparison was used.

The host report includes CheckM, embedded prophages versus whole-contig viral candidates, CheckV and terminal-repeat evidence. The terminal-repeat candidate flag requires explicit geNomad `DTR`/`ITR` topology or the circularity rule's threshold-qualified `terminal_repeat_type` (using `circularity_min_repeat_bp`). `No terminal repeats` and sub-threshold raw overlaps do not qualify. Whole-host-contig repeat evidence is not assigned to an embedded prophage. BACPHLIP is assessed only for CheckV-complete Caudoviricetes genomes; other sequences remain not assessed. Terminal repeats alone do not establish physical circularity, completeness or biological activity.

### Additional exploratory host-region notebook

The same `host_prophage_activity=True` setting also adds `NOTEBOOKS/08_host_enriched_regions.tot.ipynb`. Discovery reuses the assigned-host per-base depth files, host-region metadata and saved geNomad plasmid calls, without new mapping. Its exported loci have a separate annotation step described below. Activity thresholds and 50 kb activity plotting flanks are unchanged.

This coverage-first scan includes host contigs without geNomad predictions. It excludes scaffold-end masks, embedded geNomad regions, whole-contig viral candidates and geNomad plasmid-positive contigs. Discovery uses independent thresholds: full non-overlapping 1 kb bins must have mean depth >=5x, breadth >=80% at >=1x depth and enrichment >=2 against the pooled assigned-host background. Adjacent passing bins form candidate segments of at least 5 kb. Each segment must then pass the ratio, Cohen's d >=0.70, mean-depth and breadth settings using unbinned coverage, additionally excluding that focal segment from the pooled background. Clean portions of its parent contig still contribute. The segment must also have >=2-fold enrichment against local flanks, excluding known embedded prophages and end masks, using 10 kb on each side and clipped at contig boundaries. Zero or unavailable pooled or local background is not assessable and cannot qualify a hit. Proximity to a prophage is report-only, not an exclusion. These are exploratory `unpredicted_enriched_region` entries, not additional active-prophage calls or validated viral boundaries. The predicted-prophage activity thresholds remain unchanged.

Plots show read depth in blue, existing geNomad regions in orange and unpredicted enriched blocks in red. Orange regions include existing CheckV quality/completeness where available; no CheckV assessment is inferred for red blocks. Each block has 200 kb of context on either side, clipped at contig ends. Overlapping display windows are combined across samples for consistent coordinates; separate enriched blocks are not joined across gaps. Figures are ordered by host, parent contig, genomic window and full sample ID. Only samples with a qualifying block in that window are plotted. All assigned samples remain in the summary, including no-hit and unassessable samples; hosts with no hits receive an explanatory notebook entry.

| Parameter | Default | Use |
| --- | --- | --- |
| `host_enriched_regions_min_ratio` | 2.0 | Minimum bin/segment enrichment against pooled assigned-host background. |
| `host_enriched_regions_min_cohen_d` | 0.70 | Minimum segment/background effect size. |
| `host_enriched_regions_min_mean_depth` | 5.0x | Minimum bin and segment mean depth. |
| `host_enriched_regions_breadth_min_depth` | 1.0x | Base depth used for breadth calculations; independent of minimum mean depth. |
| `host_enriched_regions_min_breadth_percent` | 80% | Minimum bin and segment breadth at the configured base depth. |
| `host_enriched_regions_min_local_ratio` | 2.0 | Minimum segment enrichment against available local flanks. |
| `host_enriched_regions_bin_bp` | 1000 bp | Scan resolution; reported boundaries follow the bin grid. |
| `host_enriched_regions_min_length_bp` | 5000 bp | Minimum consecutive enriched segment length. |
| `host_enriched_regions_local_flank_bp` | 10000 bp | Local background on each side for the local-enrichment filter. |
| `host_enriched_regions_plot_flank_bp` | 200000 bp | Display context on each side, not the statistical background. |
| `host_enriched_regions_prophage_adjacency_bp` | 1000 bp | Maximum intervening bases for the report-only `next_to_prophage` flag on the same host contig. |

The segment table is `FIGURES_AND_TABLES/08_host_enriched_regions.tot.tsv`; coordinates are 1-based inclusive. It includes length, segment/background depth, usable background length, `background_scope`, ratio, Cohen's d, breadth, local-flank measurements and distance to the nearest embedded geNomad region. `08_host_enriched_regions_summary.tot.tsv` retains one row per assigned sample, including excluded plasmid-contig count and pooled-background length and mean depth. The dashed profile line uses the pooled baseline for bin scanning; per-segment background measurements additionally exclude that segment. PNG/SVG profiles, `coverage_windows.tsv` and `coverage_profiles.tsv` are saved under `FIGURES_AND_TABLES/08_host_enriched_regions.tot/` and displayed in the notebook. This reporting change reuses existing coverage and does not change isolate contig retention or sequential read-accounting groups.

`extract_host_enriched_regions` exports `HOST/enriched_regions/host_enriched_regions.tot.fasta` and a matching `.tsv` table. Overlapping enriched intervals on the same host/contig are merged across samples and exported once using their union coordinates. Separate, non-overlapping intervals and different hosts/contigs remain separate. The table records coordinates, length, supporting samples, sample-call count and original intervals; the original per-sample segment table is unchanged. Only the enriched sequence is extracted, without the 200 kb plotting flanks. These are coverage-enriched host loci, not confirmed viral sequences. This extraction reuses the completed discovery table and does not repeat mapping, region discovery or plotting.

FASTA descriptions and the merged-locus table include `next_to_prophage`, the nearest embedded-prophage ID, its intervening-base distance and the configured adjacency threshold. The default flag means a gap of at most 1,000 bp on the same host contig; a gap of zero means immediately adjoining or overlapping coordinates. Whole-contig viral candidates and regions on another host/contig are not neighbours for this flag. Missing neighbours are `not reported`, not a zero-distance match.

### Enriched-locus annotation

The isolate target annotates the merged enriched FASTA once and writes `HOST/enriched_regions/host_enriched_regions_metadata.tot.tsv`. Every exported locus remains in this table, regardless of classifier calls or reference hits. These annotations never change the discovery boundaries, contig selection, read accounting or sample decisions.

- Phold runs directly on the nucleotide FASTA in CPU mode, predicts coding sequences and adds protein-function evidence. Native gene-level predictions, annotated GenBank and function counts stay under `07_ANNOTATION/phold_host_enriched_regions.tot/`. The metadata includes CDS counts and reported function/product/confidence summaries. No separate Pharokka run is required. See [Phold FASTA input and CPU usage](https://phold.readthedocs.io/en/latest/run/) and [output definitions](https://phold.readthedocs.io/en/latest/output/).
- geNomad collects virus/plasmid scores, taxonomy, topology and the reported sequence ID. An internal provirus prediction is labelled `viral region` with `provirus region` evidence scope, not a whole-locus viral call. Uncalled whole-contig classification scores remain available when geNomad reports them.
- The existing CheckV rule supplies quality, completeness, contamination and other available fields; the existing TaxMyPhage rule supplies reported taxonomy.
- Existing RefSeq BLAST supplies the best hit and its alignment measurements. With `metavr_blast=True`, the existing METAVR BLAST rule adds its best hit and available reference quality, genome type and ICTV taxonomy. Best-hit ordering matches the all-assembled metadata: lowest e-value, then highest query coverage, identity and alignment length, with accession breaking ties. Missing hits or assessments say `not reported`; disabled METAVR says `not requested`.

Phold uses its own pinned `envs/phold.yaml`; installation is automatic through `downloadPholdDB`. Defaults are `phold_db=db/phold`, `phold_threads=8` and `phold_mem_mb=32768` (32 GiB scheduling reservation). The first run downloads both database and model and needs network access and disk space; subsequent runs reuse them. The metadata shell block reuses `env7.yaml`, including its existing Parquet reader for METAVR, rather than relying on optional packages in the main Snakemake environment. Empty enriched FASTAs produce header-only annotation outputs and a header-only metadata table, without classifier execution.

Coverage enrichment alone cannot establish viral identity, transduction or gene expression. Without DNase treatment before capsid lysis, it cannot distinguish protected/packaged DNA from free host DNA. Repeats, cross-mapping, strain differences and copy-number differences remain alternative explanations. This report does not change retention, purity accounting, activity calls or sample decisions.

## Recovery/purity decisions and separate QC warnings

### Configurable cutoffs

The values below are defaults in `config.yaml`; override them with `--config` when needed. Read alignment identity/length and the independent 95/85 clustering criterion are the shared workflow policies, not additional isolate-specific config switches.

| Parameter | Default | Use |
| --- | --- | --- |
| `min_len` | 1000 bp | Minimum saved SPAdes assembly length, before catalogue creation. |
| `isolate_min_fastp_pairs` | 20000 pairs | QC warning below this number of surviving fastp pairs; independent of recovery/purity decisions. |
| `isolate_min_depth` | 5x | Minimum own-assembly mean depth for retention. |
| `isolate_min_length_bp` | 4000 bp | Minimum length unless the short-contig depth requirement passes. |
| `isolate_short_min_depth` | 10x | Depth needed to retain a shorter saved contig. |
| `isolate_host_min_identity` | 90% | Inclusive assigned-host BLAST identity threshold. |
| `isolate_host_min_query_coverage` | 90% | Inclusive assigned-host BLAST query-coverage threshold; both host tests must pass. |
| `isolate_host_plasmid_min_identity` | 95% | Report-only flag; both sequences must have saved geNomad plasmid-positive calls. |
| `isolate_host_plasmid_min_query_coverage` | 75% | Minimum isolate-contig coverage for the report-only plasmid flag. |
| `isolate_host_plasmid_min_subject_coverage` | 90% | Minimum assigned-host plasmid-contig coverage for the report-only flag. |
| `isolate_min_retained_mapping_percent` | 70% | `REVIEW` below this own-retained read fraction. |
| `isolate_max_host_mapping_percent` | 10% | `REVIEW` above this combined stage-02/03 fraction. |
| `isolate_max_unexplained_percent` | 10% | `REVIEW` above this final unexplained fraction. |
| `sourmash_contig_catalogue_min_length` | 5000 bp | Inclusive representative-screen cutoff; never removes shorter clusters from the catalogue. |
| `sourmash_clean_min_shared_bp` | 50000 bp | Minimum estimated shared bases for the clean-read gather/support screen. |
| `sourmash_clean_min_reference_fraction` | 0.10 | Minimum reference k-mer containment for a supported clean-read match; not mapping breadth. |

`host_prophage_activity` defaults to `True` in isolate runs with hosts. Its ratio/effect/depth/breadth thresholds and plotting flank are described above. `sourmash_clean_reads`, `sourmash_contig_catalogue`, `host_identification_test`, `map_to_RefSeq` and `metavr_blast` default to `False`. The optional host-identification evidence cutoffs are `isolate_host_test_min_reads=100` and `isolate_host_test_min_breadth_percent=1`.

The summary `decision` describes recovery/purity, independently of read-QC warning levels:

- `PASS`: passes the configured recovery/purity checks. This is not proof of a pure single-virus isolate.
- `REVIEW`: retained-contig count !=1, <70% own-retained mapping, >10% host-associated reads, >10% unexplained reads, supported host-prophage exclusions, unassessed/unavailable host references, or an ambiguous/alternative best host in the optional identification test. These concerns do not erase a recovered genome; a supported host-prophage-associated genome excluded by policy remains `REVIEW`.
- `FAIL`: no QC-passed reads, no assembled contigs, or missing, negative/non-finite or inconsistent read accounting. Stage assignments must conserve entering reads, consecutive stages must agree and all stage denominators must match. A failure takes precedence over review warnings.

`warnings` contains recovery/purity review reasons; `failure_reasons` explains failures. `read_accounting_conserved` retains the original total-count check; `read_accounting_valid` additionally checks stage transitions and count availability. Existing biological and QC thresholds are unchanged.

Read-QC levels remain available in `qc_status` and their individual status columns. The `qc_warnings` field describes low-quality filtering, SuperDeduper PCR duplicates, Kraken Eukaryota estimates and too few surviving fastp pairs. These appear in a separate notebook/HTML table with measured values. A QC `WARN` or `FAIL` does not automatically make the isolate decision `REVIEW` or `FAIL`; missing QC measurements remain unassessed (`INFO`/missing), not evidence that QC passed. A sample can therefore have `decision=PASS` and a PCR-duplicate or low-read-count warning.

Existing project reports are not edited in place: rerun the workflow to regenerate them. Previously product-clustered outputs must be regenerated using `--forcerun vOUTclustering` when migrating. The general `07_Normalise.py.ipynb` and virome-positive filtering remain unchanged; independent original 95/85 clustering is the intentional repository-wide change.
