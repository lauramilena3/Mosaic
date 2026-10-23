#!/usr/bin/env python3
"""Draw the virome guide using the isolate manual's compact diagram layout.

Read configuration defaults only; do not execute Snakemake or project analyses.
"""

import argparse
from pathlib import Path

import yaml

from generate_isolate_diagrams import AMBER, BLUE, MUTED, PURPLE, TEAL, Diagram, number


WORKFLOW = Path(__file__).resolve().parents[1]


def diagram(title, subtitle, height=750):
    return Diagram(title, subtitle, height, workflow="viral_metagenome")


def overview(directory, cfg):
    d = diagram("MOSAIC / VIRAL METAGENOME", "Recover a supported viral catalogue and compare communities across samples", 1180)
    qc = d.card(60, 160, 510, "01 / Read QC + cleaning", [
        "fastp → Kraken → biological cleaning", "Keep trimmed reads and QC reports"], BLUE, size=22)
    assembly = d.card(645, 160, 510, "02 / Assembly", [
        "BBnorm → metaSPAdes", "Optional RNA: three extra assemblers"], TEAL, size=22)
    discovery = d.card(1230, 160, 510, "03 / Viral discovery", [
        "DNA: geNomad-positive sequences", "Beta RNA stays outside main catalogue"], TEAL, size=21)
    d.across(qc, assembly, BLUE)
    d.across(assembly, discovery, TEAL)
    cluster = d.card(1230, 415, 510, "04 / Shared catalogue", [
        "MMseqs exact dereplication → 95/85 ANI", "CheckV-guided representatives"], TEAL, size=21)
    final = d.card(645, 415, 510, "Viral + mapping support", [
        "Standard or default geNomad-only calls", "One main filtered FASTA"], TEAL, size=21)
    mapping = d.card(60, 415, 510, "05 / Map + normalise", [
        "Full cleaned pairs → final catalogue", "Counts / RPKM / breadth / depth"], BLUE, size=22)
    d.down(discovery, cluster, TEAL)
    d.arrow([(1230, 486), (1155, 486)], TEAL)
    d.arrow([(645, 486), (570, 486)], BLUE)
    reports = d.card(60, 690, 510, "07 / Community reports", [
        "Quality / abundance / hosts / prevalence", "Executed notebooks + tables + figures"], MUTED, size=21)
    annotation = d.card(645, 690, 510, "06 / Annotation", [
        "Genes / taxonomy / host / lifestyle", "RefSeq similarity; optional METAVR"], PURPLE, size=22)
    branches = d.card(1230, 690, 510, "Independent optional branches", [
        "RefSeq detection / microdiversity", "All-assembled / beta RNA exploration"], AMBER, size=21)
    d.down(mapping, reports, BLUE)
    d.down(final, annotation, PURPLE)
    d.arrow([(645, 761), (570, 761)], PURPLE)
    d.band(900, "READ SUBSETS ARE NOT INTERCHANGEABLE", "Full clean pairs: final abundance / RefSeq / variation · 2M-named subset: diagnostics / all-assembled", BLUE)
    d.note(1030, "Shared downloads: Kraken · geNomad · CheckV · VirSorter2 · VIBRANT · RefSeq · phage annotation/host databases")
    d.note(1080, "Blue: reads/mapping · teal: contigs/catalogue · purple: annotations/references · amber: optional branches · grey: reports")
    d.save(directory, "viral_metagenome_overview")


def read_qc(directory, cfg):
    d = diagram("01 / READ QC AND BIOLOGICAL CLEANING", "Virome mode removes configured biological contaminants; isolate mode does not", 880)
    raw = d.card(60, 160, 490, "Raw paired FASTQs", ["FastQC / counts / SuperDeduper", "Duplicate estimate: report only"], BLUE, size=22)
    trim = d.card(650, 160, 1090, "fastp → trimmed paired reads + merged orphans", [
        "PE overlap analysis + adapter detection; no adapter FASTA",
        f"Poly-G ≥{number(cfg['fastp_poly_g_min_len'])} nt; front/right quality trimming",
        f"Keep reads ≥{number(cfg['fastp_length_required'])} bp; save HTML/JSON reports"], BLUE)
    d.across(raw, trim, BLUE)
    kraken = d.card(650, 440, 490, "Kraken + cleaning", [
        "Remove Eukaryota by default", "BBDuk if contaminant references exist", "Empty list: bypass BBDuk"], BLUE, size=22)
    clean = d.card(1250, 440, 490, "Cleaned libraries", [
        "Paired FASTQs + cleaned orphans", "Post-QC / post-cleaning reporting", "Retain trimmed files for reuse"], BLUE, size=22)
    d.down(trim, kraken, BLUE)
    d.across(kraken, clean, BLUE)
    d.card(60, 440, 490, "QC interpretation", [
        "Filtering loss / PCR duplicates", "Kraken eukaryotic fraction", "K-mer sampling evidence"], MUTED, size=22)
    d.band(710, "DEFAULTS", "remove_euk=True · contaminants_list empty · no implicit PhiX-reference removal", BLUE)
    d.save(directory, "viral_metagenome_01_read_qc")


def assembly(directory, cfg):
    d = diagram("02 / DNA AND OPTIONAL RNA ASSEMBLY", "Saved FASTAs have persistent IDs; provenance stays beside each assembly", 1010)
    dna = d.card(60, 160, 500, "Standard DNA branch", [
        f"BBnorm target {number(cfg['max_norm'])}× → metaSPAdes",
        f"Save scaffolds ≥{int(cfg['min_len']):,} bp", "QUAST + DNA assembly notebook"], TEAL, size=22)
    spades = d.card(650, 160, 500, "RNAviralSPAdes", ["BBNorm-normalized paired reads", "SPAdes --rnaviral, not --rna"], TEAL, size=21)
    megahit = d.card(1240, 160, 500, "MEGAHIT", ["BBNorm-normalized paired reads", "Independent assembler output"], TEAL, size=21)
    trinity = d.card(650, 420, 500, "Trinity", [
        "BBNorm-normalized paired reads", "--no_salmon; --full_cleanup", f"Default workers: {int(cfg['rna_trinity_threads'])}"], TEAL, size=21)
    saved = d.card(1240, 420, 500, "Each saved RNA FASTA", [
        f"Filter length ≥{int(cfg['rna_min_contig_length']):,} nt", "Then rename / geNomad / repeats", "Flat files in 03_CONTIGS/RNA/"], TEAL, size=22)
    d.arrow([(900, spades[1] + spades[3]), (900, 375), (1390, 375), (1390, 420)], TEAL)
    d.down(megahit, saved, TEAL)
    d.across(trinity, saved, TEAL)
    d.text(650, 126, "RNA_enriched=True: beta exploration for all samples", size=23, color=AMBER)
    combined = d.card(650, 710, 1090, "Per-sample combined_derreplicated FASTA", [
        "Pool selected RNA assemblers; collapse exact duplicates / reverse complements",
        "Not yet viral-positive and not the final vOTU catalogue"], TEAL, size=22)
    d.down(saved, combined, TEAL)
    d.card(60, 710, 500, "Persistent names + provenance", [
        "sample_assembler_00001_len_length", ".ids.tsv = provenance only"], TEAL, size=21)
    d.note(930, "Beta RNA feeds optional all-assembled analysis, not the main catalogue; DNA assembly remains enabled.")
    d.save(directory, "viral_metagenome_02_assembly")


def discovery(directory, cfg):
    d = diagram("03 / VIRAL DISCOVERY AND TERMINAL REPEATS", "Main DNA candidate extraction is separate from beta RNA exploration", 860)
    dna = d.card(60, 160, 500, "Saved DNA scaffolds", ["geNomad relaxed end-to-end", "Extract virus / provirus candidates"], TEAL, size=22)
    rna = d.card(650, 160, 500, "Beta RNA exploration", [
        f"VirSorter2 groups: {cfg['rna_viral_groups']}", "Keep original candidate sequences"], TEAL, size=22)
    candidates = d.card(1240, 160, 500, "Exploratory RNA CheckV", [
        "Candidate quality / completeness", "Kept outside main quality summary"], PURPLE, size=21)
    d.across(rna, candidates, PURPLE)
    repeats = d.card(60, 435, 1090, "Every saved DNA/RNA assembly: terminal-repeat table", [
        f"Exact DTR / ITR threshold: ≥{int(cfg['circularity_min_repeat_bp'])} bp",
        "One row per contig; sequence, length and end coordinates",
        "Non-overlapping ends; geNomad topology shown alongside"], TEAL, size=22)
    individual = d.card(1240, 435, 500, "Each RNA assembler: geNomad", [
        "RNAviralSPAdes / MEGAHIT / Trinity", "Evidence before combining", "Exploratory annotation only"], PURPLE, size=21)
    d.band(690, "INTERPRETATION", "Terminal repeats are sequence-end evidence, not proof of a circular molecule or a complete genome.", MUTED)
    d.save(directory, "viral_metagenome_03_discovery")


def catalogue(directory, cfg):
    d = diagram("04 / MAIN vOTU CATALOGUE AND FINAL ACCEPTANCE", "Independent ANI and coverage thresholds; beta RNA assemblies stay separate", 1070)
    pool = d.card(60, 160, 500, "Pool viral candidates", ["DNA + optional subassemblies", "Optional user reference sequences"], TEAL, size=22)
    exact = d.card(650, 160, 500, "MMseqs dereplication", ["100% identity", "100% target coverage"], TEAL, size=22)
    cluster = d.card(1240, 160, 500, "CheckV ANI clustering", ["ANI ≥95% AND target coverage ≥85%", "Original aniclust; min_qcov=0"], TEAL, size=21)
    d.across(pool, exact, TEAL)
    d.across(exact, cluster, TEAL)
    rep = d.card(650, 420, 1090, "Select representatives and collect viral evidence", [
        "CheckV-guided selection, not simply the longest sequence",
        "geNomad (including conservative) + VIBRANT + VirSorter2"], PURPLE, size=22)
    d.down(cluster, rep, PURPLE)
    standard = d.card(60, 700, 500, "Standard acceptance", [
        "Three-tool agreement OR", "conservative geNomad support", f"Length ≥{int(cfg['min_votu_length']):,} nt OR", "CheckV MIUViG High-quality"], TEAL, size=22)
    genomad = d.card(650, 700, 500, "genomad_only=True", [
        "Default geNomad viral call", "Other tools: metadata only", f"Length ≥{int(cfg['min_votu_length']):,} nt OR",
        "CheckV MIUViG High-quality"], TEAL, size=22)
    final = d.card(1240, 700, 500, "Both modes need mapping", [
        ">1 passing mapped read", "In ≥1 configured key sample", "Save filtered_95-85_...tot.fasta", "Record the acceptance route"], BLUE, size=21)
    d.arrow([(1195, 565), (1195, 640), (310, 640), (310, 700)], TEAL)
    d.arrow([(1195, 640), (900, 640), (900, 700)], TEAL)
    d.across(genomad, final, BLUE)
    d.arrow([(310, standard[1] + standard[3]), (310, 960),
             (1490, 960), (1490, final[1] + final[3])], BLUE)
    d.note(995, "No RNA acceptance shortcut. The modified product-based helper is not used; clustering stays unchanged.")
    d.save(directory, "viral_metagenome_04_catalogue")


def mapping(directory, cfg):
    d = diagram("05 / MAPPING AND ABUNDANCE NORMALISATION", "Do not confuse BBNorm assembly inputs, diagnostic subsets and final abundance reads", 1020)
    full = d.card(60, 160, 500, "Full cleaned paired reads", [
        "Final filtered vOTUs: tot abundance", "RefSeq / variation: separate branches"], BLUE, size=22)
    diagnostic = d.card(650, 160, 500, "2M-named paired subset", [
        "BBTools reads=2000000; seed=1", "Assembly / unfiltered diagnostics", "Optional all-assembled catalogue"], BLUE, size=22)
    filter_reads = d.card(1240, 160, 500, "Bowtie2 + CoverM", [
        "≥95% read identity", "≥85% aligned read length", "All / operational unique summaries"], BLUE, size=22)
    stats = d.card(60, 455, 1090, "Coverage tables → 07_Normalise notebook", [
        "Raw counts / RPKM + ambiguity-adjusted counts / RPKM",
        "Breadth (%) / covered bases; mean depth in per-sample CoverM files",
        "Normalised RPKM denominator: adjusted mapped reads, not library size"], BLUE, size=22)
    coverage = d.card(1240, 455, 500, "Filtered abundance cells", [
        "Reference ≥6,667 bp:", "  >5,000 covered bases", "Shorter reference:", "  >75% breadth"], TEAL, size=22)
    d.down(full, stats, BLUE)
    d.down(filter_reads, coverage, TEAL)
    d.card(60, 755, 1090, "Optional per-base exports", [
        "per_base_coverage=True → *_basecov.txt.gz",
        "Nonzero positions; 0-based coordinates; summaries always saved"], BLUE, size=22)
    d.card(1240, 755, 500, "What the filter changes", [
        "Nonpassing matrix cells become 0", "FASTA records are not removed"], TEAL, size=22)
    d.note(948, "raw RPKM = before abundance adjustment, not mapping raw FASTQs. Counts are not molecule counts.")
    d.save(directory, "viral_metagenome_05_mapping")


def annotation(directory, cfg):
    d = diagram("06 / ANNOTATION AND OPTIONAL COMPARISONS", "Annotations add evidence; independent reference mapping does not replace the catalogue", 990)
    final = d.card(60, 160, 500, "Final filtered catalogue", ["Main-catalogue representatives", "The input to core annotation"], TEAL, size=22)
    core = d.card(650, 160, 1090, "Core evidence", [
        "Pharokka / Phynteny / gene and nucleotide composition",
        "PhaGCN / TaxMyPhage / RefSeq BLAST / DIAMOND / AAI",
        "iPHoP host predictions / BACPHLIP lifestyle evidence"], PURPLE, size=22)
    d.across(final, core, PURPLE)
    reference = d.card(60, 470, 500, "Reference/community context", [
        "METAVR BLAST + metadata", "vConTACT2 protein-sharing context", "DRAM-v functional annotation"], PURPLE, size=22)
    refseq = d.card(650, 470, 500, "Independent RefSeq detection", [
        "Map full cleaned pairs to RefSeq", "Confidence categories + NC comparison", "Automatic cached NCBI metadata"], AMBER, size=21)
    all_assembled = d.card(1240, 470, 500, "All-assembled exploration", [
        f"Top {int(cfg['all_assembled_top_per_sample'])}/non-NC sample +", f"top {int(cfg['all_assembled_cluster_top_n']):,} by mean raw RPKM",
        "Cluster selection; keep all reps", "No remapping or member sums"], AMBER, size=22)
    d.band(780, "OPTIONAL INPUTS", "Existing microbial spacer FASTA → SpacePHARER · selected all-assembled representatives → optional Cenote", PURPLE)
    d.note(914, "Phage-oriented tools may be uninformative for RNA/eukaryotic viruses; absence of an annotation is not absence of a virus.")
    d.save(directory, "viral_metagenome_06_annotation")


def reports(directory, cfg):
    d = diagram("07 / REPORTS AND INTERPRETATION", "Start with final metadata and abundance; return to the intermediate evidence when needed", 900)
    main = d.card(60, 160, 1090, "10_abundance_lifestyle_summary", [
        "vOTU_metadata_polished.csv: all final catalogue sequences",
        "Sample / quality / taxonomy / host / lifestyle summaries",
        "Presence, prevalence, relative abundance and sample accumulation"], MUTED, size=22)
    stages = d.card(1240, 160, 500, "Stage reports", [
        "01 QC / 03 assembly / 04 viral ID", "05 representative selection/filtering", "07 normalisation / mapping stats"], MUTED, size=21)
    optional = d.card(60, 475, 500, "Optional reports", [
        "03 RNA assembly / 08 METAVR", "07 RefSeq detection", "11 microdiversity (separate)"], AMBER, size=22)
    data = d.card(650, 475, 500, "Interpret abundance correctly", [
        "Presence: filtered normalised RPKM >0", "Relative abundance: catalogue share", "NC not automatically removed"], BLUE, size=21)
    evidence = d.card(1240, 475, 500, "Keep uncertainty visible", [
        "Similarity ≠ confirmed species", "Terminal repeats ≠ circular proof", "Read-level diversity ≠ activity"], PURPLE, size=22)
    d.band(740, "FILES", "NOTEBOOKS: executed analysis · FIGURES_AND_TABLES: reports/tables/plots · BENCHMARK: per-job measurements", MUTED)
    d.save(directory, "viral_metagenome_07_reports")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=WORKFLOW / "config.yaml")
    parser.add_argument("--output-dir", type=Path, default=WORKFLOW / "docs/figures")
    args = parser.parse_args()
    with args.config.open() as handle:
        cfg = yaml.safe_load(handle)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for draw in [overview, read_qc, assembly, discovery, catalogue, mapping, annotation, reports]:
        draw(args.output_dir, cfg)
    print(f"Saved 8 PNG/SVG diagram pairs in {args.output_dir}")


if __name__ == "__main__":
    main()
