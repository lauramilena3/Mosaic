#!/usr/bin/env python3
"""Draw the isolate guide's conceptual overview and seven step diagrams.

These are manually arranged biological methods diagrams.
Configured thresholds come from config.yaml; no workflow or project is executed.
"""

import argparse
from io import StringIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import yaml


WORKFLOW = Path(__file__).resolve().parents[1]
DPI = 150
INK = "#172b4d"
MUTED = "#526277"
BLUE = "#3A72ED"
TEAL = "#0f766e"
PURPLE = "#7c3aed"
AMBER = "#b45309"
RED = "#b42318"
FILLS = {BLUE: "#f1f6ff", TEAL: "#eff9f6", PURPLE: "#f7f3ff",
         AMBER: "#fff7ed", RED: "#fff1f0", MUTED: "#f3f5f8"}
READ_ASSIGNMENT_COLORS = {
    "01_own_retained": BLUE,
    "02_host_chromosomes": PURPLE,
    "03_host_viral": "#b69aef",
    "04_own_excluded": AMBER,
    "05_other_retained": "#93b4f7",
    "06_other_excluded": "#e6b07c",
    "07_unexplained": "#8a8a8a",
}
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                     "svg.hashsalt": "mosaic-isolate-methods"})


def number(value):
    return f"{float(value):g}"


class Diagram:
    """Pixel-based layout, with text-fitted cards and checked text boundaries."""

    def __init__(self, title, subtitle, height, width=1800, workflow="phage_isolates"):
        self.width, self.height = width, height
        self.workflow = workflow
        self.fig = plt.figure(figsize=(width / DPI, height / DPI), dpi=DPI, facecolor="white")
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set(xlim=(0, width), ylim=(height, 0))
        self.ax.axis("off")
        self.labels = []
        self.text(60, 35, title, size=34, weight="bold")
        self.text(60, 88, subtitle, size=23, color=MUTED)

    def text(self, x, y, value, size=23, color=INK, weight="normal", ha="left", bounds=None):
        label = self.ax.text(x, y, value, fontsize=size * 72 / DPI, color=color,
                             weight=weight, ha=ha, va="top", linespacing=1.35)
        self.labels.append((label, bounds))
        return label

    def rect(self, x, y, width, height, color, fill=None, radius=14):
        self.ax.add_patch(FancyBboxPatch((x, y), width, height,
                          boxstyle=f"round,pad=0,rounding_size={radius}", linewidth=1.1,
                          edgecolor=color, facecolor=fill or FILLS[color]))

    def card(self, x, y, width, title, lines, color, size=23):
        # Heights follow the actual line count: no fixed oversized card boxes.
        title_rows = title.count("\n") + 1
        body_rows = sum(line.count("\n") + 1 for line in lines)
        title_height = title_rows * 33
        body_height = body_rows * 31 if lines else 0
        height = 36 + title_height + (10 + body_height if lines else 0)
        self.rect(x, y, width, height, color, fill="white")
        self.ax.plot([x + 15, x + 15], [y + 19, y + height - 19],
                     color=color, linewidth=3, solid_capstyle="round")
        bounds = (x + 25, y + 8, x + width - 10, y + height - 8)
        heading = self.text(x + 32, y + 18, title, size=27, color=color, weight="bold", bounds=bounds)
        heading_width = heading.get_window_extent(self.fig.canvas.get_renderer()).width
        if heading_width > width - 48:
            heading.set_fontsize(max(21, 27 * (width - 48) / heading_width) * 72 / DPI)
        body_y = y + 18 + title_height + 10
        for line in lines:
            self.text(x + 32, body_y, line, size=size, bounds=bounds)
            body_y += 31 * (line.count("\n") + 1)
        return (x, y, width, height)

    def arrow(self, points, color=MUTED, dashed=False):
        style = "--" if dashed else "-"
        if len(points) > 2:
            self.ax.plot([p[0] for p in points[:-1]], [p[1] for p in points[:-1]],
                         color=color, linewidth=1.4, linestyle=style)
        self.ax.add_patch(FancyArrowPatch(points[-2], points[-1], arrowstyle="-|>",
                          mutation_scale=14, linewidth=1.4, color=color,
                          linestyle=style, shrinkA=0, shrinkB=3))

    def down(self, upper, lower, color=MUTED):
        self.arrow([(upper[0] + upper[2] / 2, upper[1] + upper[3]),
                    (lower[0] + lower[2] / 2, lower[1])], color)

    def across(self, left, right, color=MUTED):
        self.arrow([(left[0] + left[2], left[1] + left[3] / 2),
                    (right[0], right[1] + right[3] / 2)], color)

    def note(self, y, value, color=MUTED):
        self.text(60, y, value, size=22, color=color)

    def band(self, y, title, value, color, height=100):
        self.rect(60, y, self.width - 120, height, color)
        self.text(82, y + 15, title, size=24, color=color, weight="bold")
        self.text(82, y + 53, value, size=23)

    def save(self, directory, stem):
        self.note(self.height - 34, f"MOSAIC / {self.workflow} · conceptual methods diagram · config.yaml defaults")
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        for label, bounds in self.labels:
            box = label.get_window_extent(renderer).transformed(self.ax.transData.inverted())
            left, right = sorted([box.x0, box.x1])
            top, bottom = sorted([box.y0, box.y1])
            limits = bounds or (8, 8, self.width - 8, self.height - 8)
            if left < limits[0] or top < limits[1] or right > limits[2] or bottom > limits[3]:
                raise ValueError(f"Text outside layout in {stem}: {label.get_text()!r}")
        self.fig.savefig(directory / (stem + ".png"), dpi=DPI)
        svg = StringIO()
        self.fig.savefig(svg, format="svg", metadata={"Date": None})
        (directory / (stem + ".svg")).write_text(
            "\n".join(line.rstrip() for line in svg.getvalue().splitlines()) + "\n", encoding="utf-8")
        plt.close(self.fig)


def overview(directory, cfg):
    d = Diagram("MOSAIC / PHAGE ISOLATES", "Genome recovery, purity and relatedness", 1190)
    lanes = [(60, BLUE, "READS"), (650, TEAL, "ASSEMBLED CONTIGS"),
             (1240, PURPLE, "ASSIGNED HOSTS")]
    for x, color, title in lanes:
        d.rect(x, 155, 500, 620, color, fill=FILLS[color])
        d.text(x + 20, 173, title, size=24, color=color, weight="bold")
    qc = d.card(80, 220, 460, "01 / Read QC", ["fastp + contamination reporting", "Keep all QC-passed reads"], BLUE, size=22)
    reads = d.card(80, 432, 460, "Full-read mapping", ["Paired reads + orphans; never 2M", "Depth for per-sample selection"], BLUE, size=22)
    account = d.card(80, 638, 460, "06 / Read accounting", ["Six priority stages + unexplained"], BLUE, size=22)
    d.down(qc, reads, BLUE)
    d.down(reads, account, BLUE)
    assembly = d.card(670, 220, 460, "02 / Assembly + evidence", ["BBnorm → SPAdes → saved contigs", "geNomad / CheckV / terminal repeats"], TEAL, size=21)
    selection = d.card(670, 432, 460, "04 / Retained or excluded", ["Own-assembly depth + length", "Assigned-host BLAST evidence"], TEAL, size=22)
    sets = d.card(670, 638, 460, "Both sets preserved", ["No viral-positive selection gate"], TEAL, size=22)
    d.down(assembly, selection, TEAL)
    d.down(selection, sets, TEAL)
    d.arrow([(670, 697), (540, 697)], TEAL)
    d.across(qc, assembly, BLUE)
    d.across(reads, selection, BLUE)
    hosts = d.card(1260, 220, 460, "Host FASTAs + assignments", ["Explicit sample–host metadata"], PURPLE, size=22)
    evidence = d.card(1260, 432, 460, "03 / Host characterisation", ["CheckM / geNomad / CheckV", "Chromosomes + viral regions"], PURPLE, size=22)
    activity = d.card(1260, 638, 460, "Host activity evidence", ["Coverage enrichment; report only"], PURPLE, size=22)
    d.down(hosts, evidence, PURPLE)
    d.down(evidence, activity, PURPLE)
    d.arrow([(1260, 508), (1130, 508)], PURPLE)
    d.text(1195, 466, "BLAST", size=19, color=PURPLE, ha="center")
    cluster = d.card(650, 830, 500, "05 / Dereplication + clustering", ["All saved contigs: retained + excluded", "MMseqs exact reps → MOSAIC clusters"], TEAL, size=22)
    d.down(sets, cluster, TEAL)
    report = d.card(1240, 830, 500, "07 / Reports + review", ["Sample / contig / cluster / host tables", "PASS / REVIEW / FAIL + figures"], MUTED, size=22)
    d.down(activity, report, PURPLE)
    d.arrow([(1150, 900), (1240, 900)], TEAL)
    d.arrow([(310, 748), (310, 1010), (1490, 1010), (1490, 971)], BLUE)
    d.text(60, 1060, "COLOURS", size=20, color=MUTED, weight="bold")
    for x, label, color in [(220, "Reads", BLUE), (385, "Contigs", TEAL),
                            (570, "Hosts", PURPLE), (735, "Excluded / review", AMBER),
                            (1050, "Unexplained / reports", MUTED)]:
        d.text(x, 1060, label, size=22, color=color, weight="bold")
    d.note(1108, "Shared downloads: Kraken · geNomad · CheckV · CheckM · RefSeq · optional Sourmash / METAVR")
    d.save(directory, "phage_isolates_overview")


def read_qc(directory, cfg):
    d = Diagram("01 / READ QC", "Trim technical sequence; report contamination without removing biological reads", 820)
    raw = d.card(60, 165, 470, "Paired raw reads", ["R1 + R2 FASTQs", "FastQC and raw read counts"], BLUE)
    fastp = d.card(660, 165, 1080, "fastp", [
        "Paired-end overlap analysis + adapter detection; no adapter FASTA",
        f"Poly-G trimming: minimum run {number(cfg['fastp_poly_g_min_len'])} nt",
        f"Front: Q{number(cfg['fastp_cut_front_mean_quality'])} / {number(cfg['fastp_cut_front_window_size'])}-base window; right: Q{number(cfg['fastp_cut_right_mean_quality'])} / {number(cfg['fastp_cut_right_window_size'])}-base window",
        f"Keep reads ≥{number(cfg['fastp_length_required'])} bp; qualified bases Q{number(cfg['fastp_qualified_quality_phred'])}; at most {number(cfg['fastp_unqualified_percent_limit'])}% unqualified / {number(cfg['fastp_n_base_limit'])} Ns"], BLUE)
    d.across(raw, fastp, BLUE)
    clean = d.card(660, 462, 1080, "All QC-passed reads → mapping", [
        "Keep paired reads and orphans; no 2M subset in isolate mode",
        "Bypass configured-contaminant and eukaryotic read removal"], BLUE)
    d.down(fastp, clean, BLUE)
    duplicates = d.card(60, 420, 470, "SuperDeduper", ["Raw-read PCR-duplicate estimate", "QC report only; not assembly input"], BLUE, size=22)
    d.down(raw, duplicates, BLUE)
    d.band(654, "REPORTING", "Kraken once on QC-passed reads · FastQC / MultiQC · read counts · separate QC warnings", BLUE)
    d.save(directory, "phage_isolates_01_read_qc")


def assembly(directory, cfg):
    d = Diagram("02 / ASSEMBLY AND CONTIG EVIDENCE", "Assembly normalisation does not change the full-read mapping denominator", 730)
    qc = d.card(60, 165, 490, "QC-passed assembly reads", ["BBnorm normalisation", f"Configured target: {number(cfg['max_norm'])}×"], BLUE)
    spades = d.card(655, 165, 490, "SPAdes", ["DNA isolate assembly mode", "No mandatory viral classifier gate"], TEAL)
    saved = d.card(1250, 165, 490, "Saved contigs", [f"Length ≥{int(cfg['min_len']):,} bp", "Persistent IDs + .ids.tsv provenance"], TEAL, size=22)
    d.across(qc, spades, TEAL)
    d.across(spades, saved, TEAL)
    annotations = d.card(655, 410, 1085, "Run evidence tools on the saved assembly", [
        "geNomad: virus / plasmid classification, scores, taxonomy and topology",
        "CheckV: quality, completeness, contamination and provirus evidence",
        "Terminal repeats: type, sequence, length and end coordinates"], TEAL)
    d.ax.plot([1495, 1495, 900], [306, 355, 355], color=TEAL, linewidth=1.4)
    d.arrow([(900, 355), (1197, 355), (1197, 410)], TEAL)
    d.arrow([(900, 355), (305, 355), (305, 410)], BLUE)
    d.card(60, 410, 490, "Own-assembly mapping", ["All QC-passed paired + orphan reads", "≥95% identity AND ≥85% aligned length", "Per-contig depth / breadth / counts"], BLUE, size=21)
    d.note(648, "Evidence is collected for every saved contig. Unclassified contigs remain eligible for retention.")
    d.save(directory, "phage_isolates_02_assembly_evidence")


def hosts(directory, cfg):
    d = Diagram("03 / HOST CHARACTERISATION", "The assigned host drives selection and residual-read accounting; other hosts are report-only", 1100)
    references = d.card(60, 160, 500, "Host FASTAs + assignments", [
        "HOST/{host}.fasta", "host_mapping_file.tsv: sample → host", "Host FASTAs require explicit metadata"], PURPLE, size=22)
    quality = d.card(650, 160, 500, "CheckM + geNomad", [
        "Host completeness / contamination", "Embedded prophages vs whole-contig", "viral candidates + scope / coordinates"], PURPLE, size=21)
    viral = d.card(1240, 160, 500, "Predicted viral regions", [
        "CheckV + terminal-repeat evidence", "BACPHLIP: CheckV Complete", "Caudoviricetes genomes only"], PURPLE, size=22)
    d.across(references, quality, PURPLE)
    d.across(quality, viral, PURPLE)
    stage_refs = d.card(60, 425, 500, "Assigned-host references", [
        "Chromosomes: viral regions masked", "Viral reference: embedded prophages +", "whole-contig candidates; scope labelled"], PURPLE, size=21)
    d.down(references, stage_refs, PURPLE)
    activity = d.card(650, 425, 1090, "Optional host prophage activity report (enabled by default)", [
        "Use assigned-host unmasked full-read BAM; assess embedded prophages only",
        f"Depth ratio ≥{number(cfg['host_prophage_activity_min_ratio'])}; Cohen’s d ≥{float(cfg['host_prophage_activity_min_cohen_d']):.2f}",
        f"Prophage mean depth ≥{number(cfg['host_prophage_activity_min_mean_depth'])}×; breadth ≥{number(cfg['host_prophage_activity_min_breadth_percent'])}% at depth ≥1×",
        f"Region length ≥{int(cfg['host_prophage_activity_min_length_bp']):,} bp; mask {int(cfg['host_prophage_activity_mask_bp'])} bp at scaffold ends"], PURPLE, size=23)
    d.arrow([(1490, 338), (1490, 395), (1195, 395), (1195, 425)], PURPLE)
    d.band(720, "ACTIVITY BACKGROUND", "Whole parent scaffold, excluding all embedded prophages and masked ends; valid background required.", PURPLE)
    d.note(844, f"Profiles show ±{int(cfg['host_prophage_activity_plot_flank_bp']):,} bp, clipped to the scaffold. The display window does not set the statistics.")
    d.note(886, "Coverage enrichment is compatible with replication, not confirmed induction; it does not change retention.")
    d.band(947, "OPTIONAL ALL-HOST IDENTIFICATION TEST", "Separate masked / unmasked comparison against every host; report evidence, never automatically relabel samples.", PURPLE)
    d.save(directory, "phage_isolates_03_hosts")


def selection(directory, cfg):
    d = Diagram("04 / RETAINED AND EXCLUDED CONTIGS", "A length/depth and assigned-host decision — not a geNomad-positive filter", 790)
    d.card(60, 160, 780, "Own assembly + full-read mapping", [
        f"Saved contigs ≥{int(cfg['min_len']):,} bp; use own-assembly mean depth",
        f"Depth ≥{number(cfg['isolate_min_depth'])}× AND (length ≥{int(cfg['isolate_min_length_bp']):,} bp OR depth ≥{number(cfg['isolate_short_min_depth'])}×)"], TEAL)
    d.card(960, 160, 780, "Assigned-host BLAST screen", [
        f"Identity ≥{number(cfg['isolate_host_min_identity'])}% AND query coverage ≥{number(cfg['isolate_host_min_query_coverage'])}%",
        "Query coverage is the fraction of the isolate contig",
        "Match assigned chromosome OR predicted host viral regions"], PURPLE, size=22)
    d.arrow([(450, 307), (450, 355), (900, 355), (900, 395)], TEAL)
    d.arrow([(1350, 338), (1350, 355), (900, 355)], PURPLE)
    d.band(395, "RETAIN ONLY WHEN", "Coverage rule passes AND there is no qualifying assigned-host chromosome / viral-region match.", TEAL)
    d.card(60, 555, 780, "Retained FASTA", ["Coverage-supported; no qualifying assigned-host match", "geNomad-negative / unclassified contigs can be retained"], TEAL, size=23)
    d.card(960, 555, 780, "Excluded FASTA + reason", ["Length/depth fails OR qualifying assigned-host match", "Supported host-prophage matches are excluded + flagged"], AMBER, size=23)
    d.arrow([(500, 495), (450, 555)], TEAL)
    d.arrow([(1300, 495), (1350, 555)], AMBER)
    d.note(717, "Both sets are saved and enter the global catalogue. Other-host matches do not drive exclusion.")
    d.save(directory, "phage_isolates_04_selection")


def clustering(directory, cfg):
    d = Diagram("05 / DEREPLICATION AND CLUSTERING", "All saved contigs, including excluded sequences, contribute to relatedness", 765)
    original = d.card(60, 160, 490, "Global contig set", ["Retained + excluded, across samples", "No top-N or viral-positive restriction"], TEAL, size=22)
    exact = d.card(655, 160, 490, "MMseqs exact representatives", ["100% identity", "100% target coverage"], TEAL)
    cluster = d.card(1250, 160, 490, "MOSAIC clusters", ["ANI ≥95% AND target coverage ≥85%", "Independent tests; query cutoff = 0"], TEAL, size=21)
    d.across(original, exact, TEAL)
    d.across(exact, cluster, TEAL)
    d.band(375, "PROVENANCE", "original_contig → exact_rep → mosaic_cluster → cluster_rep", TEAL)
    d.card(60, 535, 780, "Representative sequence + similarity", ["Extract one representative per cluster; RefSeq BLAST", f"Optional: Sourmash reps ≥{int(cfg['sourmash_contig_catalogue_min_length']):,} bp / METAVR BLAST"], TEAL)
    d.card(960, 535, 780, "Cluster metadata", ["Original / exact-rep / sample counts; retained + excluded", "Member annotations kept distinct from representative quality"], TEAL, size=22)
    d.arrow([(1495, 307), (1495, 350), (900, 350), (900, 375)], TEAL)
    d.arrow([(450, 475), (450, 535)], TEAL)
    d.arrow([(1350, 475), (1350, 535)], TEAL)
    d.note(696, "Use the original installed aniclust: --min_ani 95 --min_tcov 85 --min_qcov 0; never the modified product helper.")
    d.save(directory, "phage_isolates_05_clustering")


def accounting(directory, cfg):
    d = Diagram("06 / WHERE DO THE READS GO?", "All QC-passed paired and orphan reads · ≥95% identity AND ≥85% of read length aligned", 660, width=2240)
    d.band(145, "PRIORITY ASSIGNMENT", "Only reads without a passing alignment advance. Each read is counted once; the first stage reuses the own-assembly BAM.", BLUE)
    groups = [("01", "Own", "retained", BLUE), ("02", "Masked host", "chromosome", PURPLE),
              ("03", "Host viral", "regions", PURPLE), ("04", "Own", "excluded", AMBER),
              ("05", "Other samples’", "retained", BLUE), ("06", "Other samples’", "excluded", AMBER),
              ("REST", "Unexplained", "reads", MUTED)]
    for index, ((stage, first, second, family), color) in enumerate(zip(groups, READ_ASSIGNMENT_COLORS.values())):
        x = 60 + 310 * index
        d.rect(x, 300, 260, 139, color, fill=FILLS[family])
        d.ax.plot([x + 16, x + 244], [430, 430], color=color, linewidth=6, solid_capstyle="round")
        d.text(x + 130, 316, stage, size=20, color=family, weight="bold", ha="center")
        d.text(x + 130, 357, first, size=23, color=family, weight="bold", ha="center")
        d.text(x + 130, 389, second, size=23, color=family, weight="bold", ha="center")
        if index:
            d.arrow([(x - 50, 369), (x, 369)])
        if index < 6:
            d.arrow([(x + 130, 439), (x + 130, 483)], color)
            d.text(x + 130, 495, "Assigned here", size=22, color=family, ha="center")
    d.note(553, "Host stages use the assigned host only. Other-sample matches indicate sequence sharing, not proof of origin.")
    d.note(595, "Preserve alternative / discordant / non-unique evidence and single mates; final unexplained reads remain available.")
    d.save(directory, "phage_isolates_06_read_accounting")


def reports(directory, cfg):
    d = Diagram("07 / REPORTS AND REVIEW", "Recovery/purity decisions and QC warnings are separate assessments", 900)
    d.card(60, 160, 500, "PASS", ["Recovery/purity checks pass", "Not proof of a pure single-virus isolate"], TEAL, size=22)
    d.card(650, 160, 500, "REVIEW", ["Recovery or purity concern", "Keep genome and measured reasons"], AMBER, size=22)
    d.card(1240, 160, 500, "FAIL", ["No QC reads / no assembled contigs", "Missing or invalid read accounting"], RED, size=22)
    d.card(60, 385, 780, "Quantitative review triggers", [
        "Retained-contig count ≠1",
        f"Own-retained mapping <{number(cfg['isolate_min_retained_mapping_percent'])}% of QC-passed reads",
        f"Host-associated reads >{number(cfg['isolate_max_host_mapping_percent'])}% (chromosome + viral regions)",
        f"Unexplained reads >{number(cfg['isolate_max_unexplained_percent'])}% of QC-passed reads"], AMBER)
    d.card(960, 385, 780, "Also requires review", [
        "Coverage-supported host-prophage genome excluded",
        "No assessed host reference; host origin remains unknown",
        "Optional host test: alternative best host",
        "Optional host test: ambiguous best host"], AMBER)
    d.band(665, "MAIN OUTPUTS", "08_phage_isolates_summary.tot.csv / .html · executed notebooks · all_contig_metadata.tsv · cluster_metadata.tsv", MUTED)
    d.note(790, "Figures: mapping heatmap (% and reference counts) · cluster/sample breadth · Complete representative plots")
    d.note(832, "QC warnings stay separate. Host activity is supporting evidence; neither automatically rejects a recovered genome.")
    d.save(directory, "phage_isolates_07_reports")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "figures")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (WORKFLOW / "config.yaml").open() as handle:
        cfg = yaml.safe_load(handle)
    for draw in [overview, read_qc, assembly, hosts, selection, clustering, accounting, reports]:
        draw(args.output_dir, cfg)
    print(f"Wrote overview + seven step diagrams (PNG and editable SVG) to {args.output_dir}")


if __name__ == "__main__":
    main()
