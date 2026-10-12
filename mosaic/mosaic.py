from __future__ import annotations

import logging
import multiprocessing
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Tuple

import click

try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:  # pragma: no cover - Python <3.8 fallback
    from importlib_metadata import PackageNotFoundError, version

try:
    from .config import set_logger
except ImportError:  # Allows `python mosaic.py ...` from a source checkout.
    from config import set_logger


def get_version() -> str:
    try:
        return version("mosaic")
    except PackageNotFoundError:
        return "0+local"


__version__ = get_version()

WORKFLOW_DIR = Path(__file__).resolve().parent
DEFAULT_SNAKEFILE = WORKFLOW_DIR / "Snakefile"

DEFAULT_FLAGS: Dict[str, object] = {
    "subsampling": False,
    "subassembly": False,
    "assembly_stats": False,
    "cross_assembly": False,
    "long_index": False,
    "run_vcontact": False,
    "run_DRAM": False,
    "metavr_blast": False,
    "extract_mapped": False,
    "visualization_tool": "lovis4u",
    "visualization_max_contigs": 50,
}

PRESETS: Dict[str, Dict[str, object]] = {
    "rna_viruses": {
        **DEFAULT_FLAGS,
        "metagenome": True,
        "isolates": False,
        "microbial": False,
        "sourmash_clean_reads": False,
        "remove_euk": False,
        "RNA_enriched": True,
    },
    "viral_metagenome": {
        **DEFAULT_FLAGS,
        "metagenome": True,
        "isolates": False,
        "microbial": False,
        "sourmash_clean_reads": False,
        "remove_euk": True,
    },
    "phage_isolates": {
        **DEFAULT_FLAGS,
        "metagenome": False,
        "isolates": True,
        "microbial": False,
        "sourmash_clean_reads": False,
        "remove_euk": False,
    },
    "microbial_metagenome": {
        **DEFAULT_FLAGS,
        "metagenome": True,
        "isolates": False,
        "microbial": True,
        "sourmash_clean_reads": True,
        "remove_euk": True,
    },
    "bacteria_isolate": {
        **DEFAULT_FLAGS,
        "metagenome": False,
        "isolates": True,
        "microbial": True,
        "sourmash_clean_reads": False,
        "remove_euk": False,
    },
    "fungi_isolate": {
        **DEFAULT_FLAGS,
        "metagenome": False,
        "isolates": True,
        "microbial": True,
        "sourmash_clean_reads": False,
        "remove_euk": False,
    },
}

MODES: Dict[str, Dict[str, str]] = {
    "rna_viruses": {
        "target": "runRNA",
        "preset": "rna_viruses",
        "reads": "illumina",
        "description": "RNA assembly, viral candidates, and optional microdiversity",
    },
    "rna_assembly": {
        "target": "runRNAAssembly",
        "preset": "rna_viruses",
        "reads": "illumina",
        "description": "RNA assembly with selected assemblers",
    },
    "microdiversity": {
        "target": "runMicrodiversity",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "Within-sample nucleotide variation against the final vOTU catalogue",
    },
    "refseq_mapping": {
        "target": "map_to_RefSeq",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "Cleaned paired-end read mapping and abundance against RefSeq Viral",
    },
    "qc": {
        "target": "runQC",
        "preset": "viral_metagenome",
        "reads": "any",
        "description": "QC only",
    },
    "only_qc": {
        "target": "runQC",
        "preset": "viral_metagenome",
        "reads": "any",
        "description": "QC only",
    },
    "assembly": {
        "target": "runAssembly",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "QC plus short-read assembly",
    },
    "viral_id": {
        "target": "runViralID",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "QC, assembly, and viral identification",
    },
    "votu": {
        "target": "runvOTUClustering",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "Viral identification and vOTU clustering",
    },
    "abundance": {
        "target": "runAbundance",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "Read mapping and abundance tables",
    },
    "annotation": {
        "target": "runAnnotation",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "Viral annotation, taxonomy, and host prediction",
    },
    "end_to_end": {
        "target": "runWorkflow",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "Viral metagenome workflow",
    },
    "viral_metagenome": {
        "target": "runWorkflow",
        "preset": "viral_metagenome",
        "reads": "illumina",
        "description": "Alias for end_to_end",
    },
    "phage_isolates": {
        "target": "phage_isolates",
        "preset": "phage_isolates",
        "reads": "illumina",
        "description": "Phage isolate workflow",
    },
    "microbial_metagenome": {
        "target": "microbial_metagenome",
        "preset": "microbial_metagenome",
        "reads": "illumina",
        "description": "Microbial metagenome workflow",
    },
    "short_read_bacteria": {
        "target": "assembly_short_bacteria",
        "preset": "bacteria_isolate",
        "reads": "illumina",
        "description": "Short-read bacterial isolate workflow",
    },
    "short_read_fungi": {
        "target": "assembly_short_fungi",
        "preset": "fungi_isolate",
        "reads": "illumina",
        "description": "Short-read fungal isolate workflow",
    },
    "nanopore_assembly": {
        "target": "assembly_nanopore",
        "preset": "viral_metagenome",
        "reads": "nanopore",
        "description": "Nanopore-only assembly",
    },
    "nanopore_bacteria": {
        "target": "assembly_nanopore_only_bacteria",
        "preset": "bacteria_isolate",
        "reads": "nanopore",
        "description": "Nanopore-only bacterial isolate workflow",
    },
    "hybrid_nanopore_bacteria": {
        "target": "assembly_nanopore_hybrid_bacteria",
        "preset": "bacteria_isolate",
        "reads": "illumina+nanopore",
        "description": "Illumina plus Nanopore bacterial isolate workflow",
    },
    "long_read_phage": {
        "target": "assembly_long_only_phage",
        "preset": "phage_isolates",
        "reads": "nanopore",
        "description": "Long-read phage isolate workflow",
    },
    "pacbio_bacteria": {
        "target": "assembly_pacbio_only_bacteria",
        "preset": "bacteria_isolate",
        "reads": "pacbio",
        "description": "PacBio-only bacterial isolate workflow",
    },
    "hybrid_pacbio_bacteria": {
        "target": "assembly_pacbio_hybrid_bacteria",
        "preset": "bacteria_isolate",
        "reads": "illumina+pacbio",
        "description": "Illumina plus PacBio bacterial isolate workflow",
    },
    "pacbio_fungi": {
        "target": "assembly_pacbio_only_fungi",
        "preset": "fungi_isolate",
        "reads": "pacbio",
        "description": "PacBio-only fungal isolate workflow",
    },
    "hybrid_pacbio_fungi": {
        "target": "assembly_pacbio_hybrid_fungi",
        "preset": "fungi_isolate",
        "reads": "illumina+pacbio",
        "description": "Illumina plus PacBio fungal isolate workflow",
    },
}

MODE_ALIASES = {
    "onlyqc": "qc",
    "only-qc": "qc",
    "end-to-end": "end_to_end",
    "viral-id": "viral_id",
    "short-read-bacteria": "short_read_bacteria",
    "short-read-fungi": "short_read_fungi",
    "nanopore-assembly": "nanopore_assembly",
    "nanopore-bacteria": "nanopore_bacteria",
    "hybrid-nanopore-bacteria": "hybrid_nanopore_bacteria",
    "long-read-phage": "long_read_phage",
    "pacbio-bacteria": "pacbio_bacteria",
    "hybrid-pacbio-bacteria": "hybrid_pacbio_bacteria",
    "pacbio-fungi": "pacbio_fungi",
    "hybrid-pacbio-fungi": "hybrid_pacbio_fungi",
    "phage-isolates": "phage_isolates",
    "microbial-metagenome": "microbial_metagenome",
    "viral-metagenome": "viral_metagenome",
}

STAGE_MODES = {
    "qc",
    "only_qc",
    "assembly",
    "viral_id",
    "votu",
    "abundance",
    "annotation",
}


def normalize_name(name: str) -> str:
    lowered = name.strip().lower()
    return MODE_ALIASES.get(lowered, lowered.replace("-", "_"))


def format_config_value(value: object) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    return str(value)


def quote_command(command: Iterable[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def parse_key_value(items: Iterable[str]) -> Dict[str, str]:
    parsed = {}
    for item in items:
        if "=" not in item:
            raise click.ClickException(
                f"Invalid --config value '{item}'. Use KEY=VALUE."
            )
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise click.ClickException(
                f"Invalid --config value '{item}'. The key is empty."
            )
        parsed[key] = value
    return parsed


def read_tags(config_values: Mapping[str, object]) -> Tuple[str, str, str, str]:
    return (
        str(config_values.get("forward_tag", "R1")),
        str(config_values.get("reverse_tag", "R2")),
        str(config_values.get("nanopore_tag", "nanopore")),
        str(config_values.get("pacbio_tag", "pacbio")),
    )


def discover_reads(
    raw_dir: Path,
    forward_tag: str,
    reverse_tag: str,
    nanopore_tag: str,
    pacbio_tag: str,
) -> Dict[str, object]:
    forward_suffix = f"_{forward_tag}.fastq.gz"
    reverse_suffix = f"_{reverse_tag}.fastq.gz"

    forward_samples = {
        path.name[: -len(forward_suffix)]
        for path in raw_dir.glob(f"*{forward_suffix}")
    }
    reverse_samples = {
        path.name[: -len(reverse_suffix)]
        for path in raw_dir.glob(f"*{reverse_suffix}")
    }
    nanopore_files = list(raw_dir.glob(f"*_{nanopore_tag}.fastq.gz"))
    pacbio_files = list(raw_dir.glob(f"*_{pacbio_tag}.fastq.gz"))

    return {
        "illumina_samples": sorted(forward_samples & reverse_samples),
        "forward_only": sorted(forward_samples - reverse_samples),
        "reverse_only": sorted(reverse_samples - forward_samples),
        "nanopore_count": len(nanopore_files),
        "pacbio_count": len(pacbio_files),
    }


def validate_raw_dir(raw_dir: Path, requirement: str, config_values: Mapping[str, object]) -> None:
    if not raw_dir.exists():
        raise click.ClickException(f"Raw data folder does not exist: {raw_dir}")
    if not raw_dir.is_dir():
        raise click.ClickException(f"Raw data path is not a folder: {raw_dir}")

    forward_tag, reverse_tag, nanopore_tag, pacbio_tag = read_tags(config_values)
    reads = discover_reads(raw_dir, forward_tag, reverse_tag, nanopore_tag, pacbio_tag)
    has_illumina = bool(reads["illumina_samples"])
    has_nanopore = int(reads["nanopore_count"]) > 0
    has_pacbio = int(reads["pacbio_count"]) > 0

    if reads["forward_only"] or reads["reverse_only"]:
        logging.warning(
            "Found unmatched Illumina files. Forward-only samples: %s. "
            "Reverse-only samples: %s.",
            ", ".join(reads["forward_only"]) or "none",
            ", ".join(reads["reverse_only"]) or "none",
        )

    checks = {
        "illumina": has_illumina,
        "nanopore": has_nanopore,
        "pacbio": has_pacbio,
    }

    if requirement == "any":
        if has_illumina or has_nanopore or has_pacbio:
            return
        raise click.ClickException(
            "No recognized raw FASTQ files found. Expected files like "
            f"*_{forward_tag}.fastq.gz / *_{reverse_tag}.fastq.gz, "
            f"*_{nanopore_tag}.fastq.gz, or *_{pacbio_tag}.fastq.gz."
        )

    missing = [part for part in requirement.split("+") if not checks.get(part, False)]
    if missing:
        raise click.ClickException(
            "Raw data folder is missing required read type(s) for this mode: "
            + ", ".join(missing)
        )


def set_optional_bool(
    config_values: MutableMapping[str, object],
    key: str,
    value: Optional[bool],
) -> None:
    if value is not None:
        config_values[key] = value


def build_command(
    target: str,
    config_values: Mapping[str, object],
    jobs: int,
    snakefile: Path,
    workflow_dir: Path,
    use_conda: bool,
    printshellcmds: bool,
    keep_going: bool,
    rerun_incomplete: bool,
    dry_run: bool,
    reason: bool,
    latency_wait: int,
    show_failed_logs: bool,
    profile: Optional[str],
    extra_args: Tuple[str, ...],
) -> List[str]:
    command = [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--directory",
        str(workflow_dir),
    ]
    if use_conda:
        command.append("--use-conda")
    if printshellcmds:
        command.append("-p")
    if profile:
        command.extend(["--profile", profile])

    command.append(target)

    if config_values:
        command.append("--config")
        command.extend(
            f"{key}={format_config_value(value)}"
            for key, value in config_values.items()
        )

    command.extend(["--jobs", str(jobs)])

    if keep_going:
        command.append("--keep-going")
    if rerun_incomplete:
        command.append("--rerun-incomplete")
    if latency_wait:
        command.extend(["--latency-wait", str(latency_wait)])
    if dry_run:
        command.append("--dry-run")
    if reason:
        command.append("--reason")
    if show_failed_logs:
        command.append("--show-failed-logs")
    command.extend(extra_args)
    return command


def build_unlock_command(snakefile: Path, workflow_dir: Path) -> List[str]:
    return [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--directory",
        str(workflow_dir),
        "--unlock",
    ]


def run_snakemake_command(command: List[str]) -> None:
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError:
        raise click.ClickException("Could not find 'snakemake' on PATH.")
    except subprocess.CalledProcessError as exc:
        raise SystemExit(exc.returncode)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__)
def cli() -> None:
    """mosaic - modular viral and microbial genomics workflows."""


@cli.command("modes")
def list_modes() -> None:
    """Show wrapper modes and their Snakemake targets."""
    click.echo("Modes:")
    for name in sorted(MODES):
        mode = MODES[name]
        click.echo(
            f"  {name:28} -> {mode['target']:30} "
            f"preset={mode['preset']} reads={mode['reads']}"
        )


@cli.command("unlock")
@click.option(
    "--snakefile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_SNAKEFILE,
    show_default=True,
    help="Snakefile whose Snakemake lock should be removed.",
)
@click.option(
    "--workflow-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=WORKFLOW_DIR,
    show_default=True,
    help="Workflow working directory containing the .snakemake lock folder.",
)
@click.option(
    "--print-only",
    is_flag=True,
    help="Print the resolved Snakemake unlock command and exit.",
)
def unlock(snakefile: Path, workflow_dir: Path, print_only: bool) -> None:
    """Remove a stale Snakemake lock for the MOSAIC working directory."""
    set_logger()
    logging.info("mosaic %s", __version__)
    logging.info("%s", quote_command(sys.argv))

    snakefile = snakefile.resolve()
    workflow_dir = workflow_dir.resolve()
    if not snakefile.exists():
        raise click.ClickException(f"Snakefile does not exist: {snakefile}")
    if not workflow_dir.exists():
        raise click.ClickException(f"Workflow directory does not exist: {workflow_dir}")

    command = build_unlock_command(snakefile=snakefile, workflow_dir=workflow_dir)

    click.echo("Snakemake command:")
    click.echo(quote_command(command))

    if print_only:
        return

    sys.stdout.flush()
    run_snakemake_command(command)


@cli.command(
    "run",
    context_settings={"ignore_unknown_options": True, "allow_extra_args": True},
)
@click.argument("mode")
@click.option(
    "--raw",
    "--input-dir",
    "raw_dir",
    required=True,
    type=click.Path(exists=False, file_okay=False, path_type=Path),
    help="Raw FASTQ folder. Required for every MOSAIC mode.",
)
@click.option(
    "--results",
    "results_dir",
    type=click.Path(file_okay=False, path_type=Path),
    help="Optional results directory. If omitted, MOSAIC infers it from --raw.",
)
@click.option(
    "--preset",
    help=(
        "Preset for stage modes. Choices: "
        + ", ".join(sorted(PRESETS))
        + ". Direct biological modes choose this automatically."
    ),
)
@click.option(
    "-j",
    "--jobs",
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help="Maximum number of Snakemake jobs.",
)
@click.option(
    "-n",
    "--dry-run",
    "dry_run",
    is_flag=True,
    help="Show the DAG and commands without executing jobs.",
)
@click.option(
    "--print-only",
    is_flag=True,
    help="Print the resolved Snakemake command and exit.",
)
@click.option(
    "--use-conda/--no-use-conda",
    default=True,
    show_default=True,
    help="Pass --use-conda to Snakemake.",
)
@click.option(
    "--printshellcmds/--no-printshellcmds",
    default=True,
    show_default=True,
    help="Print shell commands from Snakemake.",
)
@click.option(
    "-k",
    "--keep-going/--no-keep-going",
    default=True,
    show_default=True,
    help="Keep running independent jobs after failures.",
)
@click.option(
    "--rerun-incomplete/--no-rerun-incomplete",
    default=True,
    show_default=True,
    help="Ask Snakemake to rerun outputs marked incomplete.",
)
@click.option(
    "--latency-wait",
    default=600,
    type=int,
    show_default=True,
    help="Snakemake latency wait in seconds.",
)
@click.option(
    "--reason",
    is_flag=True,
    help="Pass --reason to Snakemake.",
)
@click.option(
    "--show-failed-logs/--no-show-failed-logs",
    default=True,
    show_default=True,
    help="Ask Snakemake to print failed job logs.",
)
@click.option(
    "--profile",
    help="Optional Snakemake profile.",
)
@click.option(
    "--snakefile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_SNAKEFILE,
    show_default=True,
    help="Snakefile to run.",
)
@click.option(
    "--workflow-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=WORKFLOW_DIR,
    show_default=True,
    help="Workflow working directory containing config.yaml, rules/, db/, and tools/.",
)
@click.option("--kraken-db", help="Set config kraken_db.")
@click.option("--ecc-memory", type=int, help="Set config ecc_memory.")
@click.option("--complexity", type=click.Choice(["low", "medium", "high", "extreme"], case_sensitive=False),
              help="Memory profile; medium by default. Does not change the biological workflow mode.")
@click.option("--sourmash-clean-reads/--no-sourmash-clean-reads", default=None,
              help="Report per-sample and pooled bacterial/archaeal profiles from all final clean reads.")
@click.option("--remove-euk/--no-remove-euk", default=None, help="Override config remove_euk.")
@click.option("--assembly-stats/--no-assembly-stats", default=None, help="Override config assembly_stats.")
@click.option("--subsampling/--no-subsampling", default=None, help="Override config subsampling.")
@click.option("--subassembly/--no-subassembly", default=None, help="Override config subassembly.")
@click.option("--cross-assembly/--no-cross-assembly", default=None, help="Override config cross_assembly.")
@click.option("--vcontact/--no-vcontact", default=None, help="Override config run_vcontact.")
@click.option("--dram/--no-dram", default=None, help="Override config run_DRAM.")
@click.option("--metavr-blast/--no-metavr-blast", default=None, help="Override config metavr_blast.")
@click.option("--map-to-refseq/--no-map-to-refseq", default=None, help="Also map cleaned paired-end reads to RefSeqViral_db.")
@click.option("--extract-mapped/--no-extract-mapped", default=None, help="Override config extract_mapped.")
@click.option("--rna-enriched/--no-rna-enriched", default=None, help="Run additional RNA assemblers for all paired-end samples.")
@click.option("--microdiversity/--no-microdiversity", default=None, help="Calculate within-sample nucleotide variation against the final vOTU catalogue.")
@click.option(
    "--visualization-tool",
    type=click.Choice(["lovis4u", "clinker"], case_sensitive=False),
    help="Genome visualization tool for small filtered vOTU sets.",
)
@click.option(
    "--visualization-max-contigs",
    type=int,
    help="Render filtered vOTU visualization only when the FASTA has fewer contigs than this value.",
)
@click.option(
    "--config",
    "config_items",
    multiple=True,
    help="Additional Snakemake config override as KEY=VALUE. Can be repeated.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run(
    mode: str,
    raw_dir: Path,
    results_dir: Optional[Path],
    preset: Optional[str],
    jobs: int,
    dry_run: bool,
    print_only: bool,
    use_conda: bool,
    printshellcmds: bool,
    keep_going: bool,
    rerun_incomplete: bool,
    latency_wait: int,
    reason: bool,
    show_failed_logs: bool,
    profile: Optional[str],
    snakefile: Path,
    workflow_dir: Path,
    kraken_db: Optional[str],
    ecc_memory: Optional[int],
    complexity: Optional[str],
    sourmash_clean_reads: Optional[bool],
    remove_euk: Optional[bool],
    assembly_stats: Optional[bool],
    subsampling: Optional[bool],
    subassembly: Optional[bool],
    cross_assembly: Optional[bool],
    vcontact: Optional[bool],
    dram: Optional[bool],
    metavr_blast: Optional[bool],
    map_to_refseq: Optional[bool],
    extract_mapped: Optional[bool],
    rna_enriched: Optional[bool],
    microdiversity: Optional[bool],
    visualization_tool: Optional[str],
    visualization_max_contigs: Optional[int],
    config_items: Tuple[str, ...],
    snakemake_args: Tuple[str, ...],
) -> None:
    """Run a MOSAIC workflow mode."""
    set_logger()
    logging.info("mosaic %s", __version__)
    logging.info("%s", quote_command(sys.argv))

    mode_name = normalize_name(mode)
    if mode_name not in MODES:
        allowed = ", ".join(sorted(MODES))
        raise click.ClickException(f"Unknown mode '{mode}'. Available modes: {allowed}")

    mode_info = MODES[mode_name]
    preset_name = normalize_name(preset) if preset else mode_info["preset"]
    if preset_name not in PRESETS:
        allowed = ", ".join(sorted(PRESETS))
        raise click.ClickException(f"Unknown preset '{preset}'. Available presets: {allowed}")
    if preset and mode_name not in STAGE_MODES and preset_name != mode_info["preset"]:
        raise click.ClickException(
            f"Mode '{mode}' uses preset '{mode_info['preset']}'. "
            "Use stage modes such as qc/assembly/annotation if you need --preset."
        )

    config_values: Dict[str, object] = dict(PRESETS[preset_name])
    config_values["input_dir"] = str(raw_dir.resolve())
    if results_dir:
        config_values["results_dir"] = str(results_dir.resolve())
    if kraken_db:
        config_values["kraken_db"] = kraken_db
    if ecc_memory is not None:
        config_values["ecc_memory"] = ecc_memory
    if complexity is not None:
        config_values["complexity"] = complexity.lower()
    if visualization_tool:
        config_values["visualization_tool"] = visualization_tool.lower()
    if visualization_max_contigs is not None:
        if visualization_max_contigs < 0:
            raise click.ClickException("--visualization-max-contigs must be zero or greater.")
        config_values["visualization_max_contigs"] = visualization_max_contigs

    set_optional_bool(config_values, "sourmash_clean_reads", sourmash_clean_reads)
    set_optional_bool(config_values, "remove_euk", remove_euk)
    set_optional_bool(config_values, "assembly_stats", assembly_stats)
    set_optional_bool(config_values, "subsampling", subsampling)
    set_optional_bool(config_values, "subassembly", subassembly)
    set_optional_bool(config_values, "cross_assembly", cross_assembly)
    set_optional_bool(config_values, "run_vcontact", vcontact)
    set_optional_bool(config_values, "run_DRAM", dram)
    set_optional_bool(config_values, "metavr_blast", metavr_blast)
    set_optional_bool(config_values, "map_to_RefSeq", map_to_refseq)
    set_optional_bool(config_values, "extract_mapped", extract_mapped)
    set_optional_bool(config_values, "RNA_enriched", rna_enriched)
    set_optional_bool(config_values, "microdiversity", microdiversity)

    config_values.update(parse_key_value(config_items))

    raw_path = raw_dir.resolve()
    validate_raw_dir(raw_path, mode_info["reads"], config_values)

    snakefile = snakefile.resolve()
    workflow_dir = workflow_dir.resolve()
    if not snakefile.exists():
        raise click.ClickException(f"Snakefile does not exist: {snakefile}")
    if not workflow_dir.exists():
        raise click.ClickException(f"Workflow directory does not exist: {workflow_dir}")

    command = build_command(
        target=mode_info["target"],
        config_values=config_values,
        jobs=jobs,
        snakefile=snakefile,
        workflow_dir=workflow_dir,
        use_conda=use_conda,
        printshellcmds=printshellcmds,
        keep_going=keep_going,
        rerun_incomplete=rerun_incomplete,
        dry_run=dry_run,
        reason=reason,
        latency_wait=latency_wait,
        show_failed_logs=show_failed_logs,
        profile=profile,
        extra_args=snakemake_args,
    )

    click.echo("Snakemake command:")
    click.echo(quote_command(command))

    if print_only:
        return

    sys.stdout.flush()
    run_snakemake_command(command)


if __name__ == "__main__":
    cli()
