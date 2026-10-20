"""Command line interface for PUFprobe.

Eight pipeline subcommands (prepare, embed, fold, search, context, rank,
report, annotate) share global options for the working directory, threads, CPU
fallback, configuration and logging.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

import click

from pufprobe import __version__
from pufprobe.config import load_config
from pufprobe.workdir import WorkDir


class Context:
    def __init__(self, workdir: WorkDir, config: Dict[str, Any]):
        self.workdir = workdir
        self.config = config


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="pufprobe")
@click.option(
    "--workdir",
    "-w",
    required=False,
    type=click.Path(file_okay=False, path_type=Path),
    envvar="PUFPROBE_WORKDIR",
    help="Working directory for intermediate and final outputs.",
)
@click.option("--threads", "-t", default=1, show_default=True, help="Threads for external tools.")
@click.option("--cpu", is_flag=True, default=False, help="Force CPU fallback models.")
@click.option("--config", "config_path", type=click.Path(exists=True, dir_okay=False, path_type=Path), help="JSON config overriding defaults.")
@click.option("--log-level", default="INFO", show_default=True, help="Logging level.")
@click.pass_context
def cli(
    ctx: click.Context,
    workdir: Path,
    threads: int,
    cpu: bool,
    config_path: Optional[Path],
    log_level: str,
) -> None:
    """Second-stage functional inference for hypothetical proteins."""
    import sys

    help_requested = any(arg in ("-h", "--help") for arg in sys.argv[1:])
    if workdir is None:
        if help_requested:
            ctx.obj = None
            return
        if ctx.invoked_subcommand == "doctor":
            ctx.obj = Context(None, load_config(config_path))
            return
        raise click.UsageError("Missing option '--workdir' / '-w'.")
    config = load_config(config_path)
    wd = WorkDir(workdir, cpu=cpu, threads=threads, log_level=log_level)
    wd.write_metadata()
    ctx.obj = Context(wd, config)


def _resume_guard(wd: WorkDir, stage: str, force: bool) -> bool:
    if wd.is_complete(stage) and not force:
        click.echo(f"[skip] {stage} already complete (use --force to rerun).")
        return True
    return False


@cli.command()
@click.option("--input", "inputs", multiple=True, type=click.Path(exists=True, path_type=Path), help="Input GFF3/GBK/FAA file or directory; repeatable.")
@click.option("--input-list", type=click.Path(exists=True, dir_okay=False, path_type=Path), help="Text file with one input path per line.")
@click.option("--include-smorf", is_flag=True, default=False, help="Keep sub-30 aa candidates via the permissive smORF path.")
@click.option("--force", is_flag=True, default=False, help="Ignore the completion marker.")
@click.pass_obj
def prepare(obj: Context, inputs: tuple, input_list: Optional[Path], include_smorf: bool, force: bool) -> None:
    """Parse Bakta/Prokka outputs and run true-gene quality control."""
    from pufprobe.prepare.pipeline import run

    if _resume_guard(obj.workdir, "prepare", force):
        return
    paths = _collect_inputs(list(inputs), input_list)
    if not paths:
        raise click.UsageError("Provide at least one --input path or --input-list.")
    run(obj.workdir, obj.config, paths, include_smorf=include_smorf)


@cli.command()
@click.option("--model", default=None, help="ESM-2 checkpoint name (defaults from config).")
@click.option("--precomputed", type=click.Path(exists=True, path_type=Path), help="Directory or .npz with precomputed embeddings.")
@click.option("--fallback-only", is_flag=True, default=False, help="Use sequence feature fallback without torch.")
@click.option("--force", is_flag=True, default=False)
@click.pass_obj
def embed(obj: Context, model: Optional[str], precomputed: Optional[Path], fallback_only: bool, force: bool) -> None:
    """Compute embeddings, cluster families and select representatives."""
    from pufprobe.embed.pipeline import run

    if _resume_guard(obj.workdir, "embed", force):
        return
    run(obj.workdir, obj.config, model=model, precomputed=precomputed, fallback_only=fallback_only)


@cli.command()
@click.option("--engine", type=click.Choice(["esmfold", "af2"]), default=None, help="Structure prediction engine.")
@click.option("--precomputed", type=click.Path(exists=True, file_okay=False, path_type=Path), help="Directory of PDB files named by family id.")
@click.option("--force", is_flag=True, default=False)
@click.pass_obj
def fold(obj: Context, engine: Optional[str], precomputed: Optional[Path], force: bool) -> None:
    """Fold family representatives and apply the pLDDT quality gate."""
    from pufprobe.fold.pipeline import run

    if _resume_guard(obj.workdir, "fold", force):
        return
    run(obj.workdir, obj.config, engine=engine, precomputed=precomputed)


@cli.command()
@click.option("--channels", default=None, help="Comma-separated channels: mmseqs,foldseek,hhblits,tmvec.")
@click.option("--dbdir", type=click.Path(exists=True, file_okay=False, path_type=Path), help="Directory holding versioned search databases.")
@click.option("--force", is_flag=True, default=False)
@click.pass_obj
def search(obj: Context, channels: Optional[str], dbdir: Optional[Path], force: bool) -> None:
    """Run structure/sequence remote-homology searches into unified evidence."""
    from pufprobe.search.pipeline import run

    if _resume_guard(obj.workdir, "search", force):
        return
    channel_list = [c.strip() for c in channels.split(",")] if channels else None
    run(obj.workdir, obj.config, channels=channel_list, dbdir=dbdir)


@cli.command()
@click.option("--card", type=click.Path(exists=True, path_type=Path), help="CARD reference (FASTA or ARO table) for neighbourhood checks.")
@click.option("--vfdb", type=click.Path(exists=True, path_type=Path), help="VFDB reference for neighbourhood checks.")
@click.option("--reference-fna", "reference_fna", type=click.Path(exists=True, file_okay=False, path_type=Path), help="Directory of genome FASTA files for micro-synteny.")
@click.option("--force", is_flag=True, default=False)
@click.pass_obj
def context(obj: Context, card: Optional[Path], vfdb: Optional[Path], reference_fna: Optional[Path], force: bool) -> None:
    """Compute operon, phylogenetic-profile and micro-synteny evidence."""
    from pufprobe.context.pipeline import run

    if _resume_guard(obj.workdir, "context", force):
        return
    run(obj.workdir, obj.config, card=card, vfdb=vfdb, reference_fna=reference_fna)


@cli.command()
@click.option("--weights", type=click.Path(exists=True, path_type=Path), help="Trained fusion model weights; absent uses rule weights.")
@click.option("--force", is_flag=True, default=False)
@click.pass_obj
def rank(obj: Context, weights: Optional[Path], force: bool) -> None:
    """Fuse evidence, calibrate confidence and rank wet-lab priorities."""
    from pufprobe.rank.pipeline import run

    if _resume_guard(obj.workdir, "rank", force):
        return
    run(obj.workdir, obj.config, weights=weights)


@cli.command()
@click.option("--formats", default="json,html,md", show_default=True, help="Comma-separated output formats.")
@click.option("--force", is_flag=True, default=False)
@click.pass_obj
def report(obj: Context, formats: str, force: bool) -> None:
    """Write JSON and auditable HTML/Markdown evidence reports."""
    from pufprobe.report.pipeline import run

    if _resume_guard(obj.workdir, "report", force):
        return
    run(obj.workdir, obj.config, formats=[f.strip() for f in formats.split(",")])


@cli.command()
@click.option("--min-confidence", default=None, type=float, help="Override the annotation confidence threshold.")
@click.option("--force", is_flag=True, default=False)
@click.pass_obj
def annotate(obj: Context, min_confidence: Optional[float], force: bool) -> None:
    """Backfill high-confidence predictions into new GFF3/GBK/FAA files."""
    from pufprobe.annotate.pipeline import run

    if _resume_guard(obj.workdir, "annotate", force):
        return
    run(obj.workdir, obj.config, min_confidence=min_confidence)


@cli.command()
@click.pass_obj
def doctor(obj: Context) -> None:
    """Report detected external tools, Python modules and GPU availability."""
    from pufprobe.externals import summary
    import json

    click.echo(json.dumps(summary(), indent=2, ensure_ascii=False))


def _collect_inputs(inputs: List[Path], input_list: Optional[Path]) -> List[Path]:
    paths = list(inputs)
    if input_list is not None:
        for line in input_list.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(Path(line))
    expanded: List[Path] = []
    for path in paths:
        if not path.exists():
            raise click.UsageError(f"Input path does not exist: {path}")
        if path.is_dir():
            for pattern in ("*.gff", "*.gff3"):
                expanded.extend(sorted(path.rglob(pattern)))
            for pattern in ("*.gbk", "*.gbff", "*.genbank"):
                expanded.extend(sorted(path.rglob(pattern)))
            for faa in sorted(path.rglob("*.faa")):
                # Bakta also writes <prefix>.hypothetical.faa, a subset holding
                # only hypothetical proteins. Indexing that subset would make
                # every annotated protein appear unaligned, so skip it and keep
                # the complete .faa.
                if faa.name.endswith(".hypothetical.faa"):
                    continue
                expanded.append(faa)
        else:
            expanded.append(path)
    return sorted(set(expanded))


if __name__ == "__main__":
    cli()
