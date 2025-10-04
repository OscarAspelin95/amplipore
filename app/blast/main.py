from .config import BlastnConfig
from .parse import parse_blast_tsv
from sh import blastn, makeblastdb, RunningCommand
from pathlib import Path
import logging
from common.decorator import with_yaspin

log = logging.getLogger(__name__)


def run_makeblastdb(asv_fasta: Path, outdir: Path) -> Path:
    nucl_db = outdir / "nucl_db"

    makeblastdb("-in", asv_fasta, "-dbtype", "nucl", "-out", nucl_db)
    return nucl_db


def run_blastn(database: Path, nucl_db: Path, outdir: Path, cfg: BlastnConfig) -> Path:
    blast_tsv = outdir / "blast.tsv"
    rc: RunningCommand = blastn(
        "-query",
        database,
        "-db",
        nucl_db,
        "-num_threads",
        cfg.threads,
        "-qcov_hsp_perc",
        cfg.qcov,
        "-word_size",
        cfg.word_size,
        "-perc_identity",
        cfg.pident,
        "-out",
        blast_tsv,
        "-outfmt",
        "6 qseqid qlen qstart qend qframe sseqid slen sstart send sframe length pident",
        _return_cmd=True,
    )

    # Might never reach this point if exit_code != 0.
    if rc.exit_code != 0:
        log.error(rc.stderr)

    assert blast_tsv.is_file()
    return blast_tsv


@with_yaspin(progress_text="Running BLASTn classification...")
def run_blast(asv_fasta: Path, db_fasta: Path, outdir: Path) -> Path:
    # Create sub-directory for BLAST results.
    outdir = outdir / "blast"
    outdir.mkdir(exist_ok=True)

    cfg = BlastnConfig()

    nucl_db = run_makeblastdb(asv_fasta, outdir)

    blast_tsv = run_blastn(db_fasta, nucl_db, outdir, cfg)

    return parse_blast_tsv(blast_tsv)
