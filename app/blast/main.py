from .config import BlastnConfig
from .parse import parse_blast_tsv
from sh import blastn, makeblastdb, RunningCommand
from pathlib import Path
import logging
import tempfile
import pandas as pd
from common.decorator import with_yaspin

log = logging.getLogger(__name__)


def run_makeblastdb(db_fasta: Path, db_dir: Path) -> Path:
    nucl_db = db_dir / "nucl_db"

    makeblastdb("-in", db_fasta, "-dbtype", "nucl", "-out", nucl_db)
    return nucl_db


def run_blastn(asv_fasta: Path, nucl_db: Path, outdir: Path, cfg: BlastnConfig) -> Path:
    blast_tsv = outdir / "blast.tsv"
    rc: RunningCommand = blastn(
        "-query",
        asv_fasta,
        "-db",
        nucl_db,
        "-num_threads",
        cfg.threads,
        "-qcov_hsp_perc",
        int(cfg.qcov * 100),
        "-word_size",
        cfg.word_size,
        "-perc_identity",
        int(cfg.pident * 100),
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
def run_blast(
    asv_fasta: Path,
    db_fasta: Path,
    outdir: Path,
    blast_pident: float = 0.70,
    blast_qcov: float = 0.70,
) -> pd.DataFrame:
    cfg = BlastnConfig(pident=blast_pident, qcov=blast_qcov)

    # Build BLAST DB in a temp directory (cleaned up automatically).
    with tempfile.TemporaryDirectory() as tmpdir:
        nucl_db = run_makeblastdb(db_fasta, Path(tmpdir))
        blast_tsv = run_blastn(asv_fasta, nucl_db, outdir, cfg)

    return parse_blast_tsv(blast_tsv, blast_pident, blast_qcov)
