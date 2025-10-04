from pydantic import BaseModel


class BlastnConfig(BaseModel):
    """Settings for BLASTn."""

    threads: int = 8
    pident: int = 90
    qcov: int = 90
    word_size: int = 11


class ParseConfig(BaseModel):
    """Settings for parsing BLASTn results."""

    pident: float = 0.90
    perc_aln: float = 0.90
