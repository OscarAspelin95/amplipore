from pydantic import BaseModel


class BlastnConfig(BaseModel):
    """Settings for BLASTn."""

    threads: int = 8
    pident: float = 0.70
    qcov: float = 0.70
    word_size: int = 11


class ParseConfig(BaseModel):
    """Settings for parsing BLASTn results."""

    pident: float = 0.70
    perc_aln: float = 0.70
