from enum import Enum, unique


@unique
class Levels(Enum):
    SPECIES = "species"
    GENUS = "genus"
    FAMILY = "family"
    ORDER = "order"
    CLASS = "class"
    PHYLUM = "phylum"
    KINGDOM = "kingdom"

    @classmethod
    def as_list(cls) -> list[str]:
        return [c.value for c in cls]

    @classmethod
    def pandas_column_rename(cls) -> dict[int, str]:
        return {i: c for i, c in enumerate(cls.as_list()[::-1])}


TAXONOMY_PREFIXES = {
    "d": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
