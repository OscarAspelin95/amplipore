from pathlib import Path


class InvalidFileEndingError(Exception):
    pass


def get_file_base(f: Path, allowed_file_endings: tuple[str]):
    f_name = f.name
    for ending in allowed_file_endings:
        if f_name.endswith(ending):
            return f_name.removesuffix(ending)

    msg = f"{f_name} has an invalid file ending. Must be {', '.join(allowed_file_endings)}"
    raise InvalidFileEndingError(msg)


def _file(f: str, allowed_file_endings: tuple[str] | None) -> Path:
    if not (f_path := Path(f)).is_file():
        raise FileNotFoundError(f)

    if allowed_file_endings and not f.endswith(allowed_file_endings):
        msg = f"File {f} has an invalid file ending. Must be {', '.join(allowed_file_endings)}"
        raise InvalidFileEndingError(msg)

    return f_path
