from yaspin import yaspin
from typing import Callable, TypeVar, ParamSpec
from functools import wraps

T = TypeVar("Type")
P = ParamSpec("ParamSpec")


def with_yaspin(progress_text: str, color: str = "cyan") -> Callable[P, T]:
    """Function decorator that adds a progress spinner."""

    def with_progress(func: Callable[P, T]) -> Callable[P, T]:
        @wraps(func)
        def inner(*args: P.args, **kwargs: P.kwargs) -> T:
            with yaspin(text=progress_text, color=color) as sp:
                result = func(*args, **kwargs)
                sp.ok("✔")

            return result

        return inner

    return with_progress
