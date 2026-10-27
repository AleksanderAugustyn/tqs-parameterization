"""Locate and load libtqs_parameterization.so for ctypes."""
from __future__ import annotations

import ctypes
import os
from pathlib import Path
from typing import Iterator, Optional

_LIB_NAME = "libtqs_parameterization.so"


def _candidate_paths() -> Iterator[Path]:
    env = os.environ.get("TQS_PARAM_LIB")
    if env:
        yield Path(env)
    # Installed-wheel location: the .so is bundled next to this file.
    yield Path(__file__).resolve().parent / _LIB_NAME
    # Dev fallback: this file is <repo>/python/tqs_parameterization/_libloader.py
    repo_root = Path(__file__).resolve().parents[2]
    yield repo_root / "build" / "release" / _LIB_NAME
    for build_dir in ("cmake-build-release", "cmake-build-debug"):
        yield repo_root / build_dir / _LIB_NAME
    # System search path
    yield Path(_LIB_NAME)


def load_library() -> ctypes.CDLL:
    """Load the shared library: TQS_PARAM_LIB, the wheel, build dirs, then the system path."""
    last_error: Optional[OSError] = None
    for path in _candidate_paths():
        try:
            return ctypes.CDLL(str(path))
        except OSError as exc:
            last_error = exc
    raise OSError(
        f"Could not load {_LIB_NAME}. Build it (cmake --build build/release) or set "
        f"TQS_PARAM_LIB to its full path. Last error: {last_error}"
    )
