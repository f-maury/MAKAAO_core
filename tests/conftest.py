from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def pytest_sessionstart(session: pytest.Session) -> None:
    """Resolve repository-relative test assets independently of launch CWD."""
    os.chdir(PROJECT_ROOT)


@pytest.fixture(scope="session")
def sample_core_path() -> Path:
    """Return the committed, license-safe MAKAAO sample CSV."""
    configured = os.environ.get("MAKAAO_SAMPLE_DATA")
    path = Path(configured).expanduser() if configured else Path("data/makaao_sample.csv")
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    assert path.is_file(), (
        f"Required public sample fixture is missing: {path}. "
        "Commit data/makaao_sample.csv or set MAKAAO_SAMPLE_DATA."
    )
    assert path.stat().st_size > 0, f"Sample fixture is empty: {path}"
    return path


@pytest.fixture(scope="session")
def sample_core_frame(sample_core_path: Path) -> pd.DataFrame:
    frame = pd.read_csv(sample_core_path, dtype=str, keep_default_na=False, low_memory=False)
    assert not frame.empty, f"Sample fixture contains no rows: {sample_core_path}"
    assert any(str(column).strip() for column in frame.columns), (
        f"Sample fixture contains no usable columns: {sample_core_path}"
    )
    return frame