from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture(scope="session")
def sample_core_path() -> Path:
    """Return the committed, license-safe MAKAAO sample CSV."""
    path = Path(os.environ.get("MAKAAO_SAMPLE_DATA", "data/makaao_sample.csv"))
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
