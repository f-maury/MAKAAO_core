import sys
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path("scripts/01_process_makaao_core_to_tables.py")
pytestmark = pytest.mark.skipif(not SCRIPT_PATH.is_file(), reason="scripts/01_process_makaao_core_to_tables.py is not present")



def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_tokenizers_and_normalizers():
    mod = load_module(
        SCRIPT_PATH, "process_core"
    )

    assert mod.split_tokens("a; b ;c") == ["a", "b", "c"]
    assert mod.unique_preserve_order(["a", "b", "a", "c", "b"]) == ["a", "b", "c"]
    assert mod.clean_cui("C1234567") == "C1234567"
    assert mod.clean_cui(" UMLS:C1234567 ") == "UMLS:C1234567"
    assert mod.norm_hp("HP:0003731") == "hp:0003731"
    assert mod.to_int_or_none("12") == 12
    assert mod.to_int_or_none("") is None


def test_load_core_and_write_some_outputs(tmp_path, sample_core_path):
    mod = load_module(
        SCRIPT_PATH, "process_core_outputs"
    )

    frame = mod.load_core(sample_core_path)
    assert "index" in frame.columns
    assert pd.api.types.is_integer_dtype(frame["index"])

    output_dir = tmp_path / "processed_tables"
    output_dir.mkdir()
    mod.write_index_name_en(frame, output_dir)
    mod.write_index_hpo_id(frame, output_dir)
    mod.write_index_parent_index(frame, output_dir)

    for filename in (
        "index_name_en.csv",
        "index_hpo_id.csv",
        "index_parent_index.csv",
    ):
        path = output_dir / filename
        assert path.is_file()
        with path.open("r", encoding="utf-8") as handle:
            assert sum(1 for _ in handle) >= 2