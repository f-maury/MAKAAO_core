import importlib.util
from pathlib import Path
import pandas as pd

def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod

def test_tokenizers_and_normalizers():
    mod = load_module(Path("scripts/01_process_makaao_core_to_tables.py"), "process_core")

    assert mod.split_tokens("a; b ;c") == ["a", "b", "c"]
    assert mod.unique_preserve_order(["a","b","a","c","b"]) == ["a","b","c"]

    assert mod.clean_cui("C1234567") == "C1234567"
    assert mod.clean_cui(" UMLS:C1234567 ") == "C1234567"

    assert mod.norm_hp("HP:0003731") == "HP_0003731"
    assert mod.norm_hp("http://purl.obolibrary.org/obo/HP_0003731") == "HP_0003731"

    assert mod.to_int_or_none("12") == 12
    assert mod.to_int_or_none("") is None

def test_load_core_and_write_some_outputs(tmp_path):
    mod = load_module(Path("scripts/01_process_makaao_core_to_tables.py"), "process_core")

    inp = Path("data/makaao_sample.csv")
    df = mod.load_core(inp)

    # must contain an integer-like 'index' column created by load_core
    assert "index" in df.columns
    assert pd.api.types.is_integer_dtype(df["index"])

    out_dir = tmp_path / "processed_tables"
    out_dir.mkdir()

    # write a small subset of derived tables; these functions accept explicit out_dir
    mod.write_index_name_en(df, out_dir)
    mod.write_index_hpo_id(df, out_dir)
    mod.write_index_parent_index(df, out_dir)

    for fn in ["index_name_en.csv", "index_hpo_id.csv", "index_parent_index.csv"]:
        p = out_dir / fn
        assert p.exists()
        # header + at least one row
        assert sum(1 for _ in p.open("r", encoding="utf-8")) >= 2
