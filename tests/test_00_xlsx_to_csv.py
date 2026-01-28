import importlib.util
from pathlib import Path

def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod

def test_normalize_cell_and_int_coercions():
    mod = load_module(Path("scripts/00_xlsx_to_csv.py"), "xlsx_to_csv")

    assert mod.normalize_cell("  A\nB  ") == "  A|B  "
    assert mod.as_int_str_if_numeric("001") == "1"
    assert mod.as_int_str_if_numeric("12.0") == "12"
    assert mod.as_int_str_if_numeric("12.5") == "12"

    # multi-int cells like "1; 2;3" should be normalized to "1;2;3" (implementation-dependent)
    assert mod.fix_multi_int_cell(" 01 ; 002;3 ") == "01 ; 002;3"

    assert mod.coerce_nullable_int("", "any_col") is None
    assert mod.coerce_nullable_int("  ", "any_col") is None
    assert mod.coerce_nullable_int("12", "any_col") == 12
