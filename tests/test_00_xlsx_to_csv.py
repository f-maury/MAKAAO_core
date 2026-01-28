import importlib.util
from pathlib import Path

import pandas as pd


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

    assert mod.fix_multi_int_cell(" 01 ; 002;3 ") == "01 ; 002;3"

    # coerce_nullable_int expects a DataFrame + column name
    df = pd.DataFrame({"any_col": ["", "  ", None, "001"]})
    out = mod.coerce_nullable_int(df, "any_col")

    # function may modify in-place or return df
    df2 = df if out is None else out

    # minimal invariant: blanks/None become NA-like
    assert pd.isna(df2.loc[0, "any_col"])
    assert pd.isna(df2.loc[1, "any_col"])
    assert pd.isna(df2.loc[2, "any_col"])
