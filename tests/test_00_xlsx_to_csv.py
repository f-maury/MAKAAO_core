import sys
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path("scripts/00_xlsx_to_csv.py")
pytestmark = pytest.mark.skipif(not SCRIPT_PATH.is_file(), reason="scripts/00_xlsx_to_csv.py is not present")



def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_normalize_cell_and_integer_coercions():
    mod = load_module(SCRIPT_PATH, "xlsx_to_csv")

    assert mod.normalize_cell("  A\nB  ") == "  A|B  "
    assert mod.as_int_str_if_numeric("001") == "1"
    assert mod.as_int_str_if_numeric("12.0") == "12"
    assert mod.as_int_str_if_numeric("12.5") == "12"
    assert mod.fix_multi_int_cell(" 01 ; 002;3 ") == "01 ; 002;3"

    frame = pd.DataFrame({"any_col": ["", "  ", None, "001"]})
    result = mod.coerce_nullable_int(frame, "any_col")
    checked = frame if result is None else result

    assert pd.isna(checked.loc[0, "any_col"])
    assert pd.isna(checked.loc[1, "any_col"])
    assert pd.isna(checked.loc[2, "any_col"])
