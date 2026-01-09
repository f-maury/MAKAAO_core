import pandas as pd
import csv
import re
from pathlib import Path

# ---------- CONFIG ----------
INPUT_XLSX = "../data/makaao_core.xlsx"   # path to the Excel file
SHEET = 0                                    # sheet index or name
OUTPUT_CSV = "../data/makaao_core.csv"                           # None → same base name with .csv
MULTI_VALUE_INT_COLS = ["parent_id"]         # columns with multi-value integers
# ----------------------------

pipe_space_rx = re.compile(r"\s*\|\s*")
num_rx = re.compile(r"^[+-]?\d+(?:\.\d+)?$")

def normalize_cell(val):
    if isinstance(val, str):
        s = val.replace("\r\n", "\n").replace("\r", "\n") # e
        s = s.replace("\n", "|")            # replace \n with |
        s = pipe_space_rx.sub("|", s)       # remove spaces around pipes
        return s
    return val

def as_int_str_if_numeric(tok: str) -> str:
    t = tok.strip()
    if t == "":
        return ""                           # preserve blank token from blank line
    if num_rx.match(t):
        try:
            return str(int(float(t)))       # "49.0" → "49"
        except Exception:
            return t
    return t

def fix_multi_int_cell(v): # normalize in case we have several number separated by "|"
    if pd.isna(v):
        return v
    if isinstance(v, (int, float)) and not (isinstance(v, float) and (v != v)):
        try:
            return str(int(v))
        except Exception:
            return v
    if isinstance(v, str):
        parts = v.split("|")
        parts = [as_int_str_if_numeric(p) for p in parts]
        return "|".join(parts)
    return v

def coerce_nullable_int(df: pd.DataFrame, col: str): # force column content to be int
    if col in df.columns:
        ser = pd.to_numeric(df[col], errors="coerce")
        df[col] = ser.astype("Int64")

def main():
    xlsx_path = Path(INPUT_XLSX)
    assert xlsx_path.exists(), f"File not found: {xlsx_path}"

    df = pd.read_excel(xlsx_path, sheet_name=SHEET, engine="openpyxl") # open excel file as pandas df

    # normalize only string/object columns
    obj_cols = [c for c in df.columns if df[c].dtype == "object"] # normalize items in each column
    if obj_cols:
        df[obj_cols] = df[obj_cols].applymap(normalize_cell)

    # fix multi-value integer columns
    for col in MULTI_VALUE_INT_COLS: # normalize when there are several number separated by "|"
        if col in df.columns:
            df[col] = df[col].map(fix_multi_int_cell)

    # `aab_id` is a single integer id
    coerce_nullable_int(df, "aab_id") #turns aab id into int

    out_path = Path(OUTPUT_CSV) if OUTPUT_CSV else xlsx_path.with_suffix(".csv") # write data as csv file
    df.to_csv(out_path, index=False, quoting=csv.QUOTE_MINIMAL, lineterminator="\n", na_rep="")
    print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
