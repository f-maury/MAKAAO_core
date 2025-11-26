# =========================== IMPORTS ===========================
from __future__ import annotations
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple
import pandas as pd

# =========================== DEFINE PATHS ===========================
INP = Path("../data/makaao_core.csv")  # makaao core table
OUT_DIR = Path("../data/processed_tables/") # where we will store processed tables
OUT_DIR.mkdir(parents=True, exist_ok=True)

# =========================== DEFINE REGEX ===========================
pmid_rx = re.compile(r"^\s*PMID\s*:\s*(\d+)\s*$", re.IGNORECASE) # detect PMID
http_rx = re.compile(r"^https?://", re.IGNORECASE) # detect URL

# =========================== DEFINE FUNCTIONS ===========================
def split_tokens(s: object, primary: str = "|") -> List[str]:
    """
    Split a cell into tokens:
      - First by 'primary' (default '|'), then within each chunk by ';' or newlines.
      - Trim whitespace and drop empty tokens.
    """
    if pd.isna(s):
        return []
    text = str(s).replace("\r\n", "\n").replace("\r", "\n")
    parts: List[str] = []
    for chunk in text.split(primary):
        for sub in chunk.replace("\n", ";").split(";"):
            sub = sub.strip()
            if sub:
                parts.append(sub)
    return parts


def unique_preserve_order(items: Iterable[str]) -> List[str]:
    """Deduplicate while preserving first-seen order (skips empty strings)."""
    seen, out = set(), []
    for x in items:
        if not x:
            continue
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def to_pmid_urls(tokens: Sequence[str]) -> List[str]:
    """
    Normalize tokens into PubMed/HTTP URLs:
      - 'PMID: 12345' -> 'https://pubmed.ncbi.nlm.nih.gov/12345'
      - Keep tokens that already look like HTTP/HTTPS URLs.
      - De-duplicate, order-preserving.
    """
    out, seen = [], set()
    for tok in tokens:
        m = pmid_rx.match(tok)
        url = f"https://pubmed.ncbi.nlm.nih.gov/{m.group(1)}" if m else None
        if not url and http_rx.match(tok):
            url = tok
        if url and url not in seen:
            seen.add(url)
            out.append(url)
    return out


def clean_cui(tok: str) -> str:
    """Return uppercased CUI without the 'CUI:' prefix."""
    s = str(tok).strip()
    if s.lower().startswith("cui:"):
        s = s.split(":", 1)[1]
    return s.upper()


def norm_hp(hp: str) -> str:
    """Normalize 'HP:123' -> 'hp:0000123'; otherwise return as is.""" # add 0s to have 7 digits, because sometimes, when treated as int, the 0s are removed
    s = str(hp).strip()
    if s.upper().startswith("HP:"):
        return "hp:" + s.split(":", 1)[1].zfill(7)
    return s


def first_nonnull(series: Iterable[object]):
    """Return the first non-empty (non-NaN, non-blank) value; else None."""
    for x in series:
        if pd.notna(x) and str(x).strip() != "":
            return x
    return None


def join_col(grp: pd.DataFrame, col: str) -> str:
    """Concatenate non-null values in column 'col' with '|' (or '' if column missing), so that we get a single string"""
    return "|".join(grp[col].dropna().astype(str)) if col in grp else ""


def to_int_or_none(x) -> int | None: # get the correct int from a string
    """Safely parse int-like values such as '12' or '12.0'; return None if not parseable."""
    try:
        return int(float(str(x).strip()))
    except Exception:
        return None


def normalize_loinc_part(token: str) -> str: # get LOINC part identifier, from full URI
 
    if token is None:
        raise ValueError("Empty LOINC part token")

    tok = token.strip()
    if not tok:
        raise ValueError("Empty LOINC part token")

    m_url = re.match(r"^https?://(?:www\.)?loinc\.org/(?:part/)?(LP\d+-\d+)(?:/)?$", tok, flags=re.IGNORECASE)
    if m_url:
        return m_url.group(1).upper()

    tok_no_pref = tok.replace("loinc:", "", 1).replace("LOINC:", "", 1).strip()
    m_code = re.match(r"^(LP\d+-\d+)$", tok_no_pref, flags=re.IGNORECASE)
    if m_code:
        return m_code.group(1).upper()
    raise ValueError(f"Not a valid LOINC part URL or code: {token!r}")



# -------- Slot-alignment helpers --------
def _split_slots(cell: object) -> List[str]:
    """
    Split a cell into '|' slots and trim whitespace.
    IMPORTANT: preserves empty slots ('') for position alignment.
    """
    if pd.isna(cell):
        return []
    text = str(cell).replace("\r\n", "\n").replace("\r", "\n")
    return [slot.strip() for slot in text.split("|")] # return a list


def _items_in_slot(slot: str) -> List[str]:
    """
    Split a single slot into items by ';' or newline, trim, and drop empties.
    Empty slot => [].
    """
    if slot is None:
        return []
    s = slot.replace("\r\n", "\n").replace("\r", "\n").replace("\n", ";").strip()
    if s == "":
        return []
    return [t.strip() for t in s.split(";") if t.strip()] # return a list


def pair_values_and_sources_by_slot(
    values_cell: object,
    sources_cell: object,
    idx_for_error: object = None,
    value_col_name: str = "values",
    source_col_name: str = "sources",
) -> List[Tuple[str, str]]:
    """
    Align values and sources by '|' slots. Within each slot, split on ';'.
    - The number of '|' slots MUST match between values and sources.
    - EXCEPTIONS (padding allowed):
        * ('syn_en','syn_en_source')
        * ('syn_fr','syn_fr_source')
      For these pairs:
        - If sources have fewer slots, pad source slots with '' (empty).
        - If sources have MORE slots, pad value slots with '' (ignored) so counts match.
    - A source slot may be empty => values from that slot get source="". If a source slot
      has multiple tokens (e.g., 'PMID:1;PMID:2'), each value in the paired value slot is
      emitted with each normalized source (Cartesian product).
    """
    v_slots = _split_slots(values_cell)
    s_slots = _split_slots(sources_cell)

    PAD_ALLOWED = {
        ("syn_en", "syn_en_source"), # if mismatch in the number of items, and the number of sources for these columns, we pad. other columns: error
        ("syn_fr", "syn_fr_source"),
    }

    if len(v_slots) != len(s_slots):
        if (value_col_name, source_col_name) in PAD_ALLOWED:
            # Pad the shorter side to match the longer; empty value slots emit nothing
            if len(v_slots) < len(s_slots):
                v_slots = v_slots + [""] * (len(s_slots) - len(v_slots))
            elif len(s_slots) < len(v_slots):
                s_slots = s_slots + [""] * (len(v_slots) - len(s_slots))
        else:
            def _preview(slots: List[str], n=3) -> str:
                return " | ".join(slots[:n]) + ("" if len(slots) <= n else " | ...")
            raise ValueError(
                f"Slot count mismatch at index={idx_for_error}: " # error if mismatch and not in allowed columns
                f"'{value_col_name}' has {len(v_slots)} slot(s) "
                f"vs '{source_col_name}' has {len(s_slots)} slot(s). "
                f"First slots — {value_col_name}: [{_preview(v_slots)}] ; "
                f"{source_col_name}: [{_preview(s_slots)}]"
            )

    pairs: List[Tuple[str, str]] = []
    for v_slot, s_slot in zip(v_slots, s_slots):
        values = _items_in_slot(v_slot)
        src_tokens = _items_in_slot(s_slot)
        srcs = to_pmid_urls(src_tokens)
        if not srcs:
            pairs.extend((v, "") for v in values)
        else:
            for v in values:
                for s in srcs:
                    pairs.append((v, s)) # return list of matched item-source pairs
    return pairs


# =========================== LOAD / PREP ===========================
def load_core(inp: Path) -> pd.DataFrame:
    """
    Load core CSV. Drop rows with missing/invalid aab_id and warn.
    Set 'index' == aab_id (int64).
    """
    import warnings
    t = pd.read_csv(inp)

    if "aab_id" not in t.columns:
        raise KeyError("Missing 'aab_id' column.") # load maakao ore rows that are not empty

    s = t["aab_id"]

    # Drop blank/missing aab_id
    missing = s.isna() | (s.astype(str).str.strip() == "")
    if missing.any():
        idxs = t.index[missing].tolist()[:20]
        warnings.warn(f"load_core: dropping {missing.sum()} row(s) with blank aab_id; sample indices: {idxs}")
        t = t.loc[~missing].copy()
        s = t["aab_id"]

    # Coerce to numeric and drop non-numeric
    num = pd.to_numeric(s, errors="coerce") # if problem with an aab_id, display it
    badnum = num.isna()
    if badnum.any():
        examples = t.loc[badnum, "aab_id"].astype(str).head(20).tolist()
        warnings.warn(f"load_core: dropping {badnum.sum()} row(s) with non-numeric aab_id; examples: {examples}")
        t = t.loc[~badnum].copy()
        num = num.loc[~badnum]

    # Require integers (allow '12.0' but reject '12.3')
    nonint = (num % 1 != 0)
    if nonint.any():
        examples = t.loc[nonint, "aab_id"].astype(str).head(20).tolist() # if other type of problem with an aab_id, also display it
        warnings.warn(f"load_core: dropping {nonint.sum()} row(s) with non-integer aab_id; examples: {examples}")
        t = t.loc[~nonint].copy()
        num = num.loc[~nonint]

    t["index"] = num.astype("int64").values # return filtered data
    return t




# =========================== WRITERS ===========================
def write_index_name_en(df: pd.DataFrame, out_dir: Path) -> None:
    """Write [index, name_en] by taking the first non-empty English name per index.""" # write table with aab_id - name_en
    name_en = (
        df.sort_values("index")
        .groupby("index")["name_en"]
        .apply(first_nonnull)
        .reset_index()
        .rename(columns={"name_en": "name_en"})
    ).dropna(subset=["name_en"])
    name_en.to_csv(out_dir / "index_name_en.csv", index=False)


def write_index_hpo_id(df: pd.DataFrame, out_dir: Path) -> None:
    """Write [index, hpo_id] flattening multi-valued cells and normalizing IDs.""" # write table with aab_id - hpo_id pairs
    rows = []
    for idx, grp in df.groupby("index", sort=True):
        toks = [norm_hp(x) for x in split_tokens(join_col(grp, "hpo_id"))]
        for hp in unique_preserve_order(toks):
            rows.append((idx, hp))
    pd.DataFrame(rows, columns=["index", "hpo_id"]).to_csv(
        out_dir / "index_hpo_id.csv", index=False
    )


def write_index_parent_index(df: pd.DataFrame, out_dir: Path) -> None:
    """Write [index, parent_index] after parsing multi-parents; remove self-links and dups (except index 18)."""
    rows = []
    for idx, grp in df.groupby("index", sort=True):
        raw = split_tokens(join_col(grp, "parent_id"))
        seen = set()
        for p in raw:
            pi = to_int_or_none(p)
            # allow self-link only for index 18
            if pi is None or (pi == idx and idx != 18) or pi in seen:
                continue
            seen.add(pi)
            rows.append((idx, pi))
    parents_df = pd.DataFrame(rows, columns=["index", "parent_index"]).astype(
        {"index": "int64", "parent_index": "int64"}
    )
    parents_df.to_csv(out_dir / "index_parent_index.csv", index=False)



def write_index_syn_en(df: pd.DataFrame, out_dir: Path) -> None:
    """
    Write [index, syns_en, syns_en_source] aligning 'syn_en' with 'syn_en_source' by slots.
    Within a slot: values are ';'-split; sources are normalized URLs and can be multiple.
    """
    rows = []
    for idx, grp in df.groupby("index", sort=True):
        pairs = pair_values_and_sources_by_slot(
            join_col(grp, "syn_en"),
            join_col(grp, "syn_en_source"),
            idx_for_error=idx,
            value_col_name="syn_en",
            source_col_name="syn_en_source",
        )
        # de-duplicate (value,source) while preserving order
        seen = set()
        for v, s in pairs:
            key = (v, s)
            if key not in seen:
                seen.add(key)
                rows.append((idx, v, s))
    pd.DataFrame(rows, columns=["index", "syns_en", "syns_en_source"]).to_csv(
        out_dir / "index_syn_source_en.csv", index=False
    )


def write_index_syn_fr(df: pd.DataFrame, out_dir: Path) -> None:
    """
    Write [index, syns_fr, syns_fr_source] aligning 'syn_fr' with 'syn_fr_source' by slots.
    The first synonym per index is implicitly the French 'name' by convention (no separate file).
    """
    rows = []
    for idx, grp in df.groupby("index", sort=True):
        pairs = pair_values_and_sources_by_slot(
            join_col(grp, "syn_fr"),
            join_col(grp, "syn_fr_source"),
            idx_for_error=idx,
            value_col_name="syn_fr",
            source_col_name="syn_fr_source",
        )
        seen = set()
        for v, s in pairs:
            key = (v, s)
            if key not in seen:
                seen.add(key)
                rows.append((idx, v, s))
    pd.DataFrame(rows, columns=["index", "syns_fr", "syns_fr_source"]).to_csv(
        out_dir / "index_syn_source_fr.csv", index=False
    )


def _write_index_with_sources(
    df: pd.DataFrame,
    out_path: Path,
    id_col_in: str,
    src_col_in: str,
    out_cols: Tuple[str, str, str],
    id_normalizer=lambda x: x,
) -> None:
    """
    Generic helper to write (index, id, sources) style tables.
    Align IDs with sources by '|' slots (must match; no round-robin). Within each slot:
      - IDs are ';'-split and normalized with 'id_normalizer'.
      - Sources are ';'-split and normalized to PubMed/HTTP URLs.
      - Emit Cartesian product of values and sources; if no sources in slot, emit (value, "").
    """
    rows: List[Tuple[int, str, str]] = []
    have_cols = {id_col_in, src_col_in}.issubset(df.columns)
    if not have_cols:
        pd.DataFrame(columns=list(out_cols)).to_csv(out_path, index=False)
        return

    for idx, grp in df.groupby("index", sort=True):
        pairs = pair_values_and_sources_by_slot(
            join_col(grp, id_col_in),
            join_col(grp, src_col_in),
            idx_for_error=idx,
            value_col_name=id_col_in,
            source_col_name=src_col_in,
        )
        # normalize IDs and dedup (id, source)
        seen = set()
        for v, s in pairs:
            v_norm = id_normalizer(v)
            key = (v_norm, s)
            if v_norm and key not in seen:
                seen.add(key)
                rows.append((idx, v_norm, s))

    pd.DataFrame(rows, columns=list(out_cols)).to_csv(out_path, index=False)


def write_index_cui_source(df: pd.DataFrame, out_dir: Path) -> None:
    """Write [index, umls_target_cui, umls_pmids] with slot-aligned pairing; CUIs normalized."""
    _write_index_with_sources(
        df=df,
        out_path=out_dir / "index_cui_source.csv",
        id_col_in="umls_id",
        src_col_in="umls_source",
        out_cols=("index", "umls_target_cui", "umls_pmids"),
        id_normalizer=clean_cui,
    )


def write_index_disease_source(df: pd.DataFrame, out_dir: Path) -> None:
    """Write [index, related_diseases_id, diseases_pmids] with slot-aligned pairing."""
    _write_index_with_sources(
        df=df,
        out_path=out_dir / "index_disease_source.csv",
        id_col_in="disease_id",
        src_col_in="disease_source",
        out_cols=("index", "related_diseases_id", "diseases_pmids"),
    )


def write_index_uniprot_source(df: pd.DataFrame, out_dir: Path) -> None:
    """Write [index, uniprot_target_id, uniprot_pmids] with slot-aligned pairing."""
    _write_index_with_sources(
        df=df,
        out_path=out_dir / "index_uniprot_source.csv",
        id_col_in="uniprot_id",
        src_col_in="uniprot_source",
        out_cols=("index", "uniprot_target_id", "uniprot_pmids"),
    )


def write_index_chebi_source(df: pd.DataFrame, out_dir: Path) -> None:
    """Write [index, chebi_target_id, chebi_pmids] with slot-aligned pairing."""
    _write_index_with_sources(
        df=df,
        out_path=out_dir / "index_chebi_source.csv",
        id_col_in="chebi_id",
        src_col_in="chebi_source",
        out_cols=("index", "chebi_target_id", "chebi_pmids"),
    )


def write_index_loinc(df: pd.DataFrame, out_dir: Path) -> None:
    """
    Write [aab_id, loinc_id] by splitting 'loinc_part_id' into tokens and normalizing.
    Deduplicate per index.
    """
    rows = []
    if "loinc_part_id" in df.columns:
        for idx, grp in df.groupby("index", sort=True):
            tokens = split_tokens(join_col(grp, "loinc_part_id"))
            seen = set()
            for tok in tokens:
                uri = normalize_loinc_part(tok)
                if uri and uri not in seen:
                    seen.add(uri)
                    rows.append((idx, uri))
    if rows:
        pd.DataFrame(rows, columns=["aab_id", "loinc_id"]).to_csv(
            out_dir / "index_loinc.csv", index=False
        )
    else:
        pd.DataFrame(columns=["aab_id", "loinc_id"]).to_csv(
            out_dir / "index_loinc.csv", index=False
        )


def print_report(out_dir: Path) -> None:
    """Print a short report of the number of data rows per generated CSV."""
    files = [
        "index_name_en.csv",
        "index_hpo_id.csv",
        "index_parent_index.csv",
        "index_syn_source_en.csv",
        "index_syn_source_fr.csv",
        "index_cui_source.csv",
        "index_disease_source.csv",
        "index_uniprot_source.csv",
        "index_chebi_source.csv",
        "index_loinc.csv",
    ]
    for fn in files:
        p = out_dir / fn
        try:
            rows = max(0, sum(1 for _ in open(p, "rb")) - 1)
        except Exception:
            rows = "n/a"
        print(f"{p}: rows={rows}")


# =========================== MAIN ===========================
def main():
    """
    Orchestrate the normalization pipeline:
      1) Load core CSV and derive integer 'index'.
      2) Write per-entity/per-relation CSVs into OUT_DIR (slot-aligned sources; no round-robin).
      3) Print a summary report.
    """
    df = load_core(INP)

    write_index_name_en(df, OUT_DIR)
    write_index_hpo_id(df, OUT_DIR)
    write_index_parent_index(df, OUT_DIR)

    write_index_syn_en(df, OUT_DIR)
    write_index_syn_fr(df, OUT_DIR)

    write_index_cui_source(df, OUT_DIR)
    write_index_disease_source(df, OUT_DIR)
    write_index_uniprot_source(df, OUT_DIR)
    write_index_chebi_source(df, OUT_DIR)

    write_index_loinc(df, OUT_DIR)

    print_report(OUT_DIR)


if __name__ == "__main__":
    main()
