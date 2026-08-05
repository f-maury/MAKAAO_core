#!/usr/bin/env python3
import sys
import csv
import pandas as pd
import xml.etree.ElementTree as ET
import json
import os
import re
import time
import requests
from pathlib import Path
from typing import Optional, Tuple, List
from urllib.parse import quote
from collections import Counter, defaultdict

from concurrent.futures import ThreadPoolExecutor, as_completed

# --- INPUT PATHS---
# Paths are anchored to the script location, not to the current working directory.
# Expected project layout:
#   project/scripts/02_create_enrichment_tables.py
#   project/data/makaao_core.csv
#   project/data/enrichment_tables/
#   project/kg/
SCRIPT_DIR = Path(__file__).resolve().parent


def _find_project_dir(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "data").is_dir():
            return candidate
    raise RuntimeError(
        f"Could not locate project root above {start}; expected a project data/ directory."
    )


PROJECT_DIR = _find_project_dir(SCRIPT_DIR)
DATA_DIR = str(PROJECT_DIR / "data")

# NOTE: Update IN_PATH to point to your local MRCONSO.RRF file.
IN_PATH = "/mnt/d/umls-2024AB-full_metamor/2024AB-full/2024AB/2024AB/META/MRCONSO.RRF"
XML_PATH = os.path.join(DATA_DIR, "en_product4.xml")
ENRICH_DIR = os.path.join(DATA_DIR, "enrichment_tables")
INPUT_CSV_CORE = os.path.join(DATA_DIR, "makaao_core.csv")
LOINC_PART_LINK_CSV = os.path.join(DATA_DIR, "LoincPartLink_Primary.csv")

# Any English, unsuppressed MRCONSO atom may provide a label for a UMLS CUI.
# ORPHANET atoms are also used to derive ORPHA↔UMLS identifier mappings.
UMLS_MAPPING_SABS = {"ORPHANET"}


# UMLS  we prioritize some type of UMLS terms
TTY_PRIORITY = {
    "PN": 0,   # Preferred name
    "PT": 1,   # Preferred term
    "MH": 2,   # MeSH heading
    "FN": 3,   # Full form
    "SY": 4,   # Synonym
}
DEFAULT_TTY_RANK = 100

# Final output files
OUT_ORPHA_LINKS = os.path.join(ENRICH_DIR, "orphanet_hpo_links.csv")
OUT_ORPHA_UMLS_MAPPINGS = os.path.join(ENRICH_DIR, "orpha_umls_mappings.json")
OUTPUT_CSV_FINAL = os.path.join(ENRICH_DIR, "code_names.csv")
OUT_LOINC_PART_TESTS = os.path.join(ENRICH_DIR, "loinc_part_test_dict.json")
OUT_LOINC_LABELS = os.path.join(ENRICH_DIR, "loinc_labels.json")
REPORT_PATH = os.path.join(ENRICH_DIR, "enrichment_report.md")

# API URLs
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{acc}"
OLS4_TERM_API = "https://www.ebi.ac.uk/ols4/api/ontologies/{onto}/terms"
CHEBI_COMPOUND_API = (
    "https://www.ebi.ac.uk/chebi/backend/api/public/compound/{chebi_id}/"
)
HEADERS = {"Accept": "application/json"}

# Global storage for in-memory mapping and tracking
umls_names_map = {}
# Metadata for the chosen UMLS name per CUI (e.g., SAB, preferred vs synonym)
umls_name_meta = {}
# Direct HPO-ID -> HPO label extracted from MRCONSO (SAB=HPO)
hpo_names_map = {}
# Local labels from MRCONSO vocabulary atoms.
orpha_names_map = {}
orpha_xml_names_map = {}
chebi_names_map = {}
loinc_part_names_map = {}
loinc_term_names_map = {}
results_tracker = {
    "UniProt": {"success": 0, "fail": 0},
    "ORPHA": {"success": 0, "fail": 0},
    "ChEBI": {"success": 0, "fail": 0},
    "UMLS": {"success": 0, "fail": 0},
    "HPO": {"success": 0, "fail": 0},
    "LOINC": {"success": 0, "fail": 0}
}

# Increase CSV field size limit (MRCONSO.RRF is a big file with long fields)
_lim = sys.maxsize
while True:
    try:
        csv.field_size_limit(_lim)
        break
    except OverflowError:
        _lim //= 10


# --- PART 1: PROCESS MRCONSO.RRF (Direct to Memory) ---
def process_mrconso(
    input_cuis: set[str],
    input_orphas: set[str],
    input_chebi_ids: set[str],
    required_hpo_ids: set[str],
    required_loinc_parts: set[str],
    required_loinc_terms: set[str],
) -> None:
    """Scan MRCONSO once for relevant labels and ORPHA↔UMLS mappings."""
    print(f"Processing MRCONSO.RRF from {IN_PATH}...")
    if not os.path.isfile(IN_PATH):
        raise FileNotFoundError(f"Required MRCONSO.RRF not found: {IN_PATH}")

    os.makedirs(ENRICH_DIR, exist_ok=True)
    orpha_to_umls = defaultdict(set)
    umls_to_orpha = defaultdict(set)

    # Candidate ordering, from strongest to weakest evidence:
    #   1) ISPREF=Y
    #   2) preferred TTY rank
    #   3) TS=P
    #   4) shorter label
    #   5) stable lexical/source fields for reproducibility
    best: dict[str, tuple] = {}
    best_hpo: dict[str, tuple] = {}
    best_orpha: dict[str, tuple] = {}
    best_chebi: dict[str, tuple] = {}
    best_loinc_parts: dict[str, tuple] = {}
    best_loinc_terms: dict[str, tuple] = {}

    total_rows = 0
    relevant_label_rows = 0
    mapping_rows = 0
    relevant_loinc_rows = 0

    with open(IN_PATH, "r", encoding="utf-8", errors="replace", newline="") as fin:
        reader = csv.reader(fin, delimiter="|")
        for parts in reader:
            total_rows += 1
            if len(parts) < 18:
                continue

            cui_raw = parts[0]
            lat = parts[1]
            ts = parts[2]
            ispref = parts[6]
            scui = parts[9]
            sdui = parts[10]
            sab = parts[11]
            tty = parts[12]
            code = parts[13]
            label = re.sub(r"\s+", " ", parts[14] or "").strip()
            suppress = parts[16]

            cui = norm_umls(cui_raw)

            # Mapping extraction is independent of language and label selection.
            if suppress == "N" and sab.upper() in UMLS_MAPPING_SABS:
                orpha = orpha_code_from_mrconso_atom(code, scui, sdui)
                if cui and orpha and (cui in input_cuis or orpha in input_orphas):
                    orpha_to_umls[orpha].add(cui)
                    umls_to_orpha[cui].add(orpha)
                    mapping_rows += 1

            if lat != "ENG" or suppress != "N" or not label:
                continue

            if cui and cui in input_cuis:
                relevant_label_rows += 1
                cand = (
                    0 if ispref == "Y" else 1,
                    TTY_PRIORITY.get(tty, DEFAULT_TTY_RANK),
                    0 if ts == "P" else 1,
                    len(label),
                    label.casefold(),
                    label,
                    sab,
                    tty,
                    code,
                )
                current = best.get(cui)
                if current is None or cand < current:
                    best[cui] = cand

            if sab.upper() == "ORPHANET":
                orpha_id = orpha_code_from_mrconso_atom(code, scui, sdui)
                if orpha_id and orpha_id in input_orphas:
                    local_cand = (
                        0 if ispref == "Y" else 1,
                        TTY_PRIORITY.get(tty, DEFAULT_TTY_RANK),
                        0 if ts == "P" else 1,
                        len(label),
                        label.casefold(),
                        label,
                        code,
                    )
                    current = best_orpha.get(orpha_id)
                    if current is None or local_cand < current:
                        best_orpha[orpha_id] = local_cand

            if sab.upper() == "CHEBI":
                chebi_id = norm_chebi(code) or norm_chebi(scui) or norm_chebi(sdui)
                if chebi_id and chebi_id in input_chebi_ids:
                    local_cand = (
                        0 if ispref == "Y" else 1,
                        TTY_PRIORITY.get(tty, DEFAULT_TTY_RANK),
                        0 if ts == "P" else 1,
                        len(label),
                        label.casefold(),
                        label,
                        code,
                    )
                    current = best_chebi.get(chebi_id)
                    if current is None or local_cand < current:
                        best_chebi[chebi_id] = local_cand

            if sab == "HPO":
                hpo_id = norm_hpo(code)
                if hpo_id and hpo_id in required_hpo_ids:
                    hpo_cand = (
                        0 if ispref == "Y" else 1,
                        TTY_PRIORITY.get(tty, DEFAULT_TTY_RANK),
                        0 if ts == "P" else 1,
                        len(label),
                        label.casefold(),
                        label,
                        code,
                    )
                    current = best_hpo.get(hpo_id)
                    if current is None or hpo_cand < current:
                        best_hpo[hpo_id] = hpo_cand

            if sab.upper() == "LNC":
                loinc_code = loinc_code_from_mrconso_atom(code, scui, sdui)
                if loinc_code in required_loinc_parts or loinc_code in required_loinc_terms:
                    relevant_loinc_rows += 1
                    loinc_cand = (
                        0 if ispref == "Y" else 1,
                        TTY_PRIORITY.get(tty, DEFAULT_TTY_RANK),
                        0 if ts == "P" else 1,
                        len(label),
                        label.casefold(),
                        label,
                        code,
                    )
                    target = (
                        best_loinc_parts
                        if loinc_code in required_loinc_parts
                        else best_loinc_terms
                    )
                    current = target.get(loinc_code)
                    if current is None or loinc_cand < current:
                        target[loinc_code] = loinc_cand

    for cui, cand in best.items():
        label, sab, tty, code = cand[5], cand[6], cand[7], cand[8]
        umls_names_map[cui] = label
        umls_name_meta[cui] = {
            "sab": sab,
            "tty": tty,
            "code": code,
            "ispref": cand[0] == 0,
            "ts_preferred": cand[2] == 0,
        }

    for hpo_id, cand in best_hpo.items():
        hpo_names_map[hpo_id] = cand[5]
    for orpha_id, cand in best_orpha.items():
        orpha_names_map[orpha_id] = cand[5]
    for chebi_id, cand in best_chebi.items():
        chebi_names_map[chebi_id] = cand[5]
    for loinc_code, cand in best_loinc_parts.items():
        loinc_part_names_map[loinc_code] = cand[5]
    for loinc_code, cand in best_loinc_terms.items():
        loinc_term_names_map[loinc_code] = cand[5]

    write_orpha_umls_mappings(
        OUT_ORPHA_UMLS_MAPPINGS,
        input_cuis,
        input_orphas,
        orpha_to_umls,
        umls_to_orpha,
    )
    mapping_pair_count = sum(len(values) for values in orpha_to_umls.values())
    mapped_input_orphas = set(input_orphas) & set(orpha_to_umls)
    mapped_input_cuis = set(input_cuis) & set(umls_to_orpha)
    unmapped_input_orphas = set(input_orphas) - set(orpha_to_umls)
    additional_orphas = set(orpha_to_umls) - set(input_orphas)
    print(
        f"Part 1: Wrote {mapping_pair_count} ORPHA/UMLS pairs to "
        f"{OUT_ORPHA_UMLS_MAPPINGS} "
        f"({len(mapped_input_orphas)}/{len(input_orphas)} input ORPHA codes mapped; "
        f"{len(mapped_input_cuis)}/{len(input_cuis)} input CUIs mapped; "
        f"{len(additional_orphas)} additional ORPHA counterparts discovered from input CUIs; "
        f"{len(unmapped_input_orphas)} input ORPHA codes unmapped)."
    )
    print(
        "Part 1: Loaded "
        f"{len(umls_names_map)}/{len(input_cuis)} requested UMLS labels and "
        f"{len(hpo_names_map)}/{len(required_hpo_ids)} requested HPO labels "
        f"(relevant UMLS label rows: {relevant_label_rows}; "
        f"ORPHANET mapping rows: {mapping_rows}; "
        f"relevant LOINC label rows: {relevant_loinc_rows}; scanned: {total_rows})."
    )
    print(
        "Part 1: Loaded "
        f"{len(loinc_part_names_map)}/{len(required_loinc_parts)} requested LOINC Part labels and "
        f"{len(loinc_term_names_map)}/{len(required_loinc_terms)} requested LOINC Term labels."
    )
    print(
        "Part 1: Loaded local labels for "
        f"ORPHA MRCONSO {len(orpha_names_map)}/{len(input_orphas)}; "
        f"Orphanet XML {len(orpha_xml_names_map)}/{len(input_orphas)}; "
        f"ChEBI {len(chebi_names_map)}/{len(input_chebi_ids)}."
    )
    print(
        "Part 1: UMLS label selection accepts every MRCONSO source "
        "(English, unsuppressed atoms only)."
    )


# --- PART 2: PROCESS ORPHANET XML ---
def local_tag(tag: str) -> str:
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag

def get_lang(elem):
    return elem.attrib.get("lang") or elem.attrib.get("{http://www.w3.org/XML/1998/namespace}lang")

def freq_rank(txt: str) -> Tuple[int, str]:
    if not txt: return (-1, "Unknown")
    t = txt.strip().lower()
    if "obligate" in t or "always present" in t or "100%" in t: return (5, txt.strip())
    if "very frequent" in t or ("99" in t and "80" in t): return (4, txt.strip())
    if "frequent" in t and "very" not in t: return (3, txt.strip())
    return (-1, txt.strip())

def process_orphanet(input_orphas: set[str]) -> set[str]:
    """Write phenotype links only for Orphanet codes used by MAKAAO."""
    print(f"Processing Orphanet XML from {XML_PATH}...")
    if not os.path.isfile(XML_PATH):
        raise FileNotFoundError(f"Required Orphanet XML not found: {XML_PATH}")

    rows = []
    required_hpo_ids: set[str] = set()
    try:
        root = ET.parse(XML_PATH).getroot()
        for disorder in root.iter():
            if local_tag(disorder.tag) != "Disorder":
                continue
            oc, assoc_list, disorder_name = None, None, None
            for child in disorder:
                name = local_tag(child.tag)
                if name == "OrphaCode" and (child.text or "").strip():
                    oc = str(int(child.text.strip()))
                elif name == "Name" and (child.text or "").strip():
                    language = (get_lang(child) or "").lower()
                    if disorder_name is None or language == "en":
                        disorder_name = re.sub(r"\s+", " ", child.text).strip()
                elif name == "HPODisorderAssociationList":
                    assoc_list = child
            if not oc or oc not in input_orphas:
                continue
            if disorder_name:
                orpha_xml_names_map[oc] = disorder_name
            if assoc_list is None:
                continue

            for assoc in assoc_list:
                hpo_id, hpo_term = None, None
                hpo = next((node for node in assoc if local_tag(node.tag) == "HPO"), None)
                if hpo is not None:
                    hid = next((node for node in hpo if local_tag(node.tag) == "HPOId"), None)
                    htm = next((node for node in hpo if local_tag(node.tag) == "HPOTerm"), None)
                    if hid is not None:
                        hpo_id = norm_hpo((hid.text or "").strip())
                    if htm is not None:
                        hpo_term = (htm.text or "").strip()
                if not hpo_id:
                    continue

                frequency_name = None
                frequency = next(
                    (node for node in assoc if local_tag(node.tag) == "HPOFrequency"),
                    None,
                )
                if frequency is not None:
                    names = [
                        node for node in frequency
                        if local_tag(node.tag) == "Name" and (node.text or "").strip()
                    ]
                    english = next(
                        (node for node in names if (get_lang(node) or "").lower() == "en"),
                        None,
                    )
                    chosen = english or (names[0] if names else None)
                    if chosen is not None:
                        frequency_name = chosen.text.strip()

                rank, frequency_name = freq_rank(frequency_name)
                if rank >= 3:
                    required_hpo_ids.add(hpo_id)
                    rows.append(
                        (
                            f"http://www.orpha.net/ORDO/Orphanet_{oc}",
                            f"http://purl.obolibrary.org/obo/{hpo_id.replace(':', '_')}",
                            hpo_term or "",
                            frequency_name,
                            rank,
                        )
                    )

        df = pd.DataFrame(
            rows,
            columns=["orpha_code", "HPOId", "HPOTerm", "frequency", "rank"],
        )
        if not df.empty:
            df = (
                df.sort_values(
                    ["orpha_code", "rank", "HPOId", "HPOTerm"],
                    ascending=[True, False, True, True],
                    kind="mergesort",
                )
                .drop_duplicates(["orpha_code", "HPOId"], keep="first")
                .drop(columns=["rank"])
            )

        os.makedirs(ENRICH_DIR, exist_ok=True)
        temporary = OUT_ORPHA_LINKS + ".tmp"
        df.to_csv(temporary, index=False, lineterminator="\n")
        os.replace(temporary, OUT_ORPHA_LINKS)
        print(
            f"Part 2: Wrote {len(df)} rows for {len(input_orphas)} input "
            f"Orphanet codes to {OUT_ORPHA_LINKS}"
        )
        return required_hpo_ids
    except Exception as exc:
        temporary = OUT_ORPHA_LINKS + ".tmp"
        if os.path.exists(temporary):
            os.remove(temporary)
        raise RuntimeError(f"Failed to process Orphanet XML {XML_PATH}") from exc


# --- PART 3: ENRICHMENT USING APIS FOR WHA WAS NOT IN MRCONSO -----
def req_get(url: str, params: Optional[dict] = None) -> Optional[requests.Response]:
    for i in range(3):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=20)
            if r.status_code in (200, 404): return r
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(0.5 * (2**i))
                continue
            return r
        except requests.RequestException:
            time.sleep(0.5 * (2**i))
    return None

def split_items(cell: str) -> List[str]:
    if not cell: return []
    return [s for s in re.split(r"[ \t\r\n|,;]+", cell.strip()) if s]

def split_items_pipe(cell: str) -> List[str]:
    if not cell: return []
    return [s.strip() for s in re.split(r"[|\r\n]+", cell.strip()) if s.strip()]

# Normalizers
def norm_uniprot(x: str) -> Optional[str]:
    x = re.sub(r"(?i)^UP:?", "", (x or "").strip())
    return x.upper() if re.fullmatch(r"[A-Za-z0-9]{6,10}", x) else None

def norm_umls(x: str) -> Optional[str]:
    m = re.fullmatch(r"(?i)(?:CUI:)?(C\d{7,8})", (x or "").strip())
    return m.group(1).upper() if m else None

def norm_orpha(x: str) -> Optional[str]:
    s = (x or "").strip()
    if not s:
        return None
    if s.lower().startswith(("http://", "https://")):
        tail = s.split("?", 1)[0].split("#", 1)[0].rstrip("/").rsplit("/", 1)[-1]
        m = re.fullmatch(r"(?i)(?:ORPHA|ORPHANET)[_:\-]?(\d+)", tail)
        if m:
            return str(int(m.group(1)))
        if "orpha" in s.lower() and re.fullmatch(r"\d+", tail):
            return str(int(tail))
        return None
    m = re.fullmatch(r"(?i)(?:ORPHA|ORPHANET)[_:\-]?(\d+)", s)
    return str(int(m.group(1))) if m else None


def orpha_code_from_mrconso_atom(*values: str) -> Optional[str]:
    """Extract an Orphanet code from MRCONSO CODE, SCUI, or SDUI fields."""
    for value in values:
        s = (value or "").strip()
        if not s:
            continue
        code = norm_orpha(s)
        if code:
            return code
        tail = s.split("?", 1)[0].split("#", 1)[0].rstrip("/").rsplit("/", 1)[-1]
        code = norm_orpha(tail)
        if code:
            return code
        if re.fullmatch(r"\d+", tail):
            return str(int(tail))
    return None


def collect_core_orpha_umls_identifiers() -> Tuple[set, set]:
    """Collect CUIs and explicit Orphanet codes from the MAKAAO core CSV."""
    cuis, orphas = set(), set()
    if not os.path.isfile(INPUT_CSV_CORE):
        raise FileNotFoundError(f"Required core input not found: {INPUT_CSV_CORE}")

    with open(INPUT_CSV_CORE, newline="", encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        cols = {
            re.sub(r"[\s\ufeff]+", "", c.lower()): c
            for c in (rdr.fieldnames or [])
            if c is not None
        }
        for row in rdr:
            if "umls_id" in cols:
                for token in split_items(row.get(cols["umls_id"], "")):
                    cui = norm_umls(token)
                    if cui:
                        cuis.add(cui)
            if "disease_id" in cols:
                for token in split_items(row.get(cols["disease_id"], "")):
                    cui = norm_umls(token)
                    if cui:
                        cuis.add(cui)
                        continue
                    orpha = norm_orpha(token)
                    if orpha:
                        orphas.add(orpha)

    print(
        "Part 1: Core identifiers for ORPHA/UMLS mapping — "
        f"CUIs:{len(cuis)} ORPHA:{len(orphas)}"
    )
    return cuis, orphas


def collect_core_chebi_identifiers() -> set[str]:
    """Collect ChEBI identifiers used by the MAKAAO core CSV."""
    identifiers: set[str] = set()
    with open(INPUT_CSV_CORE, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = {
            re.sub(r"[\s\ufeff]+", "", column.lower()): column
            for column in (reader.fieldnames or [])
            if column is not None
        }
        column = columns.get("chebi_id")
        if column:
            for row in reader:
                for token in split_items(row.get(column, "")):
                    identifier = norm_chebi(token)
                    if identifier:
                        identifiers.add(identifier)
    return identifiers


def atomic_write_json(output_path: str, payload) -> None:
    """Write JSON deterministically and atomically."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    temporary = output_path + ".tmp"
    try:
        with open(temporary, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        os.replace(temporary, output_path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def write_orpha_umls_mappings(
    output_path: str,
    input_cuis: set,
    input_orphas: set,
    orpha_to_umls: dict,
    umls_to_orpha: dict,
) -> None:
    """Write deterministic, explicitly bidirectional ORPHA/UMLS dictionaries."""
    forward_pairs = {
        (orpha, cui)
        for orpha, cuis in orpha_to_umls.items()
        for cui in cuis
    }
    reverse_pairs = {
        (orpha, cui)
        for cui, orphas in umls_to_orpha.items()
        for orpha in orphas
    }
    if forward_pairs != reverse_pairs:
        raise RuntimeError(
            "Internal ORPHA/UMLS mapping error: forward and reverse pairs differ"
        )

    payload = {
        "schema_version": 1,
        "source": {
            "file": os.path.basename(IN_PATH),
            "umls_source_abbreviation": "ORPHANET",
        },
        "input": {
            "file": os.path.basename(INPUT_CSV_CORE),
            "orpha_codes": sorted(input_orphas, key=lambda value: int(value)),
            "umls_cuis": sorted(input_cuis),
        },
        "orpha_to_umls": {
            orpha: sorted(orpha_to_umls[orpha])
            for orpha in sorted(orpha_to_umls, key=lambda value: int(value))
        },
        "umls_to_orpha": {
            cui: sorted(umls_to_orpha[cui], key=lambda value: int(value))
            for cui in sorted(umls_to_orpha)
        },
    }
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    temporary = output_path + ".tmp"
    try:
        with open(temporary, "w", encoding="utf-8", newline="\n") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
            f.write("\n")
        os.replace(temporary, output_path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def norm_chebi(x: str) -> Optional[str]:
    m = re.fullmatch(r"(?i)(?:CHEBI:|CHE:)?(\d+)", (x or "").strip())
    return m.group(1) if m else None

def norm_loinc_part(x: str) -> Optional[str]:
    value = (x or "").strip().upper()
    value = re.sub(r"(?i)^https?://(?:www\.)?loinc\.org/", "", value)
    value = value.split("?", 1)[0].split("#", 1)[0].strip("/")
    value = re.sub(r"(?i)^LOINC\s*:\s*", "", value).strip()
    return value if re.fullmatch(r"LP\d+(?:-\d+)?", value) else None

def norm_loinc_term(x: str) -> Optional[str]:
    value = (x or "").strip().upper()
    value = re.sub(r"(?i)^https?://(?:www\.)?loinc\.org/", "", value)
    value = value.split("?", 1)[0].split("#", 1)[0].strip("/")
    value = re.sub(r"(?i)^LOINC\s*:\s*", "", value).strip()
    return value if re.fullmatch(r"\d+-\d+", value) else None


def loinc_code_from_mrconso_atom(*values: str) -> Optional[str]:
    """Return a canonical LOINC Part or Term code from MRCONSO identifiers."""
    for value in values:
        if part := norm_loinc_part(value):
            return part
        if term := norm_loinc_term(value):
            return term
    return None


def _clean_label(value: str) -> str:
    """Normalize whitespace without changing the label's capitalization."""
    return re.sub(r"\s+", " ", value or "").strip()


def _select_stable_label(current: str, candidate: str) -> str:
    """Select a deterministic non-empty label when duplicate rows disagree."""
    candidate = _clean_label(candidate)
    if not candidate:
        return current
    current = _clean_label(current)
    if not current:
        return candidate
    return min((current, candidate), key=lambda value: (len(value), value.casefold(), value))


def collect_core_loinc_parts_and_labels() -> tuple[set[str], dict[str, str]]:
    """Return core LOINC Parts and any labels paired with them in the core CSV."""
    if not os.path.isfile(INPUT_CSV_CORE):
        raise FileNotFoundError(f"Required core input not found: {INPUT_CSV_CORE}")
    parts: set[str] = set()
    labels: dict[str, str] = {}
    with open(INPUT_CSV_CORE, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = {
            re.sub(r"[\s\ufeff]+", "", column.lower()): column
            for column in (reader.fieldnames or [])
            if column is not None
        }
        part_column = columns.get("loinc_part_id")
        label_column = columns.get("loinc_part")
        if not part_column:
            return parts, labels
        for row in reader:
            raw_parts = split_items_pipe(row.get(part_column, ""))
            raw_labels = (
                split_items_pipe(row.get(label_column, "")) if label_column else []
            )
            for index, raw_part in enumerate(raw_parts):
                part = norm_loinc_part(raw_part)
                if not part:
                    continue
                parts.add(part)
                candidate = raw_labels[index] if index < len(raw_labels) else ""
                labels[part] = _select_stable_label(labels.get(part, ""), candidate)
    return parts, labels


def parse_loinc_primary_links(
    input_parts: set[str],
) -> tuple[dict[str, list[str]], set[str], dict[str, str], dict[str, str]]:
    """Parse relevant Primary links and their local Part/Term labels."""
    if not os.path.isfile(LOINC_PART_LINK_CSV):
        raise FileNotFoundError(
            f"Required LOINC Part-link table not found: {LOINC_PART_LINK_CSV}"
        )
    links: dict[str, set[str]] = defaultdict(set)
    part_labels: dict[str, str] = {}
    term_labels: dict[str, str] = {}
    required_headers = {
        "LoincNumber",
        "LongCommonName",
        "PartNumber",
        "PartName",
        "PartCodeSystem",
        "LinkTypeName",
    }
    with open(LOINC_PART_LINK_CSV, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        headers = set(reader.fieldnames or [])
        missing = sorted(required_headers - headers)
        if missing:
            raise ValueError(
                f"{LOINC_PART_LINK_CSV} is missing required columns: {', '.join(missing)}"
            )
        for row_number, row in enumerate(reader, start=2):
            if (row.get("LinkTypeName") or "").strip().casefold() != "primary":
                continue
            code_system = (row.get("PartCodeSystem") or "").strip().rstrip("/").casefold()
            if code_system not in {"http://loinc.org", "https://loinc.org"}:
                continue
            part = norm_loinc_part(row.get("PartNumber", ""))
            if not part or part not in input_parts:
                continue
            term = norm_loinc_term(row.get("LoincNumber", ""))
            if not term:
                raise ValueError(
                    f"{LOINC_PART_LINK_CSV}: row {row_number} has an invalid LoincNumber "
                    f"for relevant Part {part}: {row.get('LoincNumber')!r}"
                )
            links[part].add(term)
            part_labels[part] = _select_stable_label(
                part_labels.get(part, ""), row.get("PartName", "")
            )
            term_labels[term] = _select_stable_label(
                term_labels.get(term, ""), row.get("LongCommonName", "")
            )
    part_tests = {
        part: sorted(links.get(part, set()), key=lambda value: (value.casefold(), value))
        for part in sorted(input_parts, key=lambda value: (value.casefold(), value))
    }
    terms = {term for values in part_tests.values() for term in values}
    return part_tests, terms, part_labels, term_labels



def apply_loinc_local_label_fallbacks(
    input_parts: set[str],
    input_terms: set[str],
    link_part_labels: dict[str, str],
    link_term_labels: dict[str, str],
    core_part_labels: dict[str, str],
) -> dict[str, int]:
    """Fill missing MRCONSO LOINC labels from local CSV files only."""
    mrconso_parts = len(set(input_parts) & set(loinc_part_names_map))
    mrconso_terms = len(set(input_terms) & set(loinc_term_names_map))
    part_link_fallbacks = 0
    core_fallbacks = 0
    term_link_fallbacks = 0

    for part in sorted(input_parts):
        if loinc_part_names_map.get(part):
            continue
        if label := _clean_label(link_part_labels.get(part, "")):
            loinc_part_names_map[part] = label
            part_link_fallbacks += 1
        elif label := _clean_label(core_part_labels.get(part, "")):
            loinc_part_names_map[part] = label
            core_fallbacks += 1

    for term in sorted(input_terms):
        if loinc_term_names_map.get(term):
            continue
        if label := _clean_label(link_term_labels.get(term, "")):
            loinc_term_names_map[term] = label
            term_link_fallbacks += 1

    return {
        "total_parts": len(input_parts),
        "total_terms": len(input_terms),
        "mrconso_parts": mrconso_parts,
        "mrconso_terms": mrconso_terms,
        "part_link_fallbacks": part_link_fallbacks,
        "core_part_fallbacks": core_fallbacks,
        "term_link_fallbacks": term_link_fallbacks,
        "resolved_parts": sum(bool(loinc_part_names_map.get(code)) for code in input_parts),
        "resolved_terms": sum(bool(loinc_term_names_map.get(code)) for code in input_terms),
    }

def write_loinc_enrichment_files(
    part_tests: dict[str, list[str]],
    input_parts: set[str],
    input_terms: set[str],
    label_stats: dict[str, int],
) -> tuple[list[str], list[str]]:
    """Write the two deterministic JSON files consumed by script 03."""
    missing_parts = sorted(
        code for code in input_parts if not loinc_part_names_map.get(code)
    )
    missing_terms = sorted(
        code for code in input_terms if not loinc_term_names_map.get(code)
    )
    labels = {
        "parts": {
            code: loinc_part_names_map.get(code, "")
            for code in sorted(input_parts, key=lambda value: (value.casefold(), value))
        },
        "tests": {
            code: loinc_term_names_map.get(code, "")
            for code in sorted(input_terms, key=lambda value: (value.casefold(), value))
        },
    }
    atomic_write_json(OUT_LOINC_PART_TESTS, part_tests)
    atomic_write_json(OUT_LOINC_LABELS, labels)
    print(
        f"LOINC enrichment: wrote {len(part_tests)} Parts and {len(input_terms)} linked Terms "
        f"to {OUT_LOINC_PART_TESTS}."
    )
    print(
        "LOINC enrichment: local-first labels — "
        f"MRCONSO Parts {label_stats['mrconso_parts']}/{label_stats['total_parts']}, "
        f"Part-link CSV fallback {label_stats['part_link_fallbacks']}, "
        f"core CSV fallback {label_stats['core_part_fallbacks']}; "
        f"MRCONSO Terms {label_stats['mrconso_terms']}/{label_stats['total_terms']}, "
        f"Part-link CSV fallback {label_stats['term_link_fallbacks']}; "
        f"unresolved Parts {len(missing_parts)}, Terms {len(missing_terms)}; "
        f"output={OUT_LOINC_LABELS}."
    )
    return missing_parts, missing_terms

def norm_hpo(x: str) -> Optional[str]:
    s = (x or "").strip()
    if not s:
        return None

    # URI form
    if "purl.obolibrary.org/obo/" in s:
        s = s.rsplit("/", 1)[-1]

    # Canonicalize separators
    s = s.replace("_", ":").replace("-", ":").upper()
    s = re.sub(r"\s+", "", s)

    # Exact forms
    m = re.fullmatch(r"HP:(\d{7})", s)
    if m:
        return f"HP:{m.group(1)}"
    m = re.fullmatch(r"HP(\d{7})", s)
    if m:
        return f"HP:{m.group(1)}"
    m = re.fullmatch(r"(\d{7})", s)
    if m:
        return f"HP:{m.group(1)}"

    # Fallback: look for HP + 7 digits anywhere in the string
    m = re.search(r"(?i)HP\s*[:_\- ]?\s*(\d{7})", x or "")
    if m:
        return f"HP:{m.group(1)}"
    return None

def extract_hpo_ids(cell: str) -> List[str]:
    """Extract one or multiple HPO IDs from a free-text cell (supports '|', ',', ';', URI, etc.)."""
    raw = (cell or "").strip()
    if not raw:
        return []
    found = []

    # First pass: split on common separators
    for tok in re.split(r"[|,;\r\n]+", raw):
        tok = tok.strip()
        if not tok:
            continue
        h = norm_hpo(tok)
        if h:
            found.append(h)

    # Second pass: regex scan over full cell (captures embedded IDs)
    for mm in re.finditer(r"(?i)HP\s*[:_\- ]?\s*(\d{7})", raw):
        found.append(f"HP:{mm.group(1)}")

    # de-duplicate while preserving order
    uniq = []
    seen = set()
    for h in found:
        if h not in seen:
            seen.add(h)
            uniq.append(h)
    return uniq


# Resolvers
def uniprot_name(acc: str) -> Tuple[Optional[str], str]:
    page = f"https://www.uniprot.org/uniprotkb/{acc}"
    r = req_get(UNIPROT_JSON_URL.format(acc=acc))
    if not r or r.status_code != 200: return None, page
    try:
        d = r.json()
        pd = d.get("proteinDescription") or {}
        name = ((((pd.get("recommendedName") or {}).get("fullName") or {}).get("value")) or 
                next((x.get("fullName", {}).get("value") for x in pd.get("submissionNames", []) if x.get("fullName")), None) or 
                next((x.get("fullName", {}).get("value") for x in pd.get("alternativeNames", []) if x.get("fullName")), None) or d.get("uniProtkbId"))
        return name, page
    except: return None, page

def orpha_name(orpha_id: str) -> Tuple[Optional[str], str]:
    page = f"https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert={orpha_id}"
    def _check(params):
        r = req_get(OLS4_TERM_API.format(onto="ordo"), params=params)
        if r and r.status_code == 200:
            try:
                terms = (r.json().get("_embedded") or {}).get("terms") or []
                if terms and terms[0].get("label"): return terms[0].get("label")
            except: pass
        return None
    name = _check({"short_form": f"Orphanet_{orpha_id}"})
    if not name: name = _check({"iri": f"http://www.orpha.net/ORDO/Orphanet_{orpha_id}"})
    if not name: name = _check({"obo_id": f"Orphanet:{orpha_id}"})
    return name, page

def chebi_name(num: str) -> Tuple[Optional[str], str]:
    """Resolve a ChEBI label using the dedicated ChEBI REST API first."""
    chebi_id = f"CHEBI:{num}"
    page = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}"
    endpoint = CHEBI_COMPOUND_API.format(chebi_id=quote(chebi_id, safe=":"))
    response = req_get(endpoint)
    if response and response.status_code == 200:
        try:
            payload = response.json()
            name = payload.get("name") if isinstance(payload, dict) else None
            if isinstance(name, str) and _clean_label(name):
                return _clean_label(name), page
        except (ValueError, TypeError, AttributeError):
            pass

    # Secondary fallback for transitional or ontology-indexed records.
    response = req_get(
        OLS4_TERM_API.format(onto="chebi"), params={"obo_id": chebi_id}
    )
    if response and response.status_code == 200:
        try:
            terms = (response.json().get("_embedded") or {}).get("terms") or []
            label = terms[0].get("label") if terms else None
            if isinstance(label, str) and _clean_label(label):
                return _clean_label(label), page
        except (ValueError, TypeError, AttributeError, IndexError):
            pass
    return None, page

def write_report(
    report_path: str,
    runtime_seconds: float,
    used_umls_sabs: Counter,
    failed_ids: dict,
    results_tracker_snapshot: dict,
    output_csv: str,
    output_orpha_links: str,
    output_orpha_umls_mappings: str,
    output_loinc_part_tests: str,
    output_loinc_labels: str,
    missing_loinc_part_labels: list[str],
    missing_loinc_term_labels: list[str],
    loinc_label_stats: dict[str, int],
):
    ts = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    lines = []
    lines.append("# Enrichment report")
    lines.append("")
    lines.append(f"- Run timestamp: {ts}")
    lines.append(f"- Runtime (seconds): {runtime_seconds:.2f}")
    lines.append("")
    lines.append("## Inputs and outputs")
    lines.append("")
    lines.append(f"- MRCONSO.RRF: `{IN_PATH}`")
    lines.append(f"- Core CSV: `{INPUT_CSV_CORE}`")
    lines.append(f"- Orphanet XML: `{XML_PATH}`")
    lines.append(f"- Output (Orphanet-HPO links): `{output_orpha_links}`")
    lines.append(f"- Output (ORPHA-UMLS mappings): `{output_orpha_umls_mappings}`")
    lines.append(f"- LOINC Part links: `{LOINC_PART_LINK_CSV}`")
    lines.append(f"- Output (code names): `{output_csv}`")
    lines.append(f"- Output (LOINC Part→Term dictionary): `{output_loinc_part_tests}`")
    lines.append(f"- Output (LOINC labels): `{output_loinc_labels}`")
    lines.append("")
    lines.append("## LOINC dictionary generation")
    lines.append("")
    lines.append(
        "- Label priority: MRCONSO `SAB=LNC`, then "
        "`LoincPartLink_Primary.csv`; core CSV is a final fallback for Part labels"
    )
    lines.append(
        f"- MRCONSO labels — Parts: {loinc_label_stats['mrconso_parts']}/"
        f"{loinc_label_stats['total_parts']}; Terms: "
        f"{loinc_label_stats['mrconso_terms']}/{loinc_label_stats['total_terms']}"
    )
    lines.append(
        f"- Local CSV fallbacks used — Part-link Part labels: "
        f"{loinc_label_stats['part_link_fallbacks']}; core Part labels: "
        f"{loinc_label_stats['core_part_fallbacks']}; Part-link Term labels: "
        f"{loinc_label_stats['term_link_fallbacks']}"
    )
    lines.append(
        f"- Unresolved after all local sources — Parts: "
        f"{len(missing_loinc_part_labels)}; Terms: {len(missing_loinc_term_labels)}"
    )
    if missing_loinc_part_labels:
        lines.append("")
        lines.append("### LOINC Parts without any local label")
        lines.append("")
        lines.append("```")
        lines.extend(missing_loinc_part_labels[:200])
        if len(missing_loinc_part_labels) > 200:
            lines.append(f"... ({len(missing_loinc_part_labels) - 200} more)")
        lines.append("```")
    if missing_loinc_term_labels:
        lines.append("")
        lines.append("### LOINC Terms without any local label")
        lines.append("")
        lines.append("```")
        lines.extend(missing_loinc_term_labels[:200])
        if len(missing_loinc_term_labels) > 200:
            lines.append(f"... ({len(missing_loinc_term_labels) - 200} more)")
        lines.append("```")
    lines.append("")
    lines.append("## Vocabulary policy (UMLS CUI naming)")
    lines.append("")
    lines.append("- UMLS label selection: every MRCONSO source; English, unsuppressed atoms only")
    lines.append("")
    if used_umls_sabs:
        lines.append("### UMLS vocabularies actually used (SAB → #CUIs named)")
        lines.append("")
        for sab, n in sorted(used_umls_sabs.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- {sab}: {n}")
        lines.append("")
    else:
        lines.append("### UMLS vocabularies actually used")
        lines.append("")
        lines.append("- (none)")
        lines.append("")

    lines.append("## Resolution summary")
    lines.append("")
    lines.append("| Source | Success | Failed | Total |")
    lines.append("|---|---:|---:|---:|")
    for src, stats in results_tracker_snapshot.items():
        total = stats["success"] + stats["fail"]
        lines.append(f"| {src} | {stats['success']} | {stats['fail']} | {total} |")
    lines.append("")

    lines.append("## Missing names (by source)")
    lines.append("")
    for src, ids in failed_ids.items():
        if not ids:
            continue
        lines.append(f"### {src} ({len(ids)} identifiers)")
        lines.append("")
        # Keep the report compact
        max_show = 200
        shown = ids[:max_show]
        lines.append("```")
        lines.extend(shown)
        if len(ids) > max_show:
            lines.append(f"... ({len(ids) - max_show} more)")
        lines.append("```")
        lines.append("")

    os.makedirs(os.path.dirname(report_path) or ".", exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def enrich_data(
    missing_loinc_part_labels: list[str],
    missing_loinc_term_labels: list[str],
    loinc_label_stats: dict[str, int],
):
    print("Starting enrichment process...")
    if not os.path.isfile(INPUT_CSV_CORE):
        raise FileNotFoundError(f"Required core input not found: {INPUT_CSV_CORE}")

    uni, cui_list, chebi_list = [], [], []
    orpha_from_dis, umls_from_dis = [], []
    hpo_from_core = []
    loinc_map = {}

    with open(INPUT_CSV_CORE, newline="", encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        cols = {re.sub(r"[\s\ufeff]+", "", c.lower()): c for c in (rdr.fieldnames or []) if c is not None}
        for row in rdr:
            if "uniprot_id" in cols: uni.extend(split_items(row[cols["uniprot_id"]]))
            if "umls_id" in cols: cui_list.extend(split_items(row[cols["umls_id"]]))
            if "chebi_id" in cols: chebi_list.extend(split_items(row[cols["chebi_id"]]))
            if "disease_id" in cols:
                for t in split_items(row[cols["disease_id"]]):
                    if o := norm_orpha(t): orpha_from_dis.append(o)
                    elif c := norm_umls(t): umls_from_dis.append(c)
            if "hpo_id" in cols and row.get(cols["hpo_id"]):
                for h in extract_hpo_ids(row[cols["hpo_id"]]):
                    hpo_from_core.append(h)
            if "loinc_part_id" in cols and row.get(cols["loinc_part_id"]):
                raw_ids = split_items_pipe(row[cols["loinc_part_id"]])
                nm_col = cols.get("loinc_part")
                names = split_items_pipe(row.get(nm_col, "")) if nm_col else []
                # Pair IDs and labels before discarding malformed IDs, so an
                # invalid earlier ID cannot shift labels onto later valid IDs.
                for i, raw_id in enumerate(raw_ids):
                    lid = norm_loinc_part(raw_id)
                    if not lid:
                        continue
                    nm = names[i] if i < len(names) else ""
                    if lid not in loinc_map or not loinc_map[lid]:
                        loinc_map[lid] = nm

    # Deduplicate IDs
    uni_ids = sorted({v for v in (norm_uniprot(x) for x in uni) if v})
    cui_ids = sorted({v for v in (norm_umls(x) for x in cui_list) if v} | set(umls_from_dis))
    chebi_ids = sorted({v for v in (norm_chebi(x) for x in chebi_list) if v})
    orpha_ids = sorted(set(orpha_from_dis))
    hpo_ids_of_interest = set(hpo_from_core)
    if os.path.exists(OUT_ORPHA_LINKS):
        df_orpha = pd.read_csv(OUT_ORPHA_LINKS)
        if "orpha_code" in df_orpha.columns:
            orpha_ids = sorted(set(orpha_ids) | set(df_orpha["orpha_code"].astype(str).str.split("_").str[-1]))
        if "HPOId" in df_orpha.columns:
            for raw in df_orpha["HPOId"].dropna().astype(str):
                hid = norm_hpo(raw)
                if hid:
                    hpo_ids_of_interest.add(hid)
    hpo_ids = sorted(hpo_ids_of_interest)
    print(f"Parsed {len(hpo_from_core)} HPO mentions from core CSV -> {len(set(hpo_from_core))} unique normalized IDs.")

    out_rows = []

    # Track identifiers with no resolved name, per source
    failed_ids = {k: [] for k in results_tracker.keys()}

    # UMLS local lookup
    print(f"Resolving {len(cui_ids)} UMLS IDs (local memory)...")
    used_umls_sabs = Counter()
    for c in cui_ids:
        name = umls_names_map.get(c)
        if name:
            results_tracker["UMLS"]["success"] += 1
            meta = umls_name_meta.get(c) or {}
            if meta.get("sab"):
                used_umls_sabs[meta["sab"]] += 1
        else:
            results_tracker["UMLS"]["fail"] += 1
            failed_ids["UMLS"].append(c)
            print(f"[UNSUCCESSFUL] UMLS ID: {c} - No English, unsuppressed label found in MRCONSO.")
        out_rows.append({"source": "UMLS", "id": c, "name": name or "", "url": f"https://uts.nlm.nih.gov/uts/umls/concept/{c}"})

    # HPO local lookup (from MRCONSO HPO entries)
    print(f"Resolving {len(hpo_ids)} HPO IDs (local memory)...")
    for hid in hpo_ids:
        name = hpo_names_map.get(hid)
        if name:
            results_tracker["HPO"]["success"] += 1
        else:
            results_tracker["HPO"]["fail"] += 1
            failed_ids["HPO"].append(hid)
            print(f"[UNSUCCESSFUL] HPO ID: {hid} - Not found in MRCONSO HPO entries.")
        out_rows.append({
            "source": "HPO",
            "id": hid,
            "name": name or "",
            "url": f"http://purl.obolibrary.org/obo/{hid.replace(':', '_')}",
        })

    # sanity check: ensure all extracted HPO IDs were written
    written_hpo_ids = {r["id"] for r in out_rows if r["source"] == "HPO"}
    missing_written = sorted(set(hpo_ids) - written_hpo_ids)
    if missing_written:
        print(f"[WARNING] {len(missing_written)} HPO IDs were parsed but not written to output. Example: {missing_written[:5]}")

    # Local-first external-resource resolution. APIs are called only when no
    # label is available from the configured local files.
    print("Resolving external IDs (local files first; APIs only for missing labels)...")

    pending_api: list[tuple[str, str]] = []
    for oid in orpha_ids:
        local_name = orpha_names_map.get(oid) or orpha_xml_names_map.get(oid)
        if local_name:
            results_tracker["ORPHA"]["success"] += 1
            out_rows.append({
                "source": "ORPHA",
                "id": oid,
                "name": local_name,
                "url": f"https://www.orpha.net/en/disease/detail/{oid}",
            })
        else:
            pending_api.append(("ORPHA", oid))

    for chebi_id in chebi_ids:
        local_name = chebi_names_map.get(chebi_id)
        raw_id = f"CHEBI:{chebi_id}"
        if local_name:
            results_tracker["ChEBI"]["success"] += 1
            out_rows.append({
                "source": "ChEBI",
                "id": raw_id,
                "name": local_name,
                "url": f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={raw_id}",
            })
        else:
            pending_api.append(("ChEBI", raw_id))

    # No local UniProt label source is configured, so all UniProt identifiers
    # remain API candidates.
    pending_api.extend(("UniProt", acc) for acc in uni_ids)

    if pending_api:
        print(f"API fallback required for {len(pending_api)} identifiers.")
        with ThreadPoolExecutor(max_workers=10) as executor:
            future_to_item = {}
            for source, raw_id in pending_api:
                if source == "UniProt":
                    future = executor.submit(uniprot_name, raw_id)
                elif source == "ORPHA":
                    future = executor.submit(orpha_name, raw_id)
                elif source == "ChEBI":
                    future = executor.submit(chebi_name, raw_id.removeprefix("CHEBI:"))
                else:  # pragma: no cover - internal invariant
                    raise RuntimeError(f"Unsupported API fallback source: {source}")
                future_to_item[future] = (source, raw_id)

            for future in as_completed(future_to_item):
                source, raw_id = future_to_item[future]
                try:
                    name, url = future.result()
                    if name:
                        results_tracker[source]["success"] += 1
                    else:
                        results_tracker[source]["fail"] += 1
                        failed_ids[source].append(raw_id)
                        print(f"[UNSUCCESSFUL] {source} ID: {raw_id} - API returned no name.")
                    out_rows.append({"source": source, "id": raw_id, "name": name or "", "url": url})
                except Exception as e:
                    results_tracker[source]["fail"] += 1
                    failed_ids[source].append(raw_id)
                    print(f"[ERROR] {source} ID: {raw_id} - {e}")
                    if source == "UniProt":
                        fallback_url = f"https://www.uniprot.org/uniprotkb/{raw_id}"
                    elif source == "ORPHA":
                        fallback_url = f"https://www.orpha.net/en/disease/detail/{raw_id}"
                    elif source == "ChEBI":
                        fallback_url = (
                            "https://www.ebi.ac.uk/chebi/searchId.do"
                            f"?chebiId={raw_id}"
                        )
                    else:
                        fallback_url = ""
                    out_rows.append({
                        "source": source,
                        "id": raw_id,
                        "name": "",
                        "url": fallback_url,
                    })
    else:
        print("No API fallback was required.")

    # LOINC Part labels use the same final local-first values written to
    # loinc_labels.json, so code_names.csv and the KG input cannot disagree.
    for lid in sorted(loinc_map, key=lambda value: (value.casefold(), value)):
        name = loinc_part_names_map.get(lid)
        if name:
            results_tracker["LOINC"]["success"] += 1
        else:
            results_tracker["LOINC"]["fail"] += 1
            failed_ids["LOINC"].append(lid)
            print(
                f"[UNSUCCESSFUL] LOINC ID: {lid} - No label found in "
                "MRCONSO, LoincPartLink_Primary.csv, or the core CSV."
            )
        out_rows.append(
            {
                "source": "LOINC",
                "id": lid,
                "name": name or "",
                "url": f"https://loinc.org/{lid}",
            }
        )

    # Stable output independent of concurrent request completion order.
    out_rows.sort(
        key=lambda row: (
            row["source"].casefold(), row["source"],
            row["id"].casefold(), row["id"],
            row["name"].casefold(), row["name"],
            row["url"],
        )
    )
    for ids in failed_ids.values():
        ids.sort(key=lambda value: (value.casefold(), value))

    os.makedirs(os.path.dirname(OUTPUT_CSV_FINAL) or ".", exist_ok=True)
    temporary = OUTPUT_CSV_FINAL + ".tmp"
    try:
        with open(temporary, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["source", "id", "name", "url"], lineterminator="\n")
            writer.writeheader()
            writer.writerows(out_rows)
        os.replace(temporary, OUTPUT_CSV_FINAL)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)

    # PRINT FINAL SUMMARY
    print("\n" + "="*45)
    print(f"{'ENRICHMENT SUMMARY':^45}")
    print("="*45)
    print(f"{'Source':<12} | {'Success':<10} | {'Failed':<10} | {'Total':<10}")
    print("-" * 45)
    for src, stats in results_tracker.items():
        total = stats["success"] + stats["fail"]
        print(f"{src:<12} | {stats['success']:<10} | {stats['fail']:<10} | {total:<10}")
    print("="*45)

    # Write a machine-readable report (Markdown) capturing the run and missing identifiers
    runtime_seconds = time.time() - start_t if "start_t" in globals() else 0.0
    write_report(
        report_path=REPORT_PATH,
        runtime_seconds=runtime_seconds,
        used_umls_sabs=used_umls_sabs,
        failed_ids=failed_ids,
        results_tracker_snapshot=results_tracker,
        output_csv=OUTPUT_CSV_FINAL,
        output_orpha_links=OUT_ORPHA_LINKS,
        output_orpha_umls_mappings=OUT_ORPHA_UMLS_MAPPINGS,
        output_loinc_part_tests=OUT_LOINC_PART_TESTS,
        output_loinc_labels=OUT_LOINC_LABELS,
        missing_loinc_part_labels=missing_loinc_part_labels,
        missing_loinc_term_labels=missing_loinc_term_labels,
        loinc_label_stats=loinc_label_stats,
    )
    print(f"Report written to {REPORT_PATH}")

    # Short on-screen report of missing identifiers (kept compact)
    any_missing = any(bool(v) for v in failed_ids.values())
    if any_missing:
        print("\nMissing identifiers (counts; first 10 shown per source):")
        for src, ids in failed_ids.items():
            if not ids:
                continue
            preview = ", ".join(ids[:10])
            more = f" (+{len(ids) - 10} more)" if len(ids) > 10 else ""
            print(f"- {src}: {len(ids)} [{preview}]{more}")

    if used_umls_sabs:
        print("\nUMLS vocabularies actually used for naming (SAB → #CUIs):")
        for sab, n in sorted(used_umls_sabs.items(), key=lambda item: (-item[1], item[0])):
            print(f"- {sab}: {n}")

def reset_global_state() -> None:
    umls_names_map.clear()
    umls_name_meta.clear()
    hpo_names_map.clear()
    orpha_names_map.clear()
    orpha_xml_names_map.clear()
    chebi_names_map.clear()
    loinc_part_names_map.clear()
    loinc_term_names_map.clear()
    for stats in results_tracker.values():
        stats["success"] = 0
        stats["fail"] = 0


def validate_required_inputs() -> None:
    missing = [
        path for path in (IN_PATH, INPUT_CSV_CORE, XML_PATH, LOINC_PART_LINK_CSV)
        if not os.path.isfile(path)
    ]
    if missing:
        raise FileNotFoundError(
            "Required enrichment inputs are missing:\n  " + "\n  ".join(missing)
        )


def main():
    global start_t
    start_t = time.time()
    reset_global_state()
    validate_required_inputs()

    input_cuis, input_orphas = collect_core_orpha_umls_identifiers()
    input_chebi_ids = collect_core_chebi_identifiers()
    input_loinc_parts, core_loinc_part_labels = collect_core_loinc_parts_and_labels()
    (
        loinc_part_tests,
        input_loinc_terms,
        link_loinc_part_labels,
        link_loinc_term_labels,
    ) = parse_loinc_primary_links(input_loinc_parts)
    orphanet_hpo_ids = process_orphanet(input_orphas)

    # Include HPO identifiers directly present in the core before scanning MRCONSO.
    core_hpo_ids: set[str] = set()
    with open(INPUT_CSV_CORE, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = {
            re.sub(r"[\s\ufeff]+", "", column.lower()): column
            for column in (reader.fieldnames or [])
            if column is not None
        }
        hpo_column = columns.get("hpo_id")
        if hpo_column:
            for row in reader:
                core_hpo_ids.update(extract_hpo_ids(row.get(hpo_column, "")))

    process_mrconso(
        input_cuis,
        input_orphas,
        input_chebi_ids,
        orphanet_hpo_ids | core_hpo_ids,
        input_loinc_parts,
        input_loinc_terms,
    )
    loinc_label_stats = apply_loinc_local_label_fallbacks(
        input_loinc_parts,
        input_loinc_terms,
        link_loinc_part_labels,
        link_loinc_term_labels,
        core_loinc_part_labels,
    )
    missing_loinc_part_labels, missing_loinc_term_labels = (
        write_loinc_enrichment_files(
            loinc_part_tests,
            input_loinc_parts,
            input_loinc_terms,
            loinc_label_stats,
        )
    )
    enrich_data(
        missing_loinc_part_labels,
        missing_loinc_term_labels,
        loinc_label_stats,
    )
    print(f"\nPipeline completed in {time.time() - start_t:.2f} seconds.")


if __name__ == "__main__":
    main()