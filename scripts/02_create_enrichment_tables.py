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
from typing import Optional, Tuple, List

from concurrent.futures import ThreadPoolExecutor, as_completed

# --- CONFIGURATION ---
# NOTE: Update IN_PATH to point to your local MRCONSO.RRF file
IN_PATH = "/mnt/d/umls-2024AB-full_metamor/2024AB-full/2024AB/2024AB/META/MRCONSO.RRF"

# File Paths
DATA_DIR = "../data"
ENRICH_DIR = os.path.join(DATA_DIR, "enrichment_tables")

# Final output files
XML_PATH = os.path.join(DATA_DIR, "en_product4.xml")
OUT_ORPHA_LINKS = os.path.join(ENRICH_DIR, "orphanet_hpo_links.csv")
INPUT_CSV_CORE = os.path.join(DATA_DIR, "makaao_core.csv")
OUTPUT_CSV_FINAL = os.path.join(ENRICH_DIR, "code_names.csv")

# API URLs
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{acc}"
OLS4_TERM_API = "https://www.ebi.ac.uk/ols4/api/ontologies/{onto}/terms"
HEADERS = {"Accept": "application/json"}

# Global storage for in-memory mapping and tracking
umls_names_map = {}
results_tracker = {
    "UniProt": {"success": 0, "fail": 0},
    "ORPHA": {"success": 0, "fail": 0},
    "ChEBI": {"success": 0, "fail": 0},
    "UMLS": {"success": 0, "fail": 0},
    "LOINC": {"success": 0, "fail": 0}
}

# Increase CSV field size limit
_lim = sys.maxsize
while True:
    try:
        csv.field_size_limit(_lim)
        break
    except OverflowError:
        _lim //= 10


# --- PART 1: PROCESS MRCONSO.RRF (Direct to Memory) ---
def process_mrconso():
    print(f"Processing MRCONSO.RRF from {IN_PATH}...")
    if not os.path.exists(IN_PATH):
        print(f"ERROR: Input file not found at {IN_PATH}.")
        return

    os.makedirs(ENRICH_DIR, exist_ok=True)
    count = 0
    with open(IN_PATH, "r", encoding="utf-8", errors="replace", newline="") as fin:
        reader = csv.reader(fin, delimiter="|")
        for parts in reader:
            if len(parts) < 18:
                continue
            # Filter for English Preferred names
            if parts[1] == "ENG" and parts[2] == "P" and parts[16] == "N":
                cui, name = parts[0], parts[14]
                if cui not in umls_names_map:
                    umls_names_map[cui] = name
                    count += 1
    print(f"Part 1: Loaded {count} UMLS names into memory.")


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

def process_orphanet():
    print(f"Processing Orphanet XML from {XML_PATH}...")
    if not os.path.exists(XML_PATH):
        print(f"Warning: XML file not found at {XML_PATH}. Skipping Step 2.")
        return
    rows = []
    try:
        root = ET.parse(XML_PATH).getroot()
        for disorder in root.iter():
            if local_tag(disorder.tag) != "Disorder": continue
            oc, assoc_list = None, None
            for ch in disorder:
                n = local_tag(ch.tag)
                if n == "OrphaCode" and (ch.text or "").strip(): oc = ch.text.strip()
                elif n == "HPODisorderAssociationList": assoc_list = ch
            if not oc or assoc_list is None: continue
            for assoc in assoc_list:
                hpo_id, hpo_term = None, None
                hpo = next((n for n in assoc if local_tag(n.tag) == "HPO"), None)
                if hpo is not None:
                    hid = next((n for n in hpo if local_tag(n.tag) == "HPOId"), None)
                    htm = next((n for n in hpo if local_tag(n.tag) == "HPOTerm"), None)
                    if hid is not None: hpo_id = (hid.text or "").strip()
                    if htm is not None: hpo_term = (htm.text or "").strip()
                if not hpo_id: continue
                freq_name = None
                freq = next((n for n in assoc if local_tag(n.tag) == "HPOFrequency"), None)
                if freq is not None:
                    names = [n for n in freq if local_tag(n.tag) == "Name" and (n.text or "").strip()]
                    en = next((n for n in names if (get_lang(n) or "").lower() == "en"), None)
                    chosen = en or (names[0] if names else None)
                    if chosen is not None: freq_name = chosen.text.strip()
                rank, freq_name = freq_rank(freq_name)
                if rank >= 3:
                    rows.append(("http://www.orpha.net/ORDO/Orphanet_" + oc, 
                                 "http://purl.obolibrary.org/obo/" + hpo_id.replace(":", "_"), 
                                 hpo_term, freq_name, rank))
        df = pd.DataFrame(rows, columns=["orpha_code", "HPOId", "HPOTerm", "frequency", "rank"])
        if not df.empty:
            df = df.sort_values("rank", ascending=False).drop_duplicates(["orpha_code", "HPOId"]).drop(columns=["rank"])
        df.to_csv(OUT_ORPHA_LINKS, index=False)
        print(f"Part 2: Wrote {len(df)} rows to {OUT_ORPHA_LINKS}")
    except Exception as e:
        print(f"Error in Orphanet processing: {e}")


# --- PART 3: ENRICHMENT & TRACKING ---
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
    m = re.fullmatch(r"(?i)(?:ORPHA:|ORPHANET:)(\d+)", (x or "").strip())
    return m.group(1) if m else None

def norm_chebi(x: str) -> Optional[str]:
    m = re.fullmatch(r"(?i)(?:CHEBI:|CHE:)?(\d+)", (x or "").strip())
    return m.group(1) if m else None

def norm_loinc_part(x: str) -> Optional[str]:
    x = (x or "").strip().upper()
    m = re.fullmatch(r"LP\d+(?:-\d+)?", x)
    return m.group(0) if m else None

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
    page = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{num}"
    r = req_get(OLS4_TERM_API.format(onto="chebi"), params={"obo_id": f"CHEBI:{num}"})
    if not r or r.status_code != 200: return None, page
    try:
        terms = (r.json().get("_embedded") or {}).get("terms") or []
        return (terms[0].get("label") if terms else None), page
    except: return None, page

def enrich_data():
    print("Starting enrichment process...")
    if not os.path.exists(INPUT_CSV_CORE):
        print(f"ERROR: Input CSV {INPUT_CSV_CORE} not found.")
        return

    uni, cui_list, chebi_list = [], [], []
    orpha_from_dis, umls_from_dis = [], []
    loinc_map = {}

    with open(INPUT_CSV_CORE, newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        cols = {c.lower(): c for c in (rdr.fieldnames or [])}
        for row in rdr:
            if "uniprot_id" in cols: uni.extend(split_items(row[cols["uniprot_id"]]))
            if "umls_id" in cols: cui_list.extend(split_items(row[cols["umls_id"]]))
            if "chebi_id" in cols: chebi_list.extend(split_items(row[cols["chebi_id"]]))
            if "disease_id" in cols:
                for t in split_items(row[cols["disease_id"]]):
                    if o := norm_orpha(t): orpha_from_dis.append(o)
                    elif c := norm_umls(t): umls_from_dis.append(c)
            if "loinc_part_id" in cols and row.get(cols["loinc_part_id"]):
                ids = [norm_loinc_part(x) for x in split_items_pipe(row[cols["loinc_part_id"]])]
                nm_col = cols.get("loinc_part")
                names = split_items_pipe(row.get(nm_col, "")) if nm_col else []
                for i, lid in enumerate([x for x in ids if x]):
                    nm = names[i] if i < len(names) else ""
                    if lid not in loinc_map or not loinc_map[lid]: loinc_map[lid] = nm

    # Deduplicate IDs
    uni_ids = sorted({v for v in (norm_uniprot(x) for x in uni) if v})
    cui_ids = sorted({v for v in (norm_umls(x) for x in cui_list) if v} | set(umls_from_dis))
    chebi_ids = sorted({v for v in (norm_chebi(x) for x in chebi_list) if v})
    orpha_ids = sorted(set(orpha_from_dis))
    if os.path.exists(OUT_ORPHA_LINKS):
        df_orpha = pd.read_csv(OUT_ORPHA_LINKS)
        orpha_ids = sorted(set(orpha_ids) | set(df_orpha['orpha_code'].str.split('_').str[-1]))

    out_rows = []

    # UMLS local lookup
    print(f"Resolving {len(cui_ids)} UMLS IDs (local memory)...")
    for c in cui_ids:
        name = umls_names_map.get(c)
        if name: results_tracker["UMLS"]["success"] += 1
        else:
            results_tracker["UMLS"]["fail"] += 1
            print(f"[UNSUCCESSFUL] UMLS ID: {c} - Not found in MRCONSO map.")
        out_rows.append({"source": "UMLS", "id": c, "name": name or "", "url": f"https://uts.nlm.nih.gov/uts/umls/concept/{c}"})

    # Concurrent remote lookup
    print("Resolving external IDs concurrently...")
    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_item = {}
        for acc in uni_ids: future_to_item[executor.submit(uniprot_name, acc)] = ("UniProt", acc)
        for oid in orpha_ids: future_to_item[executor.submit(orpha_name, oid)] = ("ORPHA", oid)
        for ch in chebi_ids: future_to_item[executor.submit(chebi_name, ch)] = ("ChEBI", f"CHEBI:{ch}")

        for future in as_completed(future_to_item):
            source, raw_id = future_to_item[future]
            try:
                name, url = future.result()
                if name: results_tracker[source]["success"] += 1
                else:
                    results_tracker[source]["fail"] += 1
                    print(f"[UNSUCCESSFUL] {source} ID: {raw_id} - API returned no name.")
                out_rows.append({"source": source, "id": raw_id, "name": name or "", "url": url})
            except Exception as e:
                results_tracker[source]["fail"] += 1
                print(f"[ERROR] {source} ID: {raw_id} - {e}")

    # LOINC local resolution
    for lid, name in loinc_map.items():
        if name: results_tracker["LOINC"]["success"] += 1
        else:
            results_tracker["LOINC"]["fail"] += 1
            print(f"[UNSUCCESSFUL] LOINC ID: {lid} - No name found in core CSV.")
        out_rows.append({"source": "LOINC", "id": lid, "name": name or "", "url": f"https://loinc.org/{lid}"})

    # Write output
    os.makedirs(os.path.dirname(OUTPUT_CSV_FINAL) or ".", exist_ok=True)
    with open(OUTPUT_CSV_FINAL, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["source", "id", "name", "url"])
        w.writeheader()
        w.writerows(out_rows)

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

def main():
    start_t = time.time()
    process_mrconso()
    process_orphanet()
    enrich_data()
    print(f"\nPipeline completed in {time.time() - start_t:.2f} seconds.")

if __name__ == "__main__":
    main()