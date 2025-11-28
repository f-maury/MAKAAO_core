#!/usr/bin/env python3
import sys, csv
import pandas as pd
import xml.etree.ElementTree as ET
import csv, json, os, re, time
from typing import Dict, Optional, Tuple, List
import requests


UMLS_API_KEY = "67bd1b8b-87f0-40de-bf20-255e6f1721a3"  # WRITE YOUR UMLS API KEY HERE

# allow very long STR fields
_lim = sys.maxsize
while True:
    try:
        csv.field_size_limit(_lim); break # allow processing of very long string in CSV cells
    except OverflowError:
        _lim //= 10

IN_PATH  = "/mnt/d/umls-2024AB-full_metamor/2024AB-full/2024AB/2024AB/META/MRCONSO.RRF" # path to MRCONSO.RRF we will use (ideally ../data/umls/MRCONSO.RRF)

OUT_HPO  = "../data/enrichment_tables/umls_hpo.csv" # enrichment tables we will produce
OUT_ORPH = "../data/enrichment_tables/umls_orphanet.csv"
OUT_SCT  = "../data/enrichment_tables/umls_snomed.csv"
OUT_LNC  = "../data/enrichment_tables/umls_loinc.csv"
OUT_CHE  = "../data/enrichment_tables/umls_chebi.csv"
OUT_NAME = "../data/enrichment_tables/umls_names.csv"

HEADER = ["CUI","LAT","TS","LUI","STT","SUI","ISPREF","AUI","SAUI","SCUI",
          "SDUI","SAB","TTY","CODE","STR","SRL","SUPPRESS","CVF"] # columns in mrconso

with open(IN_PATH, "r", encoding="utf-8", errors="replace", newline="") as fin, \
     open(OUT_HPO,  "w", encoding="utf-8", newline="") as fhpo, \
     open(OUT_ORPH, "w", encoding="utf-8", newline="") as forph, \
     open(OUT_SCT,  "w", encoding="utf-8", newline="") as fsct, \
     open(OUT_LNC,  "w", encoding="utf-8", newline="") as flnc, \
     open(OUT_CHE,  "w", encoding="utf-8", newline="") as fche, \
     open(OUT_NAME, "w", encoding="utf-8", newline="") as fnames:
    reader = csv.reader(fin, delimiter="|")
    whpo   = csv.writer(fhpo,   lineterminator="\n") # open files to write, and file to read
    worph  = csv.writer(forph,  lineterminator="\n")
    wsct   = csv.writer(fsct,   lineterminator="\n")
    wlnc   = csv.writer(flnc,   lineterminator="\n")
    wche   = csv.writer(fche,   lineterminator="\n")
    wnames = csv.writer(fnames, lineterminator="\n")

    whpo.writerow(HEADER)
    worph.writerow(HEADER)
    wsct.writerow(HEADER)
    wlnc.writerow(HEADER)
    wche.writerow(HEADER)
    wnames.writerow(["CUI","STR"])

    names_seen = set()

    for parts in reader:
        if len(parts) < 18:
            continue
        row = parts[:18]
        sab = row[11]

        if sab == "HPO":
            whpo.writerow(row)
        elif sab == "ORPHANET":
            worph.writerow(row)
        elif sab.startswith("SNOMEDCT"):
            wsct.writerow(row)
        elif sab == "LNC":                 # LOINC
            wlnc.writerow(row)
        elif sab == "CHEBI":               # CHEBI: normalize CODE to CHEBI:<num>
            code = (row[13] or "").strip()
            if code:
                u = code.upper()
                if u.startswith("CHEBI:"):                       # already CHEBI:<num>
                    code = "CHEBI:" + code.split(":", 1)[1]
                elif u.startswith("CHE:") and not u.startswith("CHEBI:"):  # CHE:<num> -> CHEBI:<num>
                    code = "CHEBI:" + code.split(":", 1)[1]
                elif re.fullmatch(r"\d+", code):                 # bare number -> CHEBI:<num>
                    code = "CHEBI:" + code
            row[13] = code
            wche.writerow(row)

        # main English preferred, not suppressed
        if row[1] == "ENG" and row[2] == "P" and row[16] == "N":
            key = (row[0], row[14])
            if key not in names_seen:
                wnames.writerow([row[0], row[14]])
                names_seen.add(key)

###########################################################################################

xml_path = "../data/en_product4.xml" # read Orphanet file with disease-HPO links
out_csv  = "../data/enrichment_tables/orphanet_hpo_links.csv" # write in new csv file

base_url_orpha = "http://www.orpha.net/ORDO/Orphanet_"
base_url_hpo   = "http://purl.obolibrary.org/obo/"

def local(tag: str) -> str:
    return tag.rsplit('}', 1)[-1] if '}' in tag else tag # extract name from XML tag

def get_lang(elem):
    return elem.attrib.get("lang") or elem.attrib.get("{http://www.w3.org/XML/1998/namespace}lang") # extract xml:lang attribute or lang attribute from XML element

def freq_rank(txt: str) -> tuple[int, str]: # translate frequency text to numerical rank in a tuple, if nothing found return a tuple with -1
    if not txt: return (-1, "Unknown")
    t = txt.strip().lower()
    if "obligate" in t or "always present" in t or "100%" in t: return (5, txt.strip())
    if "very frequent" in t or ("99" in t and "80" in t): return (4, txt.strip())
    if "frequent" in t and "very" not in t: return (3, txt.strip())
    return (-1, txt.strip())

rows = []  # (orpha_iri, hpo_iri, hpo_term, frequency, rank)

root = ET.parse(xml_path).getroot()
for disorder in root.iter():
    if local(disorder.tag) != "Disorder": # find Disorders in the XML structure
        continue

    oc, assoc_list = None, None
    for ch in disorder:
        n = local(ch.tag)
        if n == "OrphaCode" and (ch.text or "").strip(): # get orpha code
            oc = ch.text.strip()
        elif n == "HPODisorderAssociationList": # get HPO associations
            assoc_list = ch
    if not oc or assoc_list is None:
        continue

    for assoc in assoc_list:
        if local(assoc.tag) != "HPODisorderAssociation": # if no HPO association, skip
            continue

        hpo_id, hpo_term = None, None 
        hpo = next((n for n in assoc if local(n.tag) == "HPO"), None) # if there are HPO terms, for each of them, we get some info
        if hpo is not None:
            hid = next((n for n in hpo if local(n.tag) == "HPOId"), None)
            htm = next((n for n in hpo if local(n.tag) == "HPOTerm"), None)
            if hid is not None and (hid.text or "").strip():
                hpo_id = hid.text.strip()
            if htm is not None and (htm.text or "").strip():
                hpo_term = htm.text.strip()
        if not hpo_id:
            continue

        freq_name = None
        freq = next((n for n in assoc if local(n.tag) == "HPOFrequency"), None) # get the frequency of each HPO term
        if freq is not None:
            names = [n for n in freq if local(n.tag) == "Name" and (n.text or "").strip()]
            en = next((n for n in names if (get_lang(n) or "").lower() == "en"), None)
            chosen = en or (names[0] if names else None)
            if chosen is not None:
                freq_name = chosen.text.strip()

        rank, freq_name = freq_rank(freq_name)
        if rank >= 3:  # keep Frequent or higher
            rows.append((
                base_url_orpha + oc,
                base_url_hpo + hpo_id.replace(":", "_"),
                hpo_term,
                freq_name,
                rank,
            ))

df = pd.DataFrame(rows, columns=["orpha_code","HPOId","HPOTerm","frequency","rank"]) # write result in a csv file
if not df.empty:
    df = (df.sort_values("rank", ascending=False)
            .drop_duplicates(["orpha_code","HPOId"])
            .sort_values(["orpha_code","HPOId"])
            .drop(columns=["rank"]))

df.to_csv(out_csv, index=False)
print(f"Wrote {len(df):,} rows to {out_csv}")

############################################################################################


# --- config ---
INPUT_CSV   = "../data/makaao_core.csv" # read makaao_core
OUTPUT_CSV  = "../data/enrichment_tables/code_names.csv"  # write in code_names.csv; required for UMLS + SNOMED

# API URLs to retrieve info from various terminologies
UNIPROT_JSON_URL       = "https://rest.uniprot.org/uniprotkb/{acc}"
UMLS_CONCEPT_URL       = "https://uts-ws.nlm.nih.gov/rest/content/current/CUI/{cui}"
UMLS_SOURCE_CODE_URL   = "https://uts-ws.nlm.nih.gov/rest/content/current/source/{sab}/{code}"
OLS4_TERM_API          = "https://www.ebi.ac.uk/ols4/api/ontologies/{onto}/terms"
HEADERS = {"Accept": "application/json"}

SNOMED_SABS = ["SNOMEDCT_US", "SNOMEDCT", "SNOMEDCT_CORE", "SNOMEDCT_VET", "SNOMEDCT_ES", "SNOMEDCT_UK"] # select SNOMED versions we will check

# ---------- helpers ----------
def split_items(cell: str) -> List[str]: # split string around various delimiters, return a list
    """Generic split for ID fields: split on whitespace, pipe, comma, semicolon."""
    if not cell:
        return []
    return [s for s in re.split(r"[ \t\r\n|,;]+", cell.strip()) if s]

def split_items_pipe(cell: str) -> List[str]: # split string around |, return list
    """Strict split for LOINC columns: only '|' or newlines. Preserves spaces inside names."""
    if not cell:
        return []
    return [s.strip() for s in re.split(r"[|\r\n]+", cell.strip()) if s.strip()]

def req_get(url: str, params: Optional[dict] = None, headers: Optional[dict] = None, # function to query an API with retries
            retries: int = 3, backoff: float = 0.7, timeout: int = 25) -> Optional[requests.Response]:
    last = None
    for i in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers or HEADERS, timeout=timeout)
            if r.status_code in (200, 404):
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(backoff * (2 ** i)); continue
            return r
        except requests.RequestException as e:
            last = e
        time.sleep(backoff * (2 ** i))
    if last:
        raise last
    return None

# ---------- normalizers ----------
def norm_uniprot(x: str) -> Optional[str]:
    # Remove optional UniProt prefix "UP" or "UP:" (case-insensitive), then strip spaces.
    x = re.sub(r"(?i)^UP:?", "", (x or "").strip())
    # If the remaining string is 6–10 alphanumeric characters, return it uppercased;
    # otherwise return None (invalid UniProt-style ID).
    return x.upper() if re.fullmatch(r"[A-Za-z0-9]{6,10}", x) else None


def norm_umls(x: str) -> Optional[str]:
    # Match optional "CUI:" prefix, followed by a CUI of the form C + 7–8 digits (case-insensitive).
    m = re.fullmatch(r"(?i)(?:CUI:)?(C\d{7,8})", (x or "").strip())
    # If matched, return the CUI in uppercase (e.g. "C1234567"); else None.
    return m.group(1).upper() if m else None


def norm_snomed_prefixed(x: str) -> Optional[str]:
    """Accept only explicitly prefixed SNOMED codes from disease_id."""
    # Require one of the prefixes "SNOMEDCT:", "SNOMED:", or "SCTID:", case-insensitive,
    # followed by 3–18 digits. Capture only the numeric part.
    m = re.fullmatch(r"(?i)(?:SNOMEDCT:|SNOMED:|SCTID:)(\d{3,18})", (x or "").strip())
    # Return the numeric SNOMED ID (string) or None if no match.
    return m.group(1) if m else None


def norm_orpha(x: str) -> Optional[str]:
    # Require "ORPHA:" or "ORPHANET:" (case-insensitive) followed by digits; capture the digits.
    m = re.fullmatch(r"(?i)(?:ORPHA:|ORPHANET:)(\d+)", (x or "").strip())
    # Return the numeric Orphanet ID (string) or None.
    return m.group(1) if m else None


def norm_chebi(x: str) -> Optional[str]:
    # Allow optional "CHEBI:", "CHE:" prefix (case-insensitive) before digits; capture the digits.
    m = re.fullmatch(r"(?i)(?:CHEBI:|CHE:)?(\d+)", (x or "").strip())
    # Return the numeric ChEBI ID (string) or None.
    return m.group(1) if m else None


def norm_loinc_part(x: str) -> Optional[str]:
    # Normalize to uppercase and strip spaces.
    x = (x or "").strip().upper()
    # Accept LOINC part IDs of the form "LP" + digits, with optional "-digits" suffix (e.g. LP12345-6).
    m = re.fullmatch(r"LP\d+(?:-\d+)?", x)
    # Return the normalized part ID or None.
    return m.group(0) if m else None


# ---------- resolvers ----------
def uniprot_name(acc: str) -> Tuple[Optional[str], str]: # function to query Uniprot API using a Uniprot ID, return name and page URL
    r = req_get(UNIPROT_JSON_URL.format(acc=acc))
    page = f"https://www.uniprot.org/uniprotkb/{acc}"
    if not r or r.status_code != 200:
        return None, page
    try:
        d = r.json()
    except json.JSONDecodeError:
        return None, page
    pd = d.get("proteinDescription") or {}
    name = (
        (((pd.get("recommendedName") or {}).get("fullName") or {}).get("value"))
        or next((x.get("fullName", {}).get("value") for x in pd.get("submissionNames", []) if x.get("fullName")), None)
        or next((x.get("fullName", {}).get("value") for x in pd.get("alternativeNames", []) if x.get("fullName")), None)
        or d.get("uniProtkbId")
    )
    return name, page

def umls_name(cui: str) -> Tuple[Optional[str], str]: # function to get infos from UMLS API from a CUI, return name and page URL; use API key
    page = f"https://uts.nlm.nih.gov/uts/umls/concept/{cui}"
    r = req_get(UMLS_CONCEPT_URL.format(cui=cui), params={"apiKey": UMLS_API_KEY})
    if not r or r.status_code != 200:
        return None, page
    try:
        d = r.json()
    except json.JSONDecodeError:
        return None, page
    return (d.get("result") or {}).get("name"), page

def snomed_name(code: str) -> Tuple[Optional[str], str]: # function to get SNOMED name from UMLS API using SNOMED code, return name and page URL; use API key
    page = f"https://snomed.info/id/{code}"
    for sab in SNOMED_SABS:
        r = req_get(UMLS_SOURCE_CODE_URL.format(sab=sab, code=code), params={"apiKey": UMLS_API_KEY})
        if r and r.status_code == 200:
            try:
                d = r.json()
            except json.JSONDecodeError:
                continue
            name = (d.get("result") or {}).get("name")
            if name:
                return name, page
    return None, page

def _ols_label(resp_json: dict) -> Optional[str]: #extract 1st label from OLS response
    terms = (resp_json.get("_embedded") or {}).get("terms") or []
    return terms[0].get("label") if terms else None

def orpha_name(orpha_id: str) -> Tuple[Optional[str], str]: # query OLS API for Orphanet name using Orpha ID, return name and page URL
    page = f"https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert={orpha_id}"
    # Try short_form
    r = req_get(OLS4_TERM_API.format(onto="ordo"), params={"short_form": f"Orphanet_{orpha_id}"})
    if r and r.status_code == 200:
        try:
            name = _ols_label(r.json())
            if name:
                return name, page
        except json.JSONDecodeError:
            pass
    # Try IRI
    iri = f"http://www.orpha.net/ORDO/Orphanet_{orpha_id}"
    r = req_get(OLS4_TERM_API.format(onto="ordo"), params={"iri": iri})
    if r and r.status_code == 200:
        try:
            name = _ols_label(r.json())
            if name:
                return name, page
        except json.JSONDecodeError:
            pass
    # Try obo_id variants
    for obo in (f"Orphanet:{orpha_id}", f"ORPHA:{orpha_id}"):
        r = req_get(OLS4_TERM_API.format(onto="ordo"), params={"obo_id": obo})
        if r and r.status_code == 200:
            try:
                name = _ols_label(r.json())
                if name:
                    return name, page
            except json.JSONDecodeError:
                pass
    return None, page

def chebi_name(num: str) -> Tuple[Optional[str], str]: # query OLS API for ChEBI name using ChEBI ID number, return name and page URL
    page = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{num}"
    r = req_get(OLS4_TERM_API.format(onto="chebi"), params={"obo_id": f"CHEBI:{num}"})
    if not r or r.status_code != 200:
        return None, page
    try:
        d = r.json()
    except json.JSONDecodeError:
        return None, page
    return _ols_label(d), page

# ---------- main ----------
def main():
    if not UMLS_API_KEY or UMLS_API_KEY.startswith("REPLACE_"): # warning if no UMLS API key
        raise SystemExit("Set UMLS_API_KEY at top of script.")

    uni: List[str] = []
    cui: List[str] = []
    chebi: List[str] = []
    sct_from_dis: List[str] = []
    orpha_from_dis: List[str] = []
    umls_from_dis: List[str] = []
    loinc_map: Dict[str, str] = {}  # id -> name

    with open(INPUT_CSV, newline="", encoding="utf-8") as f: # read makaao_core.csv
        rdr = csv.DictReader(f)
        cols = {c.lower(): c for c in (rdr.fieldnames or [])}

        # column names
        up_col        = cols.get("uniprot_id")
        cui_col       = cols.get("umls_id")
        cheb_col      = cols.get("chebi_id")
        dis_col       = cols.get("disease_id")
        loinc_id_col  = cols.get("loinc_part_id")
        loinc_nm_col  = cols.get("loinc_part")

        for row in rdr: # parse each column of the dict read from makaao core csv
            if up_col and row.get(up_col):
                uni.extend(split_items(row[up_col]))
            if cui_col and row.get(cui_col):
                cui.extend(split_items(row[cui_col]))
            if cheb_col and row.get(cheb_col):
                chebi.extend(split_items(row[cheb_col]))

            if dis_col and row.get(dis_col): # check what type of disease ID we have
                for t in split_items(row[dis_col]):
                    sct = norm_snomed_prefixed(t)
                    if sct:
                        sct_from_dis.append(sct); continue
                    orp = norm_orpha(t)
                    if orp:
                        orpha_from_dis.append(orp); continue
                    cu = norm_umls(t)
                    if cu:
                        umls_from_dis.append(cu); continue

            if loinc_id_col and row.get(loinc_id_col): # if loinc column is there, we parse it.
                ids = [norm_loinc_part(x) for x in split_items_pipe(row[loinc_id_col])]
                ids = [x for x in ids if x]
                names = split_items_pipe(row.get(loinc_nm_col, "")) if loinc_nm_col else []
                for i, lid in enumerate(ids):
                    nm = names[i] if i < len(names) else ""
                    if lid not in loinc_map or not loinc_map[lid]: # try to map loinc id to names
                        loinc_map[lid] = nm

    # normalize and dedupe
    uni_ids    = sorted({v for v in (norm_uniprot(x) for x in uni) if v})
    cui_ids    = sorted({v for v in (norm_umls(x)    for x in cui) if v} | set(umls_from_dis))
    chebi_ids  = sorted({v for v in (norm_chebi(x)   for x in chebi) if v})
    sct_ids    = sorted(set(sct_from_dis))
    orpha_ids  = sorted(set(orpha_from_dis))
    loinc_ids  = sorted(loinc_map.keys())

    out_rows: List[Dict[str, str]] = []

    for acc in uni_ids:
        name, url = uniprot_name(acc)
        out_rows.append({"source": "UniProt", "id": acc, "name": name or "", "url": url})

    for c in cui_ids:
        name, url = umls_name(c)
        out_rows.append({"source": "UMLS", "id": c, "name": name or "", "url": url})

    for code in sct_ids:
        name, url = snomed_name(code)
        out_rows.append({"source": "SNOMEDCT", "id": code, "name": name or "", "url": url})

    for oid in orpha_ids:
        name, url = orpha_name(oid)
        out_rows.append({"source": "ORPHA", "id": oid, "name": name or "", "url": url})

    for ch in chebi_ids:
        name, url = chebi_name(ch)
        out_rows.append({"source": "ChEBI", "id": f"CHEBI:{ch}", "name": name or "", "url": url})

    for lid in loinc_ids:
        url = f"https://loinc.org/{lid}"
        name = loinc_map.get(lid, "") or ""
        out_rows.append({"source": "LOINC", "id": lid, "name": name, "url": url})

    os.makedirs(os.path.dirname(OUTPUT_CSV) or ".", exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["source", "id", "name", "url"])
        w.writeheader(); w.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows -> {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
