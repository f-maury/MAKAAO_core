#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ===================== IMPORT =====================
import csv
import re
import os
from collections import defaultdict
from rdflib import Graph, Namespace, RDF, RDFS, OWL, Literal, URIRef, BNode
from rdflib.namespace import XSD
from datetime import date

# ===================== CONFIG =====================

version = "1.0.0"
BASE_DIR = "../data/processed_tables/"
OUTPUT_DIR = "../data/"
makaao_core_name = "../data/makaao_core.csv"
OUTPUT_OWL_ENRICHED = f"../kg/makg-core_{version}.rdf"
OUTPUT_OWL_TBOX = f"../kg/makg-core_{version}_ontology.owl"

# Optional enrichment tables (all CSVs). Script tolerates their absence.
UMLS_HPO_CSV = "../data/enrichment_tables/umls_hpo.csv"
UMLS_ORPHANET_CSV = "../data/enrichment_tables/umls_orphanet.csv"
UMLS_SNOMED_CSV = "../data/enrichment_tables/umls_snomed.csv"
UMLS_NAMES_CSV = "../data/enrichment_tables/umls_names.csv"
ORPHANET_HPO_LINKS = "../data/enrichment_tables/orphanet_hpo_links.csv"
LOINC_UMLS = "../data/enrichment_tables/umls_loinc.csv"
LOINC_INDEX_CSV = os.path.join(BASE_DIR, "index_loinc.csv")
CODE_NAMES_CSV = "../data/enrichment_tables/code_names.csv"

# ===================== NAMESPACES =====================
MAKAAO = Namespace("http://makaao.inria.fr/kg/")
SKOS = Namespace("http://www.w3.org/2004/02/skos/core#")
PROV = Namespace("http://www.w3.org/ns/prov#")
HP = Namespace("http://purl.obolibrary.org/obo/")
BIOLINK = Namespace("https://w3id.org/biolink/vocab/")
BAO = Namespace("http://www.bioassayontology.org/bao#")
SIO = Namespace("http://semanticscience.org/resource/")
UMLS = Namespace("http://uts.nlm.nih.gov/uts/umls/concept/")

DCTERMS = Namespace("http://purl.org/dc/terms/")
DCAT = Namespace("http://www.w3.org/ns/dcat#")
ODRL = Namespace("http://www.w3.org/ns/odrl/2/")
SCHEMA = Namespace("http://schema.org/")
VOID = Namespace("http://rdfs.org/ns/void#")
FOAF = Namespace("http://xmlns.com/foaf/0.1/")
OBO = Namespace("http://www.geneontology.org/formats/oboInOwl#")

# ===================== GLOBALS =====================
relation_counter = 0
document_map = {}
GLOBAL_ACTIVITY = None

# ---------- Label de-duplication state ----------
# We maintain, per-node:
# - a cross-property set of seen texts ("any")
# - per-property sets for SKOS.prefLabel and RDFS.label
# Rule: block duplicates globally, EXCEPT allow a single cross-duplicate
#       between rdfs:label and skos:prefLabel.
_label_seen = defaultdict(
    lambda: {
        "any": set(),
        SKOS.prefLabel: set(),
        RDFS.label: set(),
    }
)


def ensure_uriref(val):  # small function to make sure a string is an URI
    if not val:
        return None  # Or raise an error
    return URIRef(val)


# normalize labels
def _norm_label_text(t: str) -> str:
    return (
        re.sub(r"\s+", " ", (t or "")).strip().casefold()
    )  # normalize whitespace and case


def add_label(g: Graph, node, prop, text: str) -> bool:
    """
    Add label triple iff:
      - text is non-empty
      - not already present on the same property
      - not already present on any property, except allow one cross-duplicate
        where the existing property is the other of {rdfs:label, skos:prefLabel}.
    """
    if not text:
        return False
    t = _norm_label_text(text)
    seen = _label_seen[node]
    # already on same property -> skip
    if t in seen[prop]:
        return False
    # cross-property allowance only for rdfs<->skos prefLabel pair
    cross_allowed = (prop == RDFS.label and t in seen[SKOS.prefLabel]) or (
        prop == SKOS.prefLabel and t in seen[RDFS.label]
    )
    if t in seen["any"] and not cross_allowed:
        return False

    g.add((node, prop, Literal(text)))
    seen[prop].add(t)
    seen["any"].add(t)
    return True  # returns True if label was added


def add_pref(g: Graph, node, text) -> bool:
    """Attach both prefLabel and rdfs:label with the cross-duplicate rule above."""
    ok1 = add_label(g, node, SKOS.prefLabel, text)
    ok2 = add_label(g, node, RDFS.label, text)
    return ok1 or ok2  # returns True if at least one label was added


# Regex to parse simple "PMID: 12345" tokens.
pmid_rx = re.compile(
    r"^\s*PMID\s*:\s*(\d+)\s*$", re.IGNORECASE
)  # matches "PMID: 12345" format


# ===================== IO / SMALL UTILS =====================
def read_csv_rows(path):
    """
    CSV → list[dict] with both original and lowercased keys.
    Handles UTF-8 BOM and trims whitespace.
    """
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        orig = [(h or "").lstrip("\ufeff").strip() for h in (rdr.fieldnames or [])]
        low = [h.lower() for h in orig]
        out = []
        for raw in rdr:
            row = {}
            for original, lower in zip(orig, low):
                v = (raw.get(original) or "").strip()
                row[original] = v  # preserve original case
                row[lower] = v  # lowercase alias
            out.append(row)
        return out  # returns list of dicts representing CSV rows


# groups dicts in lists by a key they share. returns a dict of lists of dict
def grouped(rows, key):
    g = defaultdict(list)
    for r in rows:
        k = (r.get(key) or "").strip()
        if k:
            g[k].append(r)
    return g  # returns dict {key: [row_dict, ...]}


# get the first occurrence of key/value pair from 2 columns
def read_one_map(rows, key, valkey):
    out = {}
    for r in rows:
        k = (r.get(key) or "").strip()
        v = (r.get(valkey) or "").strip()
        if k and v and k not in out:
            out[k] = v
    return out  # returns dict {key: val} of first occurrences


def to_pubmed_urls(text, *, src_file=None, row=None, col=None):
    """
    Split on '|' and ';', accept http(s) URLs or 'PMID:12345' and deduplicate.
    Rejects only newlines/tabs (hard separators).
    """
    text = (text or "").strip()
    if not text:
        return []
    for bad in ["\n", "\r", "\t"]:
        if bad in text:
            loc = []
            if src_file:
                loc.append(f"file={src_file}")
            if row is not None:
                loc.append(f"row={row}")
            if col:
                loc.append(f"col={col}")
            preview = text.replace("\n", "\\n").replace("\r", "\\r")
            if len(preview) > 160:
                preview = preview[:157] + "..."
            raise ValueError(
                f"Unexpected separator {bad!r}; only '|' and ';' are allowed. "
                + (" ".join(loc) + " " if loc else "")
                + f"value='{preview}'"
            )
    urls = []
    for slot in (s.strip() for s in text.split("|")):
        if not slot:
            continue
        for tok in (t.strip() for t in slot.split(";")):
            if not tok:
                continue
            m = pmid_rx.match(tok)
            if m:
                urls.append(f"https://pubmed.ncbi.nlm.nih.gov/{m.group(1)}")
            elif tok.lower().startswith("http://") or tok.lower().startswith(
                "https://"
            ):
                urls.append(tok)
    seen, keep = set(), []
    for u in urls:
        if u not in seen:
            seen.add(u)
            keep.append(u)
    return keep  # returns list of unique pubmed URLs


def make_valid(s):
    return (
        (s or "")
        .replace(" ", "_")
        .replace(":", "_")
        .replace("/", "_")
        .replace("|", "_")
    )  # make string safe for URI use


def hp_to_obo_uri(hp_code):
    if not hp_code:
        return None
    c = hp_code.strip().upper().replace(":", "_")
    return URIRef(f"http://purl.obolibrary.org/obo/{c}")  # returns HP OBO URI


# ===================== UniProt / ChEBI HELPERS =====================
# this removes uniprot prefix (UP:) and raises error if not present
def _strip_prefixes(s):
    s = str(s).strip()
    if s.upper().startswith("UP:"):
        return s[3:]
    raise ValueError(f"Expected 'UP:' prefix, got: {s!r}")


def canon_uniprot_id(val):
    s = (val or "").strip()
    alts = set()
    if not s:
        return "", "", []
    if s.lower().startswith("http"):
        s = s.split("/")[-1].split("?")[0].split("#")[0]
    if "|" in s:
        parts = s.split("|")
        for p in parts:
            if p:
                alts.add(p.strip())
        if len(parts) >= 2 and parts[1]:
            s = parts[1].strip()
        else:
            s = parts[-1].strip()
    s = _strip_prefixes(s).strip()
    s = s.split()[0]
    s = s.split(";")[0].strip()
    s = s.upper()
    base = s.split("-")[0]
    if s != base:
        alts.add(s)
    orig = (val or "").strip()
    if orig and orig.upper() != s:
        alts.add(orig)
    return (
        base,
        s,
        sorted(a for a in alts if a),
    )  # returns uniprot base, uniprot full, [alt ids]


_chebi_num_rx = re.compile(r"(\d+)$")


def canon_chebi_id(val):
    s = (val or "").strip()
    if not s:
        raise ValueError("Empty ChEBI value.")
    if s.lower().startswith("http"):
        s = s.split("/")[-1].split("?")[0].split("#")[0]
    if not s.upper().startswith("CHEBI:"):
        raise ValueError(
            f"Expected 'CHEBI:' prefix (or URL ending with it), got: {val!r}"
        )
    m = _chebi_num_rx.search(s)
    if not m:
        raise ValueError(f"ChEBI ID missing numeric part: {val!r}")
    num = m.group(1)
    return f"CHEBI:{num}", f"CHEBI_{num}"  # returns (CHEBI:nnnnn, CHEBI_nnnnn)


# ===================== code_names.csv READERS =====================


# build uniprot lookup table from reading csv
def read_code_names_uniprot(path):
    id2name, id2url = {}, {}
    rows = read_csv_rows(path)
    for r in rows:
        src = (r.get("source") or "").strip().lower()
        if src != "uniprot":
            continue
        uid = (r.get("id") or "").strip().upper()
        nm = (r.get("name") or "").strip()
        url = (r.get("url") or "").strip()
        if not uid:
            continue
        base = uid.split("-")[0]
        if nm:
            if base not in id2name:
                id2name[base] = nm
            if uid not in id2name:
                id2name[uid] = nm
        if url:
            if base not in id2url:
                id2url[base] = url
            if uid not in id2url:
                id2url[uid] = url
    return id2name, id2url  # 2 dicts: UniProt ID -> name, UniProt ID -> url


def read_code_names_umls(path):
    id2name, id2url = {}, {}
    rows = read_csv_rows(path)
    for r in rows:
        if (r.get("source") or "").strip().lower() != "umls":
            continue
        cui = (r.get("id") or "").strip().upper().replace("CUI:", "")
        nm = (r.get("name") or "").strip()
        url = (r.get("url") or "").strip()
        if cui and nm and cui not in id2name:
            id2name[cui] = nm
        if cui and url and cui not in id2url:
            id2url[cui] = url

    return id2name, id2url  # 2 dicts: CUI -> name, CUI -> url

def read_code_names_orpha(path):
    id2name, id2url = {}, {}
    rows = read_csv_rows(path)
    for r in rows:
        src = (r.get("source") or "").strip().lower()
        if "orpha" not in src and src != "che":
            continue
        raw = (r.get("id") or "").strip()
        nm = (r.get("name") or "").strip()
        url = (r.get("url") or "").strip()
        code_colon = raw.strip()
        if nm and code_colon and code_colon not in id2name:
            id2name[code_colon] = nm
        if url and code_colon and code_colon not in id2url:
            id2url[code_colon] = url

    return id2name, id2url


def read_code_names_chebi(path):
    id2name, id2url = {}, {}
    rows = read_csv_rows(path)
    for r in rows:
        src = (r.get("source") or "").strip().lower()
        if "chebi" not in src and src != "che":
            continue
        raw = (r.get("id") or "").strip()
        nm = (r.get("name") or "").strip()
        url = (r.get("url") or "").strip()
        code_colon, _ = canon_chebi_id(raw)
        if nm and code_colon and code_colon not in id2name:
            id2name[code_colon] = nm
        if url and code_colon and code_colon not in id2url:
            id2url[code_colon] = url

    return id2name, id2url  # 2 dicts: CHeBI code -> name, CHeBI code -> url


# ===================== REIFICATION =====================
def add_reified_relation(g, subj, pred, obj, prov_str):
    """
    Reify (subj,pred,obj) once; if triple already reified, only append new prov:wasDerivedFrom.
    """
    subj = ensure_uriref(subj)
    global relation_counter, document_map, GLOBAL_ACTIVITY  # use global relation counter and document map
    if not prov_str or str(prov_str).strip().lower() in {"", "nan"}:
        return
    prov = str(prov_str).strip()

    # if that document is already known, reuse it
    if prov in document_map:
        doc = document_map[prov]
    else:
        doc = MAKAAO[
            f"document_{len(document_map) + 1}"
        ]  # otherwise we create a new instance of Document, and we construct its URI based on the number of documents already known
        document_map[prov] = doc  # we add that new document to the document map
        g.add((doc, RDF.type, MAKAAO.Document))
        g.add((doc, OBO.xref, URIRef(prov)))  # the label of the docment is its provenance URI

    # check if the object of the relations is an individual or a literal, so we characterize it correctly
    obj_node = obj if isinstance(obj, (URIRef, BNode, Literal)) else Literal(obj)

    # Try to find existing Relation with same subj/pred/obj
    existing_rel = None
    for rel in g.subjects(RDF.subject, subj):
        if (rel, RDF.predicate, pred) in g and (
            rel,
            RDF.object,
            obj_node,
        ) in g:  # if the relation already exist, there is nothing to do
            existing_rel = rel
            break

    if existing_rel is None:
        rel = MAKAAO[f"r{relation_counter}"]
        relation_counter += 1  # if that relation does not exist yet, we create it, and increment the relation counter
        g.add((rel, RDF.type, MAKAAO.Relation))
        g.add((rel, RDF.type, RDF.Statement))
        g.add(
            (rel, RDF.subject, subj)
        )  # we add some triples about the new reified relation
        g.add((rel, RDF.predicate, pred))
        g.add((rel, RDF.object, obj_node))
        if GLOBAL_ACTIVITY is not None:
            g.add((rel, PROV.wasGeneratedBy, GLOBAL_ACTIVITY))
    else:
        rel = existing_rel

    g.add(
        (rel, PROV.wasDerivedFrom, doc)
    )  # we add the provenance information to the reified relation


# ===================== ENRICHMENT HELPERS =====================
def umls_names_index(path):  # read UMLS names CSV table into dict CUI -> [names]
    rows = read_csv_rows(path)
    names = defaultdict(list)
    for r in rows:
        cui = (r.get("CUI") or "").strip()
        s = (r.get("STR") or "").strip()
        if cui and s:
            names[cui].append(s)
    return names


def umls_code_group(
    path, key="CUI"
):  # read UMLS code names table into two dicts grouped by key and by CODE
    def canon_row(row: dict) -> dict:
        out = {}
        for k, v in row.items():
            ku = k.upper()
            if ku not in out or k.isupper():
                out = {**out, ku: v.strip() if isinstance(v, str) else v}
        return out

    rows = [canon_row(r) for r in read_csv_rows(path)]
    k = key.upper()
    by_key = grouped(rows, k)
    by_code = grouped(rows, "CODE")
    return by_key, by_code


# ===================== GRAPH BUILDERS =====================
def init_graph():  # start an empty knowledge graph and add a few basic things to it
    g = Graph()
    # Prefixes
    g.bind("mak", MAKAAO)
    g.bind("rdfs", RDFS)
    g.bind("owl", OWL)
    g.bind("skos", SKOS)
    g.bind("prov", PROV)
    g.bind("hp", HP)
    g.bind("biolink", BIOLINK)
    g.bind("bao", BAO)
    g.bind("sio", SIO)
    g.bind("dcterms", DCTERMS)
    g.bind("dcat", DCAT)
    g.bind("odrl", ODRL)
    g.bind("schema", SCHEMA)
    g.bind("void", VOID)
    g.bind("foaf", FOAF)
    g.bind("obo", OBO)

    # Ontology header
    onto = URIRef("http://makaao.inria.fr/kg/")
    g.add((onto, RDF.type, OWL.Ontology))
    add_pref(g, onto, "MAKAAO knowledge graph")
    g.add((onto, OWL.versionIRI, URIRef(f"http://makaao.inria.fr/kg/{version}")))
    g.add((onto, OWL.versionInfo, Literal(version)))

    # we manually define the main classes
    g.add((MAKAAO.Relation, RDF.type, OWL.Class))
    g.add((MAKAAO.Relation, RDFS.subClassOf, PROV.Entity))  # sublass of PROV Entity, used for reified relations that carry provenance information
    g.add((MAKAAO.Relation, RDFS.label, Literal("Relation")))
    g.add((MAKAAO.Document, RDFS.label, Literal("Document")))
    g.add((MAKAAO.Document, RDF.type, OWL.Class))
    g.add((MAKAAO.Document, RDFS.subClassOf, PROV.Entity))  # sublass of PROV Entity, used for reified relations that carry provenance information
    g.add((MAKAAO.Autoantibody, RDF.type, OWL.Class))
    g.add((MAKAAO.Autoantibody, RDFS.label, Literal("Autoantibody")))
    g.add((MAKAAO.AutoimmuneDisease, RDF.type, OWL.Class))
    g.add((MAKAAO.AutoimmuneDisease, RDFS.label, Literal("Autoimmune disease")))
    # g.add((MAKAAO.Autoantibody, SKOS.closeMatch, URIRef("http://snomed.info/id/30621004"))) # we define a closeMatch between our Autoantibody class and the corresponding one in SNOMED
    g.add((MAKAAO.Target, RDF.type, OWL.Class))
    g.add((MAKAAO.Target, RDFS.label, Literal("Target")))

    # Object properties we will use
    g.add((BAO.BAO_0000211, RDF.type, OWL.ObjectProperty))  # has_target
    g.add((BAO.BAO_0000211, RDFS.label, Literal("has target")))
    g.add((BAO.BAO_0000598, RDF.type, OWL.ObjectProperty))  # is_target_of
    g.add((BAO.BAO_0000598, RDFS.label, Literal("is target of")))
    g.add((SIO["SIO_001279"], RDF.type, OWL.ObjectProperty))  # has_phenotype
    g.add((SIO["SIO_001279"], RDFS.label, Literal("has phenotype")))
    g.add((SIO["SIO_001280"], RDF.type, OWL.ObjectProperty))  # is_phenotype_of
    g.add((SIO["SIO_001280"], RDFS.label, Literal("is phenotype of")))
    g.add((SIO["SIO_001403"], RDF.type, OWL.ObjectProperty))  # is_associated_with
    g.add((SIO["SIO_001403"], RDFS.label, Literal("is associated with")))
    g.add((BIOLINK.has_biomarker, RDF.type, OWL.ObjectProperty))  # has_biomarker
    g.add((BIOLINK.has_biomarker, RDFS.label, Literal("has biomarker")))
    g.add((BIOLINK.biomarker_for, RDF.type, OWL.ObjectProperty))  # biomarker_for
    g.add((BIOLINK.biomarker_for, RDFS.label, Literal("biomarker for")))

    # Global PROV activity
    global GLOBAL_ACTIVITY
    act = MAKAAO["activity_makaao_core"]
    g.add((act, RDF.type, PROV.Activity))
    add_pref(g, act, "MAKAAO Core processing activity")
    agent = MAKAAO["agent_makaao_core_processing_py"]  # prov:SoftwareAgent
    g.add((agent, RDF.type, PROV.SoftwareAgent))
    add_pref(g, agent, "03_build_kg_from_tables.py")
    g.add(
        (act, PROV.wasAssociatedWith, agent)
    )  # the makaao_processing_script.py (the script you are reading now, to be renamed) is associated with the MAKAAO processing acivity
    GLOBAL_ACTIVITY = act
    return g  # returned initialized KG


# read data from the processed tables, and return it structured in acomplex dictionary
def load_processed_tables(base_dir):
    parents_path = os.path.join(base_dir, "index_parent_index.csv")
    name_en_path = os.path.join(base_dir, "index_name_en.csv")
    syn_en_path = os.path.join(base_dir, "index_syn_source_en.csv")
    syn_fr_path = os.path.join(base_dir, "index_syn_source_fr.csv")
    hpo_path = os.path.join(base_dir, "index_hpo_id.csv")
    cui_path = os.path.join(base_dir, "index_cui_source.csv")
    uni_path = os.path.join(base_dir, "index_uniprot_source.csv")
    che_path = os.path.join(base_dir, "index_chebi_source.csv")
    dis_path = os.path.join(base_dir, "index_disease_source.csv")

    parents = read_csv_rows(parents_path)
    name_en = read_csv_rows(name_en_path)
    syn_en = read_csv_rows(syn_en_path)
    syn_fr = read_csv_rows(syn_fr_path)
    hpo_id = read_csv_rows(hpo_path)
    cui_src = read_csv_rows(cui_path)
    uni_src = read_csv_rows(uni_path)
    che_src = read_csv_rows(che_path)
    dis_src = read_csv_rows(dis_path)

    name_en_map = read_one_map(name_en, "index", "name_en")  # dict: {aab_id: name_en}
    KEEP = set(
        name_en_map.keys()
    )  # indices that we will not filter out (non-empty rows)

    syn_en_map = defaultdict(list)  # dict: {aab_id: [(syn_en, source), ...]}
    for r in syn_en:
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        se = (r.get("syns_en") or "").strip()
        src = (r.get("syns_en_source") or "").strip()
        if idx and se:
            syn_en_map[idx].append((se, src))

    syn_fr_map = defaultdict(list)  # dict: {aab_id: [(syn_fr, source), ...]}
    for r in syn_fr:
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        sf = (r.get("syns_fr") or "").strip()
        src = (r.get("syns_fr_source") or "").strip()
        if idx and sf:
            syn_fr_map[idx].append((sf, src))

    hpo_by_idx = defaultdict(list)  # dict: {aab_id: [hpo_id, ...]}
    for r in hpo_id:
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        hp = (r.get("hpo_id") or "").strip()
        if idx and hp and hp not in hpo_by_idx[idx]:
            hpo_by_idx[idx].append(hp)

    parent_map = defaultdict(set)  # dict: {aab_id: set(parent_aab_id, ...)}
    for r in parents:
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        p = (r.get("parent_index") or "").strip()
        if idx and p and p != idx:
            parent_map[idx].add(p)

    cui_by_idx = defaultdict(list)  # dict: {aab_id: [(cui, [pmid_url, ...]), ...]}
    for i, r in enumerate(cui_src, start=2):
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        cui = (r.get("umls_target_cui") or "").strip().upper().replace("CUI:", "")
        pm = to_pubmed_urls(
            r.get("umls_pmids"), src_file=cui_path, row=i, col="umls_pmids"
        )
        if idx and cui:
            cui_by_idx[idx].append((cui, pm))

    uni_by_idx = defaultdict(
        list
    )  # dict: {aab_id: [(uniprot_id, [pmid_url, ...]), ...]}
    for i, r in enumerate(uni_src, start=2):
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        uid = (r.get("uniprot_target_id") or "").strip()
        pm = to_pubmed_urls(
            r.get("uniprot_pmids"), src_file=uni_path, row=i, col="uniprot_pmids"
        )
        if idx and uid:
            uni_by_idx[idx].append((uid, pm))

    che_by_idx = defaultdict(list)  # dict: {aab_id: [(chebi_id, [pmid_url, ...]), ...]}
    for i, r in enumerate(che_src, start=2):
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        cid = (r.get("chebi_target_id") or "").strip()
        pm = to_pubmed_urls(
            r.get("chebi_pmids"), src_file=che_path, row=i, col="chebi_pmids"
        )
        if idx and cid:
            che_by_idx[idx].append((cid, pm))

    dis_by_idx = defaultdict(
        list
    )  # dict: {aab_id: [(disease_id, [pmid_url, ...]), ...]}
    for i, r in enumerate(dis_src, start=2):
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        did = (r.get("related_diseases_id") or "").strip()
        pm = to_pubmed_urls(
            r.get("diseases_pmids"), src_file=dis_path, row=i, col="diseases_pmids"
        )
        if idx and did:
            dis_by_idx[idx].append((did, pm))

    return {
        "name_en": name_en_map,
        "syn_en": syn_en_map,
        "syn_fr": syn_fr_map,
        "hpo_list": hpo_by_idx,
        "parents": parent_map,
        "cui": cui_by_idx,
        "uniprot": uni_by_idx,
        "chebi": che_by_idx,
        "diseases": dis_by_idx,
        "indices": KEEP,
    }


def build_core(
    g,
    data,
    umls_hpo_rows,
    umls_names,
    up_names,
    up_urls,
    umls_cn_names,
    chebi_cn_names,
    chebi_cn_urls,
):
    # Target class specializations
    uniprot_protein = URIRef(
        "http://purl.uniprot.org/core/Protein"
    )  # define the root classes for uniprot proteins and chebi molecules
    g.add((uniprot_protein, RDFS.subClassOf, MAKAAO.Target))
    chebi_root = URIRef("http://purl.obolibrary.org/obo/CHEBI_23367")
    g.add((chebi_root, RDF.type, OWL.Class))
    g.add((chebi_root, RDFS.subClassOf, MAKAAO.Target))

    # Declare makaao_core CSV as a Document used by the Activity
    csv_doc = MAKAAO["makaao_core.csv"]
    g.add((csv_doc, RDF.type, MAKAAO.Document))
    g.add((csv_doc, RDFS.seeAlso, URIRef("https://makaao.inria.fr/data/makaao_core.csv")))
    add_pref(g, csv_doc, os.path.basename(makaao_core_name))
    g.add((GLOBAL_ACTIVITY, PROV.used, csv_doc))

    umls_hpo_by_code = (
        grouped(umls_hpo_rows, "CODE") if umls_hpo_rows else {}
    )  # create index of umls_hpo by code
    primary_labels, aab_class_uri = {}, {}
    pos_uris_by_idx = {}

    # Build classes/instances
    for idx in sorted(data["indices"], key=lambda x: int(x)):
        cls = (
            MAKAAO.Autoantibody if idx == "18" else MAKAAO[f"aab_{idx}"]
        )  # 18 is autoimmune antibody: the index of the root class for all autoantibodies
        if idx != "18":
            g.add((cls, RDF.type, OWL.Class))
        aab_class_uri[idx] = cls

        pref = data["name_en"].get(idx, f"aab_{idx}")  # add pref label to AAb classes
        primary_labels[idx] = pref
        add_pref(g, cls, pref)

        for syn, src in data["syn_en"].get(idx, []):
            if (
                add_label(g, cls, RDFS.label, syn) and src
            ):  # add english as syns as rdfs:labels, and add corresponding reified relation with provenance (if there is provenance)
                add_reified_relation(g, cls, RDFS.label, Literal(syn), src)
        for syn, src in data["syn_fr"].get(idx, []):  # same thing with french syns
            if add_label(g, cls, RDFS.label, syn) and src:
                add_reified_relation(g, cls, RDFS.label, Literal(syn), src)

        # taxonomy
        parents = data["parents"].get(idx) or []
        if idx != "18":  # root AAb has no parent
            valid = [p for p in parents if p != idx and p in data["indices"]]
            if valid:
                for p in valid:
                    parent_uri = (
                        MAKAAO.Autoantibody if p == "18" else MAKAAO[f"aab_{p}"]
                    )  # we reconstruct the URI of the parent class from its ID
                    g.add(
                        (cls, RDFS.subClassOf, parent_uri)
                    )  # we add the subsumption triple
            else:
                g.add(
                    (cls, RDFS.subClassOf, MAKAAO.Autoantibody)
                )  # if no valid parent, we link to root AAb

        # instance
        inst = MAKAAO[f"aab_{idx}_instance"]
        g.add(
            (inst, RDF.type, cls)
        )  # we create an instance of the current AAb class with its preferred label
        add_pref(g, inst, primary_labels[idx])

        # positivity (HPO-backed if present)
        pos_list = []
        for hp_code in data["hpo_list"].get(
            idx, []
        ):  # we try to get the associated HPO positivity(ies) for that AAb, if it exits
            pos_uri = hp_to_obo_uri(hp_code)
            if not pos_uri:
                continue
            code_norm = hp_code.replace("_", ":").upper()
            pt_label = None
            for row in umls_hpo_by_code.get(
                code_norm, []
            ):  # we try to get the main enlgish name for the HPO code from our UMLS-HPO table
                s = (row.get("STR") or "").strip()
                tty = (row.get("TTY") or "").strip().upper()
                if not s:
                    continue
                if tty == "PT" and pt_label is None:
                    pt_label = s
                    add_pref(g, pos_uri, s)
                add_label(
                    g, pos_uri, RDFS.label, s
                )  # if we find it, we add it as a skos:prefLabel and rdfs:label
            if pt_label is None:
                add_pref(
                    g, pos_uri, code_norm
                )  # if we don't find a preferred name, we just add the code as prefLabel
            pos_list.append((pos_uri, pt_label or code_norm))

        if not pos_list:  # if there is no HPO positivity, we create a generic positivity class for that AAb
            pos_uri = MAKAAO[f"positivity_{idx}"]
            add_pref(g, pos_uri, primary_labels[idx] + " positivity")
            g.add((pos_uri, RDF.type, OWL.Class))
            pos_list.append((pos_uri, primary_labels[idx] + " positivity"))

        pos_uris_by_idx[idx] = [
            u for (u, _) in pos_list
        ]  # positivity URIS associated to each aab_id

        # positivity instances and biolink
        n = len(pos_list)
        for j, (pos_uri, pos_label) in enumerate(
            pos_list, start=1
        ):  # for each positivity associated to current AAb, we create an instance and link it with biolink relations
            suffix = "" if n == 1 else f"_{j}"
            pos_inst = MAKAAO[f"positivity_{idx}{suffix}_instance"]
            g.add((pos_inst, RDF.type, pos_uri))
            add_pref(g, pos_inst, pos_label)
            g.add((pos_inst, BIOLINK.has_biomarker, inst))
            g.add((inst, BIOLINK.biomarker_for, pos_inst))

    # mirror taxonomy to positivity
    for idx, parents in data["parents"].items():
        child_list = pos_uris_by_idx.get(
            idx, [MAKAAO[f"positivity_{idx}"]]
        )  # we get the positivity URIS for the current AAb
        for p in parents:
            if p not in data["indices"]:
                continue
            parent_list = pos_uris_by_idx.get(
                p, [MAKAAO[f"positivity_{p}"]]
            )  # for eahc parent of the current aab, we get the corresponding positivity URI
            for c_uri in child_list:
                for p_uri in parent_list:
                    g.add(
                        (c_uri, RDFS.subClassOf, p_uri)
                    )  # we make the subsumption link between the child and parent positivity URIs
                    

    # UMLS targets
    umls_local_names = umls_cn_names or {}  # umls_names from code_names table
    for idx, items in data["cui"].items():
        inst_aab = MAKAAO[
            f"aab_{idx}_instance"
        ]  # for each aab_id, we get the name of the corresponding instance
        for cui, pmids in items:
            cui_key = cui.upper().replace(
                "CUI:", ""
            )  # for each CUI associated to this aab_id
            cui_uri = URIRef(
                MAKAAO + "CUI_" + cui_key + "_instance"
            )  # we create the URI of the CUI instanceS
            g.add((cui_uri, RDF.type, MAKAAO.Target))
            g.add((cui_uri, OBO.xref, UMLS[cui_key])) # we create a instance of the Target class with the relevant CUI as URI
            all_umls_labels = (
                umls_names.get(cui_key) or []
            )  # we read from umls_names table all the names associated to that CUI
            preferred = (
                umls_local_names.get(cui_key)
                or (all_umls_labels[0] if all_umls_labels else None)
            )  # we try tp get preferred name from code_names table, or else from umls_name table
            if preferred:
                add_pref(g, cui_uri, preferred)
            else:
                add_pref(
                    g, cui_uri, cui_key
                )  # if we find a pref name, we add it as pref label, else we take the CUI as label
            for n in all_umls_labels:
                add_label(
                    g, cui_uri, RDFS.label, n
                )  # we add all the other names as rdfs:labels
            g.add(
                (inst_aab, BAO.BAO_0000211, cui_uri)
            )  # we add relations between target instance and AAb instance
            g.add((cui_uri, BAO.BAO_0000598, inst_aab))
            for p in pmids or [""]:
                add_reified_relation(
                    g, inst_aab, BAO.BAO_0000211, cui_uri, p
                )  # for each source we have, we add a reified relation to carry provenance
                add_reified_relation(g, cui_uri, BAO.BAO_0000598, inst_aab, p)

    # UniProt targets
    up_names_added, up_total = 0, 0
    for idx, items in data[
        "uniprot"
    ].items():  # for each aab_id, we get the associated uniprot ids
        inst_aab = MAKAAO[f"aab_{idx}_instance"]
        for upid_in, pmids in items:
            up_total += 1
            base, norm, _ = canon_uniprot_id(
                upid_in
            )  # we canonicalize the uniprot id read from data
            if not base:
                continue
            prot_ind = MAKAAO[
                f"{make_valid('UP_' + base + '_instance')}"
            ]  # we create the URI of the uniprot protein instance
            g.add(
                (prot_ind, RDF.type, URIRef("http://purl.uniprot.org/core/Protein"))
            )
            g.add(
                (prot_ind, RDF.type, MAKAAO.Target)
            )  # we declare it as an instance of uniprot protein class
            up_name = up_names.get(base) or up_names.get(
                norm
            )  # we try to get the name of the uniprot protein from code_names table
            if up_name:
                add_pref(g, prot_ind, up_name)
                add_label(
                    g, prot_ind, RDFS.label, base
                )  # we add prefLabel and rdfs:label if we find a name
                up_names_added += 1
            else:
                add_pref(g, prot_ind, base)
            up_url = (
                up_urls.get(base) or up_urls.get(norm)
            )  # else, we try to get the url of the uniprot protein from code_names table
            if up_url:
                g.add(
                    (prot_ind, OBO.xref, URIRef(up_url))
                )  # we add xref relation if we find a url (the real uniprot URI)
            g.add((inst_aab, BAO.BAO_0000211, prot_ind))
            g.add(
                (prot_ind, BAO.BAO_0000598, inst_aab)
            )  # we add relations between target instance and AAb instance
            for p in pmids or [""]:
                add_reified_relation(
                    g, inst_aab, BAO.BAO_0000211, prot_ind, p
                )  # for each source we have, we add a reified relation to carry provenance
                add_reified_relation(g, prot_ind, BAO.BAO_0000598, inst_aab, p)
    print(f"UniProt targets: named={up_names_added}/{up_total}")

    # ChEBI targets
    for idx, items in data["chebi"].items():
        inst_aab = MAKAAO[f"aab_{idx}_instance"]
        for chebi_raw, pmids in items:  # for each chebi id associated to that aab_id
            code_colon, code_obo = canon_chebi_id(chebi_raw)
            chebi_cls = URIRef(
                "http://purl.obolibrary.org/obo/" + code_obo
            )  # we reconstruct the chebi class URI from the chebi id read from data
            chebi_ind = MAKAAO[
                code_obo + "_instance"
            ]  # we create the URI of the chebi instance
            g.add(
                (chebi_cls, RDF.type, OWL.Class)
            )  # we declare the chebi class as an OWL class, subclass of chebi molecule
            g.add(
                (
                    chebi_cls,
                    RDFS.subClassOf,
                    URIRef("http://purl.obolibrary.org/obo/CHEBI_23367"),
                )
            )
            g.add(
                (chebi_ind, RDF.type, chebi_cls)
            )  # we create an instance of the chebi class with the relevant URI
            g.add(
                (chebi_ind, RDF.type, MAKAAO.Target)
            )
            chebi_name = (chebi_cn_names or {}).get(code_colon)
            if chebi_name:  # if we find a name for that chebi id from code_names table, we add it as prefLabel and rdfs:label
                add_pref(g, chebi_ind, chebi_name)
                add_pref(g, chebi_cls, chebi_name)
                add_label(g, chebi_ind, RDFS.label, code_colon)
            else:  # else we just add the chebi id as prefLabel and rdfs:label
                add_pref(g, chebi_ind, code_colon)
                add_pref(g, chebi_cls, code_colon)
            chebi_url = (chebi_cn_urls or {}).get(
                code_colon
            )  # we try to get the chebi url from code_names table
            if chebi_url:
                g.add(
                    (chebi_ind, RDFS.seeAlso, URIRef(chebi_url))
                )  # we add seeAlso relation if we find a url (the real chebi URI)
            g.add(
                (inst_aab, BAO.BAO_0000211, chebi_ind)
            )  # we add relations between target instance and AAb instance
            g.add((chebi_ind, BAO.BAO_0000598, inst_aab))
            for p in pmids or [""]:
                add_reified_relation(
                    g, inst_aab, BAO.BAO_0000211, chebi_ind, p
                )  # for each source we have, we add a reified relation to carry provenance
                add_reified_relation(g, chebi_ind, BAO.BAO_0000598, inst_aab, p)


def process_diseases(
    g,
    data,
    umls_orphanet_by_code,
    umls_names,
    orpha_hpo_links,
    umls_cn_names,
):
    # HPO phenotypic abnormality root (HP:0000118)
    pheno_root = hp_to_obo_uri("HP:0000118")

    # HPO codes used as positivity phenotypes for AAbs
    # (these correspond to the HP:0030057 subtree in your KG and should
    # not be additionally forced under HP:0000118 here).
    positivity_codes = set()
    for _idx, codes in (data.get("hpo_list") or {}).items():
        for hp_code in codes or []:
            code_norm = (hp_code or "").strip().upper().replace("_", ":")
            if code_norm:
                positivity_codes.add(code_norm)

    orpha_cn_names, orpha_cn_urls = read_code_names_orpha(CODE_NAMES_CSV)

    for idx, items in data["diseases"].items():
        aab_inst = MAKAAO[f"aab_{idx}_instance"]
        for (
            code_raw,
            prov,
        ) in items:  # loop for every disease code associated to that aab_id in the data
            code_upper = (code_raw or "").strip().upper()
            if not code_upper:
                continue
            inst = None

            # ORPHANET
            if code_upper.startswith(
                "ORPHA:"
            ):  # if current disease code is an Orpha code:
                orpha_num = code_upper.split(":", 1)[1].strip()
                d_cls = URIRef(f"http://www.orpha.net/ORDO/Orphanet_{orpha_num}")
                g.add((d_cls, RDF.type, OWL.Class))
                g.add(
                    (
                        d_cls,
                        RDFS.subClassOf,
                        URIRef("http://www.orpha.net/ORDO/Orphanet_C001"),
                    )
                )  # reconstruct Orphanet URI from Orpha code read in data, and add relevant triples to the graph
                g.add((d_cls, RDFS.subClassOf, MAKAAO.AutoimmuneDisease))
                for r in umls_orphanet_by_code.get(
                    orpha_num, []
                ):  # try to get the main english term from UMLS-Orphanet mapping to label the Orphanet code
                    s = (r.get("STR") or "").strip()
                    tty = (r.get("TTY") or "").strip().upper()
                    if not s:
                        continue
                    if tty == "PT":
                        add_pref(g, d_cls, s)
                    add_label(g, d_cls, RDFS.label, s)
                inst = MAKAAO[
                    f"orpha_{orpha_num}_instance"
                ]  # also add instance of Orphanet class
                g.add((inst, RDF.type, d_cls))
                try:
                    g.add((inst, RDFS.label, Literal(orpha_cn_names[orpha_num])))
                except Exception as e:
                    print("Error adding Orpha code name:", e)


                # Use the ORDO URI string (same as keys in orpha_hpo_links)
                for link in orpha_hpo_links.get(str(d_cls), []):
                    # HPO identifier (may be HP:..., HP_..., or full URI)
                    hpo_id_raw = (link.get("HPOId") or link.get("hpoid") or "").strip()
                    if not hpo_id_raw:
                        continue

                    # Normalize HPO code and build URI
                    if hpo_id_raw.startswith("http"):
                        last = hpo_id_raw.rsplit("/", 1)[-1]
                        if last.upper().startswith("HP_"):
                            code_norm = "HP:" + last[3:]
                        else:
                            code_norm = last.replace("_", ":").upper()
                        pos_uri = URIRef(hpo_id_raw)
                    else:
                        code_norm = hpo_id_raw.strip().upper().replace("_", ":")
                        pos_uri = hp_to_obo_uri(code_norm)

                    if pos_uri is None:
                        continue

                    # Add subclass of HP:0000118 for "regular" phenotypic abnormalities,
                    # but do NOT add it for positivity phenotypes already used for AAbs.
                    if pheno_root is not None and code_norm not in positivity_codes:
                        g.add((pos_uri, RDFS.subClassOf, pheno_root))

                    # Label from Orphadata (HPOTerm / hpoterm)
                    term = (link.get("HPOTerm") or link.get("hpoterm") or "").strip()
                    if term:
                        add_pref(g, pos_uri, term)

                    # Instance and phenotype links
                    ########################################
                    # code_norm : normalized code HP:xxx of phenotype linked to Orpha disease
                    # before adding an instance of that HP_xxx class linked to the orpha disease, 
                    # we need to check that HP_xxx class is not already used as positivity phenotype for an AAb
                    
                    target_code = code_norm.strip().upper()
                    
                    # Search for the key where the target_code exists inside the list of values (v_list)
                    key = next(
                        (
                            k for k, v_list in data["hpo_list"].items() 
                            if target_code in [str(x).strip().upper() for x in v_list]
                        ), 
                        None
                    )

                    if key is not None:
                        pos_inst = next((s for s in g.subjects(BIOLINK.has_biomarker, MAKAAO[f"aab_{key}_instance"]) if (s, RDF.type, pos_uri) in g), URIRef(MAKAAO[f"positivity_{key}_instance"]),)
                    elif key is None:
                        pos_inst = MAKAAO[f"{make_valid(str(pos_uri).split('/')[-1])}_instance"]
                    if key is not None:
                        # print(key)
                        pass
                    g.add((pos_inst, RDF.type, pos_uri))
                    g.add((inst, SIO["SIO_001279"], pos_inst))  # has_phenotype
                    g.add((pos_inst, SIO["SIO_001280"], inst))  # is_phenotype_of
                    class_label = g.value(pos_uri, RDFS.label) or g.value(pos_uri, SKOS.prefLabel)
                    if class_label:
                        g.add((pos_inst, RDFS.label, class_label))
                    add_reified_relation(
                        g,
                        inst,
                        SIO["SIO_001279"],
                        pos_inst,
                        "https://www.orphadata.com/data/xml/en_product4.xml",
                    )
                    add_reified_relation(
                        g,
                        pos_inst,
                        SIO["SIO_001280"],
                        inst,
                        "https://www.orphadata.com/data/xml/en_product4.xml",
                    )

            # UMLS CUI
            elif code_upper.startswith("C"):  # if current disease code is a CUI:
                code_norm = code_upper.replace("CUI:", "").strip()
                cui_uri = URIRef(UMLS + code_norm)  # reconstruct UMLS URI from code
                makaao_cui_uri = URIRef(MAKAAO + "CUI_" + code_norm + "_instance")
                inst = makaao_cui_uri
                preferred = (umls_cn_names or {}).get(
                    code_norm
                )  # get preferred english label from code_names table, if it's there
                all_umls_labels = (
                    umls_names.get(code_norm) or []
                )  # as 2nd option, get preferred eglish name from umls_names table, if it's there
                if preferred:
                    add_pref(g, inst, preferred)
                elif all_umls_labels:
                    add_pref(g, inst, all_umls_labels[0])
                else:
                    add_pref(
                        g, inst, code_norm
                    )  # otherwise we add the code as a prefLabel, if we didn't get an english name
                for n in (
                    all_umls_labels
                ):  # we add all the labels we found as regular labels
                    add_label(g, inst, RDFS.label, n)
                g.add((inst, OBO.xref, cui_uri))

            # last case, if code is not a properly formatted orpha or CUI code, we just create a generic instance with the code as label
            else:
                code_norm = code_upper.replace("CUI:", "").strip()
                inst = URIRef(MAKAAO + code_norm + "_instance")
                add_pref(g, inst, code_upper)

            if inst is None:
                continue
            g.add(
                (inst, RDF.type, MAKAAO.AutoimmuneDisease)
            )  # if no instance has been created yet, we create one.
            g.add(
                (aab_inst, SIO["SIO_001403"], inst)
            )  # we link that disease instance to the current aab instance
            g.add((inst, SIO["SIO_001403"], aab_inst))
            for p in prov or [""]:
                add_reified_relation(
                    g, aab_inst, SIO["SIO_001403"], inst, p
                )  # for each provenance entry, we add a reified relation to the graph, for the current instance


def process_loinc_mappings(g, loinc_umls_csv, loinc_index_csv, keep_indices):
    loinc_rows = read_csv_rows(
        loinc_umls_csv
    )  # that tables contains all MRCONOS rows associated to a LOINC code
    map_rows = read_csv_rows(
        loinc_index_csv
    )  # one of our processed dat table where we mapped LOINC codes to aab_id
    if not map_rows:
        return
    loinc2cui = {}
    for r in loinc_rows or []:
        code = (r.get("CODE") or "").strip()
        cui = (r.get("CUI") or "").strip()
        if code and cui and code not in loinc2cui:
            loinc2cui[code] = cui  # build a dic mapping LOINC codes to UMLS CUIs

    seen = set()
    for r in map_rows:
        idx = (r.get("aab_id") or "").strip()
        if not idx or idx not in keep_indices:
            continue
        code = (r.get("loinc_id") or "").strip()
        if not code:
            continue
        code = (
            code.replace("LOINC:", "").replace("loinc:", "").strip()
        )  # reconstruct LOINC URI from LOINC ID read in the  data
        loinc_uri = f"https://loinc.org/{code}"
        key = (idx, loinc_uri)
        if key in seen:
            continue
        seen.add(key)
        aab_cls = MAKAAO.Autoantibody if idx == "18" else MAKAAO[f"aab_{idx}"]
        L = URIRef(loinc_uri)
        g.add((aab_cls, OBO.xref, L))
        # g.add((L, SKOS.closeMatch, aab_cls)) # link the LOINC URI to the URI of current AAb class
        cui = loinc2cui.get(code)
        if cui:
            C = URIRef(
                UMLS + cui
            )  # URI of a UMLS concept corresopnding to a LOINC code
            g.add((aab_cls, OBO.xref, C))
            # g.add((C, SKOS.closeMatch, aab_cls)) # also look if that LOINC code has a UMLS CUI associated to it, annd if yes, also link it to the AAb


# def process_snomed_mappings : that was a function to add matches from a files with matches automatically found between SNOMED and our AAb, but it is now removed


# ===================== FAIR / DCAT METADATA =====================
# function to add some hard-coded triples to the KG
def append_fair_metadata(kg: Graph):
    ONT = URIRef("http://makaao.inria.fr/kg/")
    kg.add((ONT, RDF.type, OWL.Ontology))
    kg.add((ONT, RDF.type, DCAT.Dataset))
    kg.add((ONT, RDF.type, VOID.Dataset))

    kg.add((ONT, DCTERMS.identifier, ONT))
    kg.add((ONT, DCTERMS.title, Literal("MAKAAO Knowledge Graph", lang="en")))
    kg.add(
        (
            ONT,
            DCTERMS.description,
            Literal(
                "A FAIR-compliant RDF knowledge graph about autoantibodies, and autoimmune diseases",
                lang="en",
            ),
        )
    )

    for term in [
        "Knowledge Graph",
        "Autoantibodies",
        "Biomedical Ontology",
        "Autoimmune diseases",
    ]:
        kg.add((ONT, DCTERMS.subject, Literal(term, lang="en")))

    kg.add(
        (ONT, DCTERMS.license, URIRef("https://creativecommons.org/licenses/by/4.0/"))
    )
    kg.add((ONT, DCTERMS.accessRights, Literal("Open access")))
    kg.add((ONT, ODRL.hasPolicy, URIRef("https://makaao.inria.fr/usage_policies.html")))

    kg.add(
        (
            ONT,
            DCAT.downloadURL,
            URIRef(f"https://makaao.inria.fr/data/makaao_{version}.rdf"),
        )
    )
    kg.add((ONT, DCAT.endpointURL, URIRef("http://makaao.inria.fr/kg/")))
    kg.add((ONT, DCAT.mediaType, Literal("application/rdf+xml")))
    kg.add((ONT, DCAT.landingPage, URIRef("https://makaao.inria.fr/")))

    ACT = URIRef("http://makaao.inria.fr/kg/activity_makaao_core")
    kg.add((ONT, PROV.wasGeneratedBy, ACT))
    kg.add((ONT, PROV.wasDerivedFrom, ACT))

    author_uri = URIRef("https://heka.gitlabpages.inria.fr/team/members/maury.html")
    kg.add((ONT, DCTERMS.creator, author_uri))
    kg.add((author_uri, RDF.type, FOAF.Person))
    add_pref(kg, author_uri, "Fabien Maury")

    team_alpha = URIRef("https://heka.gitlabpages.inria.fr/")
    team_beta = URIRef("https://www.institutimagine.org/en")
    for team_uri, team_label in [
        (team_alpha, "Team HeKA, Inria Paris"),
        (team_beta, "Institut Imagine, Inserm"),
    ]:
        kg.add((ONT, DCTERMS.contributor, team_uri))
        kg.add((team_uri, RDF.type, FOAF.Organization))
        add_pref(kg, team_uri, team_label)

    kg.add((ONT, DCTERMS.created, Literal(date.today().isoformat(), datatype=XSD.date)))
    kg.add((ONT, OWL.versionInfo, Literal(version)))
    kg.add((ONT, OWL.imports, URIRef("http://purl.obolibrary.org/obo/ro.owl")))
    kg.add((ONT, VOID.uriSpace, Literal("http://makaao.inria.fr/kg/")))
    kg.add((ONT, SCHEMA.name, Literal("MAKAAO Knowledge Graph", lang="en")))
    kg.add(
        (
            ONT,
            SCHEMA.description,
            Literal("A knowledge graph for autoantibodies", lang="en"),
        )
    )
    kg.add((ONT, SCHEMA.url, URIRef("http://makaao.inria.fr/kg/")))
    kg.add(
        (ONT, SCHEMA.license, URIRef("https://creativecommons.org/licenses/by/4.0/"))
    )
    kg.add((ONT, RDFS.seeAlso, URIRef("https://makaao.inria.fr")))
    for kw in ["Autoantibodies", "Autoimmune diseases"]:
        kg.add((ONT, DCAT.keyword, Literal(kw, lang="en")))
    kg.add((ONT, DCAT.contactPoint, URIRef("mailto:contact.makaao@inria.fr")))
    kg.add((ONT, VOID.triples, Literal(len(kg)+1, datatype=XSD.integer)))




# ===================== T-BOX EXPORT =====================
def extract_tbox_local_only(source: Graph, local_ns: str) -> Graph:
    """
    Extract a T-box graph containing ONLY schema terms defined in the local namespace.

    Included:
      - triples where the subject is a local IRI (starts with local_ns) AND the subject is
        declared as owl:Ontology / owl:Class / owl:ObjectProperty / owl:DatatypeProperty / owl:AnnotationProperty
      - all outgoing triples from those local schema subjects
      - recursively, any blank-node structures reachable from those triples (e.g., OWL restrictions)

    Excluded:
      - any triples whose subject is an external IRI (even if present in the KG),
        including copied labels/types/subClassOf for external ontologies.
    """
    tbox = Graph()
    tbox.namespace_manager = source.namespace_manager

    schema_types = {
        OWL.Ontology,
        OWL.Class,
        OWL.ObjectProperty,
        OWL.DatatypeProperty,
        OWL.AnnotationProperty,
    }

    seeds = set()
    for s, o in source.subject_objects(RDF.type):
        if o in schema_types and isinstance(s, URIRef) and str(s).startswith(local_ns):
            seeds.add(s)

    # Copy all triples about local schema subjects; include connected blank nodes (e.g., restrictions)
    stack = list(seeds)
    seen = set(seeds)
    while stack:
        s = stack.pop()
        for _, p, o in source.triples((s, None, None)):
            tbox.add((s, p, o))
            if isinstance(o, BNode) and o not in seen:
                seen.add(o)
                stack.append(o)

    return tbox

# ===================== MAIN =====================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # create empty KG
    g = init_graph()

    # Optional enrichment datasets

    # read UMLS-HPO table where we can see the enlgish main names and syns associated to these code (this is an extrcat from MRCONSO.RRF)
    umls_hpo_rows = read_csv_rows(UMLS_HPO_CSV) if os.path.exists(UMLS_HPO_CSV) else []

    # read the csv files containing the UMLS main english name, associated to our UMLS CUI of interest (this was collected from UMLS API when we queried it four our UMLS concepts of interest) (not only HPO terms unlike the previous one)
    umls_names = (
        umls_names_index(UMLS_NAMES_CSV) if os.path.exists(UMLS_NAMES_CSV) else {}
    )

    # read UMLS-Orphanet csv table with their main english name and syns (obtained from MRCONSO.RRF)
    orpha_by_cui, orpha_by_code = (
        umls_code_group(UMLS_ORPHANET_CSV)
        if os.path.exists(UMLS_ORPHANET_CSV)
        else ({}, {})
    )
    # orpha_by_cui is a dict where keys are UMLS CUIs and values are list of rows (dict) from the csv where we can find the name/syns etc
    # orpha_by_code is similar but keys are Orphanet codes

    snomed_by_cui, _ = (
        umls_code_group(UMLS_SNOMED_CSV)
        if os.path.exists(UMLS_SNOMED_CSV)
        else ({}, {})
    )
    # snomed_by_cui is a dictionary where keys are UMLS CUIs and values are list of rows (dict) from the csv where we can find the name/syns etc (MRCONSO) (only rows of SNOMED items)

    # orpha_hpo_links is a dict where keys are Orphanet codes and values are [hpo, freq]; from the csv where we can find the HPO terms linked to the Orphanet diseases (from en_product4.xml an Orphanet file listing HPO terms associated to Orphanet dieases with their frequency)
    orpha_hpo_links = (
        grouped(read_csv_rows(ORPHANET_HPO_LINKS), "orpha_code")
        if os.path.exists(ORPHANET_HPO_LINKS)
        else {}
    )
    #print(orpha_hpo_links)

    # Names read from code_names.csv: a file obtained by querying UMLS for english names of our concepts of interest (targets, diseases...)
    if os.path.exists(CODE_NAMES_CSV):
        up_names, up_urls = read_code_names_uniprot(
            CODE_NAMES_CSV
        )  # 2 dictionnaries: {uniprot_id: english_name} and {uniprot_id: uniprot_url}

        umls_cn_names, _ = read_code_names_umls(
            CODE_NAMES_CSV
        )  # dictinoary: {CUI: english_name} (CUIs of differnt concept sof interest: targets, diseases...)
        chebi_cn_names, chebi_cn_urls = read_code_names_chebi(
            CODE_NAMES_CSV
        )  # 2 dictionaries: {chebi_id: english_name} and {chebi_id: chebi_url}
        print(
            f"Loaded names — UniProt:{len(up_names)} UMLS:{len(umls_cn_names)} ChEBI:{len(chebi_cn_names)}"
        )
    else:
        up_names, up_urls = {}, {}
        umls_cn_names = {}
        chebi_cn_names, chebi_cn_urls = {}, {}
        print(
            f"WARN: {CODE_NAMES_CSV} not found; labels may fall back to codes."
        )  # if code_names file not found, just diplay a warning

    # "data" is a dict where keys are column names, and values are also dict containing: different things, but they keys is alwaus aab_id. the values might be a list containing for ex [target,source]
    data = load_processed_tables(BASE_DIR)
    #print(data['hpo_list'])
    print("indices kept:", len(data["indices"]))
    if not data["indices"]:
        en = read_csv_rows(os.path.join(BASE_DIR, "index_name_en.csv"))
        print("index_name_en.csv headers seen:", list(en[0].keys()) if en else "EMPTY")
        raise SystemExit(
            "No indices kept. Ensure headers include 'index' and 'name_en'."
        )
    ################## end of data loading section ###########################################

    # this function call add AAb with their class hierarchy, and their positivity phenotypes. It also add their targets (UMLS, CHebi, or Uniprot), instantiate all these things and link them together as needed
    build_core(
        g,
        data,
        umls_hpo_rows,
        umls_names,
        up_names,
        up_urls,
        umls_cn_names,
        chebi_cn_names,
        chebi_cn_urls,
    )

    # this function call add diseases linked to AAbs, instantiate them and link them to AAb instances. Diseases can be Orphanet diseases (preferred) or UMLS CUIs (with possible SNOMED mappings)
    process_diseases(
        g,
        data,
        orpha_by_code,
        umls_names,
        orpha_hpo_links,
        umls_cn_names,
    )
    keep = data["indices"]
    process_loinc_mappings(g, LOINC_UMLS, LOINC_INDEX_CSV, keep)

    append_fair_metadata(g)

    # Export ontology (T-box) separately (no instances)
    tbox = extract_tbox_local_only(g, str(MAKAAO))
    tbox.serialize(destination=OUTPUT_OWL_TBOX, format="xml")
    print(f"Saved {OUTPUT_OWL_TBOX}  triples={len(tbox)}")

    g.serialize(
        destination=OUTPUT_OWL_ENRICHED, format="xml"
    )  #  we write the KG to an RDF file (XML syntax)
    print(f"Saved {OUTPUT_OWL_ENRICHED}  triples={len(g)}")


if __name__ == "__main__":  
    main()