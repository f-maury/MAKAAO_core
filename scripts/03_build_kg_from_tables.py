#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ===================== IMPORT =====================
import csv
import json
import re
import hashlib
import os
import importlib.util
import sys
import shutil
import shlex
import tempfile
import time
from pathlib import Path
from urllib.parse import quote, urlsplit, urlunsplit
from collections import defaultdict
from itertools import combinations
from rdflib import Graph, Namespace, RDF, RDFS, OWL, Literal, URIRef, BNode
from rdflib.namespace import XSD
from datetime import date

# ===================== SCRIPT VERSION =====================
SCRIPT_VERSION = "1.2.29"
SCRIPT_ITERATION = "2026-08-10-chebi-target-mappings-and-umls-altlabels"

# The dataset version is independent of the Python script version.
#KG_VERSION = "1.0.4test"  # from makaao_core_29-07-2026.xlsx
KG_VERSION = "sample"  # from makaao_sample.csv
version = KG_VERSION  # retained for compatibility with the existing code

# ===================== HARDCODED CONFIG =====================
# Paths are anchored to this script, so execution does not depend on the
# current working directory. Both layouts are supported:
#   project/03_build_kg*.py + project/data + project/kg
#   project/scripts/03_build_kg*.py + project/data + project/kg
SCRIPT_DIR = Path(__file__).resolve().parent


def _find_project_dir(start: Path) -> Path:
    """Find the project root from project/, scripts/, or scripts/03/.

    A pre-existing project ``kg/`` directory is not required; it is created
    during the build when necessary.
    """
    for candidate in (start, *start.parents):
        if (candidate / "data").is_dir():
            return candidate
    raise RuntimeError(
        f"Could not locate project root above {start}; expected a project data/ directory."
    )


PROJECT_DIR = _find_project_dir(SCRIPT_DIR)
DATA_DIR = PROJECT_DIR / "data"
# Write all generated KG and reasoning artifacts to the project's kg/ folder.
# This is intentionally not configurable through MAKAAO_KG_DIR.
KG_DIR = (PROJECT_DIR / "kg").resolve()

BASE_DIR = str(DATA_DIR / "processed_tables")
OUTPUT_DIR = str(DATA_DIR)
#makaao_core_name = str(DATA_DIR / "makaao_core.csv")
makaao_core_name = str(DATA_DIR / "makaao_sample.csv")
OUTPUT_OWL_ENRICHED = str(KG_DIR / f"makaao_kg_{version}.rdf")
OUTPUT_OWL_TBOX = str(KG_DIR / f"makaao_kg_{version}_ontology.owl")

# Optional enrichment tables.
ORPHANET_HPO_LINKS = str(DATA_DIR / "enrichment_tables" / "orphanet_hpo_links.csv")
LOINC_INDEX_CSV = str(DATA_DIR / "processed_tables" / "index_loinc.csv")
LOINC_PART_TEST_JSON = str(DATA_DIR / "enrichment_tables" / "loinc_part_test_dict.json")
CODE_NAMES_CSV = str(DATA_DIR / "enrichment_tables" / "code_names.csv")
ORPHA_UMLS_MAPPINGS_JSON = str(
    DATA_DIR / "enrichment_tables" / "orpha_umls_mappings.json"
)


# Pinned local OWL/RDF sources used to create minimal import modules.
EXTERNAL_ONTOLOGY_FILES = {
    "biolink": str(DATA_DIR / "external_ontologies" / "biolink-model.owl.ttl"),
    "hpo": str(DATA_DIR / "external_ontologies" / "hp.owl"),
    "chebi": str(DATA_DIR / "external_ontologies" / "chebi_lite.owl"),
    "ordo": str(DATA_DIR / "external_ontologies" / "ORDO_en_4.9.owl"),
    "sio": str(DATA_DIR / "external_ontologies" / "sio.owl"),
    "prov": str(DATA_DIR / "external_ontologies" / "prov-o.owl"),
    "uniprot": str(DATA_DIR / "external_ontologies" / "core.owl"),
    "bao": str(DATA_DIR / "external_ontologies" / "bao_complete_merged.owl"),
}
SKOS_SOURCE_FILE = str(DATA_DIR / "external_ontologies" / "skos-owl1-dl.rdf")

# ROBOT/HermiT are required for the strict reasoning release.
REASONER_TIMEOUT_SECONDS = int(
    os.environ.get("MAKAAO_REASONER_TIMEOUT_SECONDS", "1800")
)
JAVA_MAX_HEAP = os.environ.get("MAKAAO_JAVA_MAX_HEAP", "8G")
ROBOT_EXECUTABLE = "robot"
ROBOT_JAR = str(PROJECT_DIR / "robot.jar")

# Internal implementation loaded automatically by this script. Users run only
# this graph builder. The local project keeps the stable helper filename; script
# compatibility is enforced separately through the SCRIPT_VERSION check.
REASONING_RELEASE_BUILDER = str(SCRIPT_DIR / "03-1_build_reasoning_release.py")
BUILD_STRICT_REASONING_RELEASE = True
REQUIRE_STRICT_REASONING_RELEASE = True
REASONING_OUTPUT_DIRNAME = "reasoning"

# ===================== NAMESPACES =====================
MAKAAO = Namespace("http://makaao.inria.fr/kg/")
SKOS = Namespace("http://www.w3.org/2004/02/skos/core#")
SKOS_OWL_DL = URIRef("http://www.w3.org/TR/skos-reference/skos-owl1-dl.rdf")
PROV = Namespace("http://www.w3.org/ns/prov#")
HP = Namespace("http://purl.obolibrary.org/obo/")
BIOLINK = Namespace("https://w3id.org/biolink/vocab/")
BAO = Namespace("http://www.bioassayontology.org/bao#")
BAO_HAS_TARGET = BAO["BAO_0000211"]
BAO_IS_TARGET_FOR = BAO["BAO_0000598"]
BIOLINK_MODULE_IRI = URIRef("http://makaao.inria.fr/imports/biolink")
HPO_MODULE_IRI = URIRef("http://makaao.inria.fr/imports/hpo")
CHEBI_MODULE_IRI = URIRef("http://makaao.inria.fr/imports/chebi")
ORDO_MODULE_IRI = URIRef("http://makaao.inria.fr/imports/ordo")
SIO_MODULE_IRI = URIRef("http://makaao.inria.fr/imports/sio")
PROV_MODULE_IRI = URIRef("http://makaao.inria.fr/imports/prov")
BAO_MODULE_IRI = URIRef("http://makaao.inria.fr/imports/bao")
SIO = Namespace("http://semanticscience.org/resource/")
UMLS = Namespace("http://uts.nlm.nih.gov/uts/umls/concept/")
UP = Namespace("http://purl.uniprot.org/core/")
UNIPROT = Namespace("http://purl.uniprot.org/uniprot/")
UNIPROT_MODULE_IRI = URIRef("http://makaao.inria.fr/imports/uniprot-core")
EXTERNAL_IMPORT_IRIS = (
    BIOLINK_MODULE_IRI,
    HPO_MODULE_IRI,
    CHEBI_MODULE_IRI,
    ORDO_MODULE_IRI,
    SIO_MODULE_IRI,
    PROV_MODULE_IRI,
    BAO_MODULE_IRI,
    UNIPROT_MODULE_IRI,
)
LOINC = Namespace("https://loinc.org/")
LOINC_PROPERTY = Namespace("http://loinc.org/property/")
LOINC_COMPONENT = LOINC_PROPERTY["COMPONENT"]
MAKAAO_LOINC = Namespace("http://makaao.inria.fr/loinc/")

DCTERMS = Namespace("http://purl.org/dc/terms/")
DCAT = Namespace("http://www.w3.org/ns/dcat#")
ODRL = Namespace("http://www.w3.org/ns/odrl/2/")
SCHEMA = Namespace("http://schema.org/")
VOID = Namespace("http://rdfs.org/ns/void#")
FOAF = Namespace("http://xmlns.com/foaf/0.1/")

HPO_AUTOIMMUNE_ANTIBODY_POSITIVITY = HP["HP_0030057"]
LABEL_COLLISION_REPORT_FILENAME = "class-label-close-match-report.tsv"

# ===================== GLOBALS =====================
relation_counter = 0
document_map = {}
GLOBAL_ACTIVITY = None

# ---------- Label de-duplication state ----------
# We maintain, per-node:
# - a cross-property set of seen texts ("any")
# - per-property sets for SKOS.prefLabel, SKOS.altLabel and RDFS.label
# Rule: block duplicates globally, EXCEPT allow a single cross-duplicate
#       between rdfs:label and skos:prefLabel.
_label_seen = defaultdict(
    lambda: {
        "any": set(),
        SKOS.prefLabel: set(),
        SKOS.altLabel: set(),
        RDFS.label: set(),
    }
)

def ensure_uriref(val):  # small function to make sure a string is an URI
    if not val:
        return None  # Or raise an error
    if isinstance(val, URIRef):
        return val
    return require_absolute_iri(val, "RDF resource")


def _normalized_external_reference(value):
    """Return a valid external provenance IRI, or ``None`` for free text."""
    text = str(value or "").strip()
    if not text:
        return None

    pmid_match = pmid_rx.match(text)
    if pmid_match:
        return URIRef(f"https://pubmed.ncbi.nlm.nih.gov/{pmid_match.group(1)}")

    doi_match = re.fullmatch(r"(?i)(?:doi\s*:\s*)?(10\.\d{4,9}/\S+)", text)
    if doi_match:
        doi = doi_match.group(1).rstrip(".,;)")
        return URIRef("https://doi.org/" + quote(doi, safe="/:._-;()"))

    parsed = urlsplit(text)
    if parsed.scheme.lower() in {"http", "https"} and parsed.netloc:
        path = quote(parsed.path, safe="/%:@!$&'()*+,;=-._~")
        query = quote(parsed.query, safe="=&?/:@!$'()*+,;%-._~")
        fragment = quote(parsed.fragment, safe="/?@:!$&'()*+,;=%-._~")
        return URIRef(
            urlunsplit((parsed.scheme.lower(), parsed.netloc, path, query, fragment))
        )
    return None


_IRI_SCHEME_RX = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")
_IRI_FORBIDDEN_CHARS = frozenset('<>"{}|\\^`')
_INVALID_PERCENT_ESCAPE_RX = re.compile(r"%(?![0-9A-Fa-f]{2})")


def _iri_validation_error(value):
    """Return ``None`` for a usable absolute IRI, otherwise an error string."""
    text = str(value or "")
    if not text:
        return "IRI is empty"
    if any(character.isspace() or ord(character) < 0x20 for character in text):
        return "IRI contains whitespace or control characters"
    forbidden = sorted(character for character in _IRI_FORBIDDEN_CHARS if character in text)
    if forbidden:
        return "IRI contains forbidden characters: " + " ".join(forbidden)
    if _INVALID_PERCENT_ESCAPE_RX.search(text):
        return "IRI contains an invalid percent escape"
    if not _IRI_SCHEME_RX.match(text):
        return "IRI is relative or has no scheme"
    parsed = urlsplit(text)
    if parsed.scheme.lower() in {"http", "https"} and not parsed.netloc:
        return "HTTP(S) IRI has no authority/host"
    return None


def require_absolute_iri(value, context):
    """Validate and return an absolute ``URIRef``."""
    text = str(value or "").strip()
    error = _iri_validation_error(text)
    if error:
        raise ValueError(f"{context}: {error}: {text!r}")
    return URIRef(text)


def normalize_http_iri(value, context):
    """Normalize and validate an HTTP(S) URL from an input table."""
    text = str(value or "").strip()
    if not text:
        return None
    iri = _normalized_external_reference(text)
    if iri is None:
        raise ValueError(f"{context} must be an absolute HTTP(S) URL, got {text!r}")
    parsed = urlsplit(str(iri))
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise ValueError(f"{context} must be an absolute HTTP(S) URL, got {text!r}")
    error = _iri_validation_error(iri)
    if error:
        raise ValueError(f"{context}: {error}: {text!r}")
    return iri


def safe_local_fragment(value, prefix="entity", max_slug_length=48):
    """Create a deterministic, collision-resistant local IRI fragment."""
    raw = re.sub(r"\s+", " ", str(value or "").strip())
    if not raw:
        raise ValueError("Cannot create a local IRI fragment from an empty value")
    slug = re.sub(r"[^A-Za-z0-9._~-]+", "_", raw).strip("_")
    slug = re.sub(r"_+", "_", slug)[:max_slug_length].rstrip("_") or "value"
    digest = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:12]
    return f"{prefix}_{slug}_{digest}"


def validate_graph_iris(graph, context="graph"):
    """Reject relative, malformed, or whitespace-containing URIRefs."""
    problems = []
    seen = set()
    for subject, predicate, obj in graph:
        if not isinstance(predicate, URIRef):
            problems.append(f"predicate is not a URIRef in triple {(subject, predicate, obj)!r}")
        for position, node in (("subject", subject), ("predicate", predicate), ("object", obj)):
            if not isinstance(node, URIRef) or node in seen:
                continue
            seen.add(node)
            error = _iri_validation_error(node)
            if error:
                problems.append(f"{position} {str(node)!r}: {error}")
        if len(problems) >= 20:
            break
    if problems:
        raise ValueError(
            f"Invalid IRI(s) detected in {context}:\n  " + "\n  ".join(problems)
        )
    return len(seen)

def build_non_reified_graph(source: Graph) -> Graph:
    """Return a reasoner-oriented projection without reification records.

    The ordinary subject/predicate/object assertions are expected to already
    exist in ``source``. This function never reconstructs or materializes them
    from RDF reification triples. It validates that invariant and removes:

    * every instance of ``rdf:Statement`` or ``makaao:Relation``;
    * every instance of ``makaao:Document``; and
    * the RDF reification predicates ``rdf:subject``, ``rdf:predicate`` and
      ``rdf:object``.

    Schema declarations for ``makaao:Relation`` and ``makaao:Document`` remain
    part of the MAKAAO TBox; only their provenance-record instances are removed.
    The complete provenance-bearing KG is not modified and is serialized
    separately.
    """
    statement_classes = {RDF.Statement, MAKAAO.Relation}
    document_classes = {MAKAAO.Document}

    statement_nodes = set()
    for cls in statement_classes:
        statement_nodes.update(source.subjects(RDF.type, cls))

    document_nodes = set()
    for cls in document_classes:
        document_nodes.update(source.subjects(RDF.type, cls))

    removed_instances = statement_nodes | document_nodes
    removed_resources = removed_instances
    reification_predicates = {RDF.subject, RDF.predicate, RDF.object}

    # Reification is metadata about assertions already present in the KG. Do
    # not create assertions from it. Instead, fail clearly if that invariant is
    # violated so the non-reified output cannot silently acquire new facts.
    for relation in statement_nodes:
        subjects = list(source.objects(relation, RDF.subject))
        predicates = list(source.objects(relation, RDF.predicate))
        objects = list(source.objects(relation, RDF.object))
        if len(subjects) != 1 or len(predicates) != 1 or len(objects) != 1:
            raise ValueError(
                f"Malformed reified relation {relation}: expected exactly one "
                "rdf:subject, rdf:predicate and rdf:object"
            )
        subject = subjects[0]
        predicate = predicates[0]
        obj = objects[0]
        if not isinstance(predicate, URIRef):
            raise ValueError(
                f"Malformed reified relation {relation}: predicate is not an IRI"
            )
        if (
            subject in removed_resources
            or predicate in removed_resources
            or obj in removed_resources
        ):
            continue
        if predicate in reification_predicates:
            raise ValueError(
                f"Malformed reified relation {relation}: reifies RDF reification "
                f"predicate {predicate}"
            )
        if (subject, predicate, obj) not in source:
            raise ValueError(
                f"Reified relation {relation} has no corresponding ordinary triple: "
                f"{(subject, predicate, obj)!r}"
            )

    projected = Graph()
    projected.namespace_manager = source.namespace_manager
    for subject, predicate, obj in source:
        if (
            subject in removed_resources
            or predicate in removed_resources
            or obj in removed_resources
        ):
            continue
        if predicate in reification_predicates:
            continue
        projected.add((subject, predicate, obj))

    # Strong postconditions: no removed provenance-record instance may occur
    # anywhere in the reasoner-oriented projection.
    for removed in removed_resources:
        if any(projected.triples((removed, None, None))):
            raise RuntimeError(
                f"Non-reified projection still contains triples about removed resource {removed}"
            )
        if any(projected.triples((None, removed, None))):
            raise RuntimeError(
                f"Non-reified projection still uses removed resource as predicate {removed}"
            )
        if any(projected.triples((None, None, removed))):
            raise RuntimeError(
                f"Non-reified projection still references removed resource {removed}"
            )

    for predicate in reification_predicates:
        if any(projected.triples((None, predicate, None))):
            raise RuntimeError(
                f"Non-reified projection still uses reification predicate {predicate}"
            )

    for cls in statement_classes | document_classes:
        if any(projected.subjects(RDF.type, cls)):
            raise RuntimeError(
                f"Non-reified projection still contains instances of removed class {cls}"
            )

    validate_graph_iris(projected, "non-reified MAKAAO KG")
    return projected


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


def _is_loinc_part_class(graph: Graph, node: URIRef) -> bool:
    """Return whether ``node`` is one of the generated LOINC Part classes.

    LOINC Terms share the same IRI namespace, so namespace recognition alone is
    insufficient. Only classes explicitly placed under ``makaao:LoincPart`` by
    ``process_loinc_mappings`` are eligible for the collision audit.
    """
    return (
        isinstance(node, URIRef)
        and str(node).startswith(str(MAKAAO_LOINC))
        and (node, RDFS.subClassOf, MAKAAO.LoincPart) in graph
    )


def _is_chebi_class(node: URIRef) -> bool:
    """Return whether ``node`` is an external ChEBI class IRI."""
    return (
        isinstance(node, URIRef)
        and str(node).startswith("http://purl.obolibrary.org/obo/CHEBI_")
    )


def _is_target_class(graph: Graph, class_iri: URIRef) -> bool:
    """Return whether an individual of ``class_iri`` is used as a MAKAAO Target."""
    return any(
        (instance, RDF.type, MAKAAO.Target) in graph
        for instance in graph.subjects(RDF.type, class_iri)
    )


def _label_collision_resource_kind(graph: Graph, node: URIRef) -> str:
    """Return the generated/external resource family used by close-match rules."""
    text = str(node)
    local_tail = text[len(str(MAKAAO)) :] if text.startswith(str(MAKAAO)) else ""
    if re.fullmatch(r"aab_\d+", local_tail):
        return "autoantibody"
    if re.fullmatch(r"CUI_C\d+", local_tail):
        return "umls_cui"
    if re.fullmatch(r"UP_[A-Za-z0-9]+", local_tail):
        return "uniprot"
    if _is_chebi_class(node):
        return "chebi"
    if _is_loinc_part_class(graph, node):
        return "loinc_part"
    return "other"

def _stable_text_key(value: str) -> tuple[str, str]:
    """Return a total, reproducible ordering key for report strings."""
    return value.casefold(), value


def _class_label_index(graph: Graph):
    """Index all supported labels on eligible explicitly declared classes.

    The audit compares ``rdfs:label``, ``skos:prefLabel`` and
    ``skos:altLabel``. Generated LOINC Terms are excluded; only generated LOINC
    Part classes participate.
    """
    classes = {
        node
        for class_type in (OWL.Class, RDFS.Class)
        for node in graph.subjects(RDF.type, class_type)
        if isinstance(node, URIRef)
        and (
            not str(node).startswith(str(MAKAAO_LOINC))
            or _is_loinc_part_class(graph, node)
        )
    }
    labels_by_class = defaultdict(lambda: defaultdict(set))
    classes_by_label = defaultdict(set)

    for class_iri in classes:
        for predicate in (SKOS.prefLabel, RDFS.label, SKOS.altLabel):
            for value in graph.objects(class_iri, predicate):
                if not isinstance(value, Literal):
                    continue
                normalized = _norm_label_text(str(value))
                if not normalized:
                    continue
                labels_by_class[class_iri][normalized].add(
                    (str(predicate), str(value), value.language or "")
                )
                classes_by_label[normalized].add(class_iri)

    return labels_by_class, classes_by_label


def add_label_collision_close_matches(graph: Graph) -> list[dict[str, str | int]]:
    """Audit normalized class-label collisions and add permitted mappings.

    The function is called once after all input-derived graph construction and
    before TBox extraction or reasoning. Distinct explicitly declared classes
    are candidates when any ``rdfs:label``, ``skos:prefLabel`` or
    ``skos:altLabel`` matches after whitespace normalization and case folding.

    Symmetric ``skos:closeMatch`` links are added only for:

    * distinct UMLS CUI, UniProt and ChEBI classes when both classes are
      evidenced as MAKAAO targets; or
    * a MAKAAO autoantibody class and a generated LOINC Part class.

    No specific resource pairs are hard-coded. All other collisions remain
    report-only because equal labels can denote concepts with different scopes.
    """
    labels_by_class, classes_by_label = _class_label_index(graph)
    shared_by_pair = defaultdict(set)
    for normalized_label, class_iris in classes_by_label.items():
        for left, right in combinations(sorted(class_iris, key=str), 2):
            shared_by_pair[(left, right)].add(normalized_label)

    report_rows: list[dict[str, str | int]] = []
    for (left, right), shared_labels in sorted(
        shared_by_pair.items(), key=lambda item: (str(item[0][0]), str(item[0][1]))
    ):
        left_kind = _label_collision_resource_kind(graph, left)
        right_kind = _label_collision_resource_kind(graph, right)
        kind_pair = {left_kind, right_kind}

        decision = "review_required"
        rule = "same normalized label alone is insufficient evidence"
        target_kinds = {"umls_cui", "uniprot", "chebi"}
        if (
            left_kind in target_kinds
            and right_kind in target_kinds
            and left_kind != right_kind
            and _is_target_class(graph, left)
            and _is_target_class(graph, right)
        ):
            decision = "linked"
            rule = (
                "exact normalized label match across distinct MAKAAO target "
                "classes (UMLS CUI / UniProt / ChEBI)"
            )
        elif kind_pair == {"autoantibody", "loinc_part"}:
            decision = "linked"
            rule = "any exact normalized label match across autoantibody and LOINC Part"

        triples_added = 0
        if decision == "linked":
            for subject, obj in ((left, right), (right, left)):
                triple = (subject, SKOS.closeMatch, obj)
                if triple not in graph:
                    graph.add(triple)
                    triples_added += 1

        left_matching = sorted(
            {
                value
                for normalized in shared_labels
                for _, value, _ in labels_by_class[left][normalized]
            },
            key=_stable_text_key,
        )
        right_matching = sorted(
            {
                value
                for normalized in shared_labels
                for _, value, _ in labels_by_class[right][normalized]
            },
            key=_stable_text_key,
        )
        report_rows.append(
            {
                "decision": decision,
                "rule": rule,
                "class_1": str(left),
                "class_1_kind": left_kind,
                "class_1_matching_labels": " | ".join(left_matching),
                "class_2": str(right),
                "class_2_kind": right_kind,
                "class_2_matching_labels": " | ".join(right_matching),
                "shared_normalized_labels": " | ".join(
                    sorted(shared_labels, key=_stable_text_key)
                ),
                "close_match_triples_added": triples_added,
            }
        )

    return report_rows


def write_label_collision_report(
    rows: list[dict[str, str | int]], report_path: Path
) -> None:
    """Write the deterministic class-label collision and mapping report."""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "decision",
        "rule",
        "class_1",
        "class_1_kind",
        "class_1_matching_labels",
        "class_2",
        "class_2_kind",
        "class_2_matching_labels",
        "shared_normalized_labels",
        "close_match_triples_added",
    ]
    with report_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# Regex to parse simple "PMID: 12345" tokens.
pmid_rx = re.compile(
    r"^\s*PMID\s*:\s*(\d+)\s*$", re.IGNORECASE
)  # matches "PMID: 12345" format

# ===================== IO / SMALL UTILS =====================
def read_csv_rows(path):
    """
    CSV → list[dict] with normalized original-case and lowercase keys.

    UTF-8 BOMs and surrounding header/value whitespace are removed. Values are
    read using the raw DictReader field names, so headers such as `` index ``
    are normalized without losing their corresponding cells.
    """
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as f:
        rdr = csv.DictReader(f)
        raw_headers = list(rdr.fieldnames or [])
        normalized_headers = [
            (header or "").lstrip("\ufeff").strip() for header in raw_headers
        ]
        nonempty = [header for header in normalized_headers if header]
        header_counts = defaultdict(int)
        display_header = {}
        for header in nonempty:
            key = header.casefold()
            header_counts[key] += 1
            display_header.setdefault(key, header)
        duplicate_headers = sorted(
            display_header[key] for key, count in header_counts.items() if count > 1
        )
        if duplicate_headers:
            raise ValueError(
                f"Duplicate CSV headers after whitespace/case normalization in {path}: "
                + ", ".join(duplicate_headers)
            )

        out = []
        for raw in rdr:
            row = {}
            for raw_header, normalized in zip(raw_headers, normalized_headers):
                if not normalized:
                    continue
                value = (raw.get(raw_header) or "").strip()
                row[normalized] = value
                row[normalized.lower()] = value
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


def canonical_loinc_code(value, kind=None):
    """Normalize a LOINC Part or Term identifier to its bare canonical code."""
    s = str(value or "").strip()
    if not s:
        return None

    # Accept bare codes, ``LOINC:``-prefixed codes, and canonical LOINC URLs.
    s = re.sub(r"(?i)^https?://(?:www\.)?loinc\.org/", "", s)
    s = s.split("?", 1)[0].split("#", 1)[0].strip("/")
    s = re.sub(r"(?i)^LOINC\s*:\s*", "", s).strip()
    s = re.sub(r"[\u2010-\u2015\u2212\uFE63\uFF0D]", "-", s).upper()

    if kind == "part":
        pattern = r"LP\d+(?:-\d+)?"
    elif kind == "term":
        pattern = r"\d+-\d+"
    else:
        pattern = r"(?:LP\d+(?:-\d+)?|\d+-\d+)"
    return s if re.fullmatch(pattern, s) else None

def canonical_hpo_code(hp_code):
    """Normalize many HPO formats to HP:NNNNNNN."""
    s = (hp_code or "").strip()
    if not s:
        return None
    if s.lower().startswith("http"):
        s = s.rsplit("/", 1)[-1].split("?")[0].split("#")[0]
    s = s.replace("_", ":").upper()
    m = re.search(r"HP:(\d{7})", s)
    if m:
        return f"HP:{m.group(1)}"
    m = re.search(r"HP(\d{7})", s)
    if m:
        return f"HP:{m.group(1)}"
    m = re.search(r"(\d{7})", s)
    if m:
        return f"HP:{m.group(1)}"
    return None

def canonical_orpha_code(value, allow_bare=False):
    """Normalize explicit Orphanet identifiers to their numeric code."""
    text = str(value or "").strip()
    if not text:
        return None
    upper = text.upper()
    if upper.startswith(("HTTP://", "HTTPS://")):
        parsed = urlsplit(text)
        tail = parsed.path.rstrip("/").rsplit("/", 1)[-1]
        match = re.fullmatch(r"(?i)ORPHANET[_:-]?(\d+)", tail)
        if match:
            return match.group(1)
        # Orphanet's public disease pages use paths ending in a bare code,
        # for example https://www.orpha.net/en/disease/detail/999.
        if parsed.hostname and "orpha" in parsed.hostname.lower() and re.fullmatch(r"\d+", tail):
            return tail
        return None
    match = re.fullmatch(r"(?i)(?:ORPHA|ORPHANET)[_:-]?(\d+)", text)
    if match:
        return match.group(1)
    if allow_bare and re.fullmatch(r"\d+", text):
        return text
    return None


def split_hpo_values(cell):
    """Split a cell that may contain 0/1/N HPO IDs separated by | , ; or whitespace."""
    txt = (cell or "").strip()
    if not txt:
        return []
    parts = [p.strip() for p in re.split(r"[|,;\s]+", txt) if p.strip()]
    out, seen = [], set()
    for p in parts:
        c = canonical_hpo_code(p)
        if c and c not in seen:
            seen.add(c)
            out.append(c)
    return out

def hp_to_obo_uri(hp_code):
    c = canonical_hpo_code(hp_code)
    if not c:
        return None
    return URIRef(f"http://purl.obolibrary.org/obo/{c.replace(':', '_')}")  # returns HP OBO URI

# ===================== UniProt / ChEBI HELPERS =====================
# Remove an optional UniProt prefix. Bare accessions are valid input too.
def _strip_prefixes(s):
    s = str(s).strip()
    return s[3:] if s.upper().startswith("UP:") else s

# Official UniProtKB accession patterns: 6 or 10 characters, optionally
# followed by a positive isoform suffix (for example P12345-2).
_UNIPROT_ACCESSION_RX = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
    r"(?:-[1-9][0-9]*)?"
)

def canon_uniprot_id(val):
    s = (val or "").strip()
    alts = set()
    if not s:
        return "", "", []
    if s.lower().startswith(("http://", "https://")):
        parsed = urlsplit(s)
        path_candidates = []
        for part in parsed.path.split("/"):
            if not part.strip():
                continue
            candidate = _strip_prefixes(part).upper()
            for suffix in (".FASTA", ".TXT", ".JSON", ".XML", ".RDF"):
                if candidate.endswith(suffix):
                    candidate = candidate[: -len(suffix)]
                    break
            path_candidates.append(candidate)
        s = next(
            (
                candidate
                for candidate in reversed(path_candidates)
                if _UNIPROT_ACCESSION_RX.fullmatch(candidate)
            ),
            "",
        )
    if "|" in s:
        parts = [part.strip() for part in s.split("|")]
        alts.update(part for part in parts if part)
        if len(parts) >= 2 and parts[1]:
            s = parts[1]
        else:
            s = next((part for part in reversed(parts) if part), "")
    s = _strip_prefixes(s).strip()
    s = s.split()[0]
    s = s.split(";", 1)[0].strip().upper()
    if not _UNIPROT_ACCESSION_RX.fullmatch(s):
        raise ValueError(
            f"Invalid UniProtKB accession {val!r}; expected an official 6- or "
            "10-character accession, optionally followed by an isoform suffix."
        )
    base = s.split("-", 1)[0]
    if s != base:
        alts.add(s)
    original = (val or "").strip()
    if original and original.upper() != s:
        alts.add(original)
    return (
        base,
        s,
        sorted(alt for alt in alts if alt),
    )  # returns uniprot base, uniprot full, [alt ids]

_CHEBI_ACCESSION_RX = re.compile(r"^CHEBI[_:](\d+)$", re.IGNORECASE)

def canon_chebi_id(val):
    s = (val or "").strip()
    if not s:
        raise ValueError("Empty ChEBI value.")
    if s.lower().startswith(("http://", "https://")):
        parsed = urlsplit(s)
        s = parsed.path.rstrip("/").rsplit("/", 1)[-1]
    s = s.split("?", 1)[0].split("#", 1)[0].strip()
    match = _CHEBI_ACCESSION_RX.fullmatch(s)
    if not match:
        raise ValueError(
            "Invalid ChEBI identifier "
            f"{val!r}; expected CHEBI:1234, CHEBI_1234, or an OBO URL ending "
            "in CHEBI_1234."
        )
    num = match.group(1)
    return f"CHEBI:{num}", f"CHEBI_{num}"  # returns (CHEBI:nnnnn, CHEBI_nnnnn)

def canonical_cui_code(val):
    """Normalize common CUI formats to C########."""
    s = (val or "").strip().upper()
    if not s:
        return None
    if s.startswith("HTTP"):
        s = s.rsplit("/", 1)[-1].split("?")[0].split("#")[0]
    m = re.search(r"C\d{7,8}", s)
    return m.group(0) if m else None

# ===================== code_names.csv READERS =====================

# build uniprot lookup table from reading csv
def read_code_names_uniprot(path):
    id2name, id2url = {}, {}
    rows = read_csv_rows(path)
    for row_number, r in enumerate(rows, start=2):
        src = (r.get("source") or "").strip().lower()
        if src != "uniprot":
            continue
        uid = (r.get("id") or "").strip().upper()
        nm = (r.get("name") or "").strip()
        raw_url = (r.get("url") or "").strip()
        url_iri = normalize_http_iri(
            raw_url, f"{path}: row {row_number} UniProt URL"
        ) if raw_url else None
        url = str(url_iri) if url_iri is not None else ""
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
    for row_number, r in enumerate(rows, start=2):
        if (r.get("source") or "").strip().lower() != "umls":
            continue
        cui = (r.get("id") or "").strip().upper().replace("CUI:", "")
        nm = (r.get("name") or "").strip()
        raw_url = (r.get("url") or "").strip()
        url_iri = normalize_http_iri(
            raw_url, f"{path}: row {row_number} UMLS URL"
        ) if raw_url else None
        url = str(url_iri) if url_iri is not None else ""
        if cui and nm and cui not in id2name:
            id2name[cui] = nm
        if cui and url and cui not in id2url:
            id2url[cui] = url

    return id2name, id2url  # 2 dicts: CUI -> name, CUI -> url

def read_code_names_umls_synonyms(path):
    """Read CUI -> English alternative labels from ``code_names.csv``."""
    result = {}
    for row_number, r in enumerate(read_csv_rows(path), start=2):
        if (r.get("source") or "").strip().lower() != "umls":
            continue
        cui = canonical_cui_code(r.get("id"))
        if not cui:
            continue
        raw = (r.get("synonyms_en") or "").strip()
        if not raw:
            continue
        try:
            values = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"{path}: row {row_number} synonyms_en must be a JSON array: {exc}"
            ) from exc
        if not isinstance(values, list):
            raise ValueError(
                f"{path}: row {row_number} synonyms_en must be a JSON array"
            )

        preferred_norm = _norm_label_text(r.get("name") or "")
        seen = set()
        labels = []
        for value in values:
            label = re.sub(r"\s+", " ", str(value or "")).strip()
            normalized = _norm_label_text(label)
            if (
                not label
                or not normalized
                or normalized == preferred_norm
                or normalized in seen
            ):
                continue
            seen.add(normalized)
            labels.append(label)
        if labels:
            result[cui] = tuple(labels)
    return result

def read_code_names_orpha(path):
    id2name, id2url = {}, {}
    rows = read_csv_rows(path)
    for row_number, r in enumerate(rows, start=2):
        src = (r.get("source") or "").strip().lower()
        if "orpha" not in src:
            continue
        raw = (r.get("id") or "").strip()
        nm = (r.get("name") or "").strip()
        raw_url = (r.get("url") or "").strip()
        url_iri = normalize_http_iri(
            raw_url, f"{path}: row {row_number} Orphanet URL"
        ) if raw_url else None
        url = str(url_iri) if url_iri is not None else ""
        orpha_num = canonical_orpha_code(raw, allow_bare=True)
        if nm and orpha_num and orpha_num not in id2name:
            id2name[orpha_num] = nm
        if url and orpha_num and orpha_num not in id2url:
            id2url[orpha_num] = url

    return id2name, id2url

def read_code_names_chebi(path):
    id2name, id2url = {}, {}
    rows = read_csv_rows(path)
    for row_number, r in enumerate(rows, start=2):
        src = (r.get("source") or "").strip().lower()
        if "chebi" not in src and src != "che":
            continue
        raw = (r.get("id") or "").strip()
        nm = (r.get("name") or "").strip()
        raw_url = (r.get("url") or "").strip()
        url_iri = normalize_http_iri(
            raw_url, f"{path}: row {row_number} ChEBI URL"
        ) if raw_url else None
        url = str(url_iri) if url_iri is not None else ""
        code_colon, _ = canon_chebi_id(raw)
        if nm and code_colon and code_colon not in id2name:
            id2name[code_colon] = nm
        if url and code_colon and code_colon not in id2url:
            id2url[code_colon] = url

    return id2name, id2url  # 2 dicts: CHeBI code -> name, CHeBI code -> url

def read_code_names_hpo(path):
    id2name, id2url = {}, {}
    rows = read_csv_rows(path)
    for row_number, r in enumerate(rows, start=2):
        src = (r.get("source") or "").strip().lower()
        if src != "hpo":
            continue
        raw = (r.get("id") or "").strip()
        nm = (r.get("name") or "").strip()
        raw_url = (r.get("url") or "").strip()
        url_iri = normalize_http_iri(
            raw_url, f"{path}: row {row_number} HPO URL"
        ) if raw_url else None
        url = str(url_iri) if url_iri is not None else ""
        code = canonical_hpo_code(raw)
        if not code:
            continue
        if nm and code not in id2name:
            id2name[code] = nm
        if url and code not in id2url:
            id2url[code] = url
    return id2name, id2url  # 2 dicts: HPO code -> name, HPO code -> url


def read_code_names_loinc(path):
    """Return LOINC Part and Term labels from ``code_names.csv``."""
    part_names, term_names = {}, {}
    rows = read_csv_rows(path)
    for row_number, row in enumerate(rows, start=2):
        if (row.get("source") or "").strip().lower() != "loinc":
            continue

        raw_code = (row.get("id") or "").strip()
        label = (row.get("name") or "").strip()
        code = canonical_loinc_code(raw_code)
        if not code:
            raise ValueError(
                f"{path}: row {row_number} has an invalid LOINC identifier: "
                f"{raw_code!r}"
            )

        names = part_names if code.startswith("LP") else term_names
        if code in names and names[code] != label:
            raise ValueError(
                f"{path}: conflicting labels for LOINC {code}: "
                f"{names[code]!r} and {label!r}"
            )
        names.setdefault(code, label)

    return part_names, term_names

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
        external_reference = _normalized_external_reference(prov)
        if external_reference is not None:
            g.add((doc, RDFS.seeAlso, external_reference))
        else:
            # Preserve non-IRI provenance text without constructing an invalid URI.
            g.add((doc, DCTERMS.identifier, Literal(prov)))

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

def ensure_cui_class(
    g, cui_code, preferred_label=None, source_url=None, synonyms=None
):
    """Create or complete one local class for a canonical UMLS CUI.

    The selected preferred name is attached as ``skos:prefLabel`` and
    ``rdfs:label``. Distinct English MRCONSO alternatives are attached as
    ``skos:altLabel`` values.
    """
    code = canonical_cui_code(cui_code)
    if not code:
        return None
    cui_cls = MAKAAO[f"CUI_{code}"]
    g.add((cui_cls, RDF.type, OWL.Class))
    g.add((cui_cls, RDFS.subClassOf, MAKAAO.CUI))
    label = str(preferred_label or "").strip() or code
    add_pref(g, cui_cls, label)
    for synonym in synonyms or ():
        add_label(g, cui_cls, SKOS.altLabel, synonym)
    external = (
        normalize_http_iri(source_url, f"UMLS URL for {code}")
        if source_url
        else UMLS[code]
    )
    g.add((cui_cls, RDFS.seeAlso, external))
    return cui_cls


def ensure_cui_class_and_instance(
    g, cui_code, preferred_label=None, source_url=None, synonyms=None
):
    """Create one local CUI class and one local instance for a canonical CUI."""
    code = canonical_cui_code(cui_code)
    if not code:
        return None, None
    cui_cls = ensure_cui_class(
        g,
        code,
        preferred_label=preferred_label,
        source_url=source_url,
        synonyms=synonyms,
    )
    cui_inst = MAKAAO[f"CUI_{code}_instance"]
    g.add((cui_inst, RDF.type, cui_cls))
    label = str(preferred_label or "").strip() or code
    add_pref(g, cui_inst, label)
    for synonym in synonyms or ():
        add_label(g, cui_inst, SKOS.altLabel, synonym)
    external = (
        normalize_http_iri(source_url, f"UMLS URL for {code}")
        if source_url
        else UMLS[code]
    )
    g.add((cui_inst, RDFS.seeAlso, external))
    return cui_cls, cui_inst

def validate_local_cui_labels(graph):
    """Require labels on every generated local CUI class and instance."""
    prefix = str(MAKAAO) + "CUI_C"
    unlabeled = []
    candidates = {
        node
        for triple in graph
        for node in (triple[0], triple[2])
        if isinstance(node, URIRef) and str(node).startswith(prefix)
    }
    for node in sorted(candidates, key=str):
        if not any(graph.objects(node, RDFS.label)) and not any(
            graph.objects(node, SKOS.prefLabel)
        ):
            unlabeled.append(str(node))
    if unlabeled:
        preview = ", ".join(unlabeled[:20])
        suffix = "" if len(unlabeled) <= 20 else f" (+{len(unlabeled) - 20} more)"
        raise RuntimeError(
            "Generated local UMLS resources without labels: "
            + preview
            + suffix
        )
    return len(candidates)


def read_orpha_umls_mappings(path):
    """Read and validate the symmetric ORPHA/UMLS enrichment dictionary.

    Script 02 writes a deterministic JSON object containing ``orpha_to_umls``
    and ``umls_to_orpha`` dictionaries. Both directions are required and must
    describe exactly the same set of identifier pairs.
    """
    mapping_path = Path(path)
    if not mapping_path.is_file():
        raise FileNotFoundError(
            "Required ORPHA/UMLS enrichment dictionary is missing: "
            f"{mapping_path}. Run script 02 before script 03."
        )
    try:
        payload = json.loads(mapping_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(
            f"Could not read ORPHA/UMLS enrichment dictionary {mapping_path}: {exc}"
        ) from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{mapping_path} must contain a JSON object")

    def pairs_from(direction, key_kind):
        raw = payload.get(direction)
        if not isinstance(raw, dict):
            raise ValueError(f"{mapping_path}: {direction!r} must be a JSON object")
        pairs = set()
        for raw_key, raw_values in raw.items():
            values = raw_values if isinstance(raw_values, list) else [raw_values]
            if key_kind == "orpha":
                orpha = canonical_orpha_code(raw_key, allow_bare=True)
                if not orpha:
                    raise ValueError(
                        f"{mapping_path}: invalid Orphanet key in {direction}: {raw_key!r}"
                    )
                for raw_cui in values:
                    cui = canonical_cui_code(str(raw_cui or ""))
                    if not cui:
                        raise ValueError(
                            f"{mapping_path}: invalid CUI mapped from ORPHA:{orpha}: "
                            f"{raw_cui!r}"
                        )
                    pairs.add((orpha, cui))
            else:
                cui = canonical_cui_code(str(raw_key or ""))
                if not cui:
                    raise ValueError(
                        f"{mapping_path}: invalid CUI key in {direction}: {raw_key!r}"
                    )
                for raw_orpha in values:
                    orpha = canonical_orpha_code(raw_orpha, allow_bare=True)
                    if not orpha:
                        raise ValueError(
                            f"{mapping_path}: invalid Orphanet code mapped from {cui}: "
                            f"{raw_orpha!r}"
                        )
                    pairs.add((orpha, cui))
        return pairs

    forward = pairs_from("orpha_to_umls", "orpha")
    reverse = pairs_from("umls_to_orpha", "umls")
    if forward != reverse:
        only_forward = sorted(forward - reverse)[:20]
        only_reverse = sorted(reverse - forward)[:20]
        raise ValueError(
            f"{mapping_path}: ORPHA/UMLS dictionaries are not symmetric; "
            f"only in orpha_to_umls={only_forward}, "
            f"only in umls_to_orpha={only_reverse}"
        )
    return sorted(forward, key=lambda pair: (int(pair[0]), pair[1]))


def add_orpha_umls_close_matches(
    graph,
    mapping_pairs,
    umls_names=None,
    umls_urls=None,
    orpha_names=None,
    umls_synonyms=None,
):
    """Materialize symmetric class-level ORPHA/UMLS ``skos:closeMatch`` links.

    Mapping-only CUIs are labelled from their mapped Orphanet concept when no
    direct UMLS label was generated. The canonical CUI remains the final
    fallback through :func:`ensure_cui_class`.
    """
    umls_names = umls_names or {}
    umls_urls = umls_urls or {}
    orpha_names = orpha_names or {}
    umls_synonyms = umls_synonyms or {}
    triples_added = 0
    for orpha, cui in mapping_pairs:
        orpha_cls = URIRef(f"http://www.orpha.net/ORDO/Orphanet_{orpha}")
        graph.add((orpha_cls, RDF.type, OWL.Class))
        graph.add((orpha_cls, RDFS.subClassOf, MAKAAO.AutoimmuneDisease))
        cui_cls = ensure_cui_class(
            graph,
            cui,
            preferred_label=umls_names.get(cui) or orpha_names.get(orpha),
            source_url=umls_urls.get(cui),
            synonyms=umls_synonyms.get(cui, ()),
        )
        # ORDO and UMLS disease concepts use the same class-level MAKAAO role.
        graph.add((cui_cls, RDFS.subClassOf, MAKAAO.AutoimmuneDisease))
        for subject, obj in ((orpha_cls, cui_cls), (cui_cls, orpha_cls)):
            triple = (subject, SKOS.closeMatch, obj)
            if triple not in graph:
                graph.add(triple)
                triples_added += 1
    return triples_added


# ===================== GRAPH BUILDERS =====================
def init_graph():  # start an empty knowledge graph and add a few basic things to it
    global relation_counter, document_map, GLOBAL_ACTIVITY, _label_seen
    relation_counter = 0
    document_map = {}
    GLOBAL_ACTIVITY = None
    _label_seen = defaultdict(
        lambda: {
            "any": set(),
            SKOS.prefLabel: set(),
            SKOS.altLabel: set(),
            RDFS.label: set(),
        }
    )
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
    g.bind("up", UP)
    g.bind("uniprot", UNIPROT)
    g.bind("dcterms", DCTERMS)
    g.bind("dcat", DCAT)
    g.bind("odrl", ODRL)
    g.bind("schema", SCHEMA)
    g.bind("void", VOID)
    g.bind("foaf", FOAF)
    g.bind("loinc", LOINC)
    g.bind("loinc_property", LOINC_PROPERTY)
    g.bind("mloinc", MAKAAO_LOINC)

    # Ontology header
    onto = URIRef("http://makaao.inria.fr/kg/")
    g.add((onto, RDF.type, OWL.Ontology))
    add_pref(g, onto, "MAKAAO knowledge graph")
    g.add((onto, OWL.versionIRI, URIRef(f"http://makaao.inria.fr/kg/{version}")))
    g.add((onto, OWL.versionInfo, Literal(version)))
    g.add((onto, OWL.imports, SKOS_OWL_DL))
    for module_iri in EXTERNAL_IMPORT_IRIS:
        g.add((onto, OWL.imports, module_iri))

    # we manually define the main classes
    g.add((MAKAAO.Relation, RDF.type, OWL.Class))
    g.add((MAKAAO.Relation, RDFS.subClassOf, PROV.Entity))  # sublass of PROV Entity, used for reified relations that carry provenance information
    g.add((MAKAAO.Relation, RDFS.label, Literal("Relation")))
    g.add((MAKAAO.Document, RDFS.label, Literal("Document")))
    g.add((MAKAAO.Document, RDF.type, OWL.Class))
    g.add((MAKAAO.Document, RDFS.subClassOf, PROV.Entity))  # sublass of PROV Entity, used for reified relations that carry provenance information
    g.add((MAKAAO.BiomolecularEntity, RDF.type, OWL.Class))
    g.add((MAKAAO.BiomolecularEntity, RDFS.label, Literal("Biomolecular entity")))
    g.add(
        (
            MAKAAO.BiomolecularEntity,
            RDFS.subClassOf,
            BIOLINK.ChemicalEntityOrGeneOrGeneProduct,
        )
    )
    g.add((MAKAAO.Autoantibody, RDF.type, OWL.Class))
    g.add((MAKAAO.Autoantibody, RDFS.label, Literal("Autoantibody")))
    g.add((MAKAAO.Autoantibody, RDFS.subClassOf, MAKAAO.BiomolecularEntity))
    g.add((MAKAAO.AutoantibodyPositivity, RDF.type, OWL.Class))
    add_pref(
        g,
        MAKAAO.AutoantibodyPositivity,
        "Autoimmune antibody positivity",
    )
    g.add(
        (
            MAKAAO.AutoantibodyPositivity,
            RDFS.subClassOf,
            BIOLINK.PhenotypicFeature,
        )
    )
    g.add(
        (
            MAKAAO.AutoantibodyPositivity,
            SKOS.closeMatch,
            HPO_AUTOIMMUNE_ANTIBODY_POSITIVITY,
        )
    )
    g.add(
        (
            HPO_AUTOIMMUNE_ANTIBODY_POSITIVITY,
            SKOS.closeMatch,
            MAKAAO.AutoantibodyPositivity,
        )
    )
    g.add((MAKAAO.AutoimmuneDisease, RDF.type, OWL.Class))
    g.add((MAKAAO.AutoimmuneDisease, RDFS.label, Literal("Autoimmune disease")))
    g.add((MAKAAO.Target, RDF.type, OWL.Class))
    g.add((MAKAAO.Target, RDFS.label, Literal("Target")))
    g.add((MAKAAO.CUI, RDF.type, OWL.Class))
    g.add((MAKAAO.CUI, RDFS.label, Literal("UMLS concept unique identifier class")))

    # Molecular-target relationships use the two official BAO predicates only.
    # Their declarations, labels, domains, ranges and inverse axiom are supplied
    # by the pinned BAO reasoning module extracted from bao_complete_merged.owl.

    # Global PROV activity
    act = MAKAAO["activity_makaao_core"]
    g.add((act, RDF.type, PROV.Activity))
    add_pref(g, act, "MAKAAO Core processing activity")
    version_token = SCRIPT_VERSION.replace(".", "_")
    agent = MAKAAO[f"agent_03_build_kg_v{version_token}"]  # prov:SoftwareAgent
    g.add((agent, RDF.type, PROV.SoftwareAgent))
    add_pref(g, agent, Path(__file__).name)
    g.add((agent, DCTERMS.hasVersion, Literal(SCRIPT_VERSION)))
    g.add((agent, DCTERMS.identifier, Literal(SCRIPT_ITERATION)))
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
    for row_number, r in enumerate(syn_en, start=2):
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        se = (r.get("syns_en") or "").strip()
        src = (r.get("syns_en_source") or "").strip()
        if idx and se:
            normalized = to_pubmed_urls(
                src, src_file=syn_en_path, row=row_number, col="syns_en_source"
            )
            # PMID tokens become canonical PubMed URLs. Existing non-PMID
            # source strings are preserved, so the reification/provenance model
            # itself is unchanged.
            for source in (normalized or [src]):
                syn_en_map[idx].append((se, source))

    syn_fr_map = defaultdict(list)  # dict: {aab_id: [(syn_fr, source), ...]}
    for row_number, r in enumerate(syn_fr, start=2):
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        sf = (r.get("syns_fr") or "").strip()
        src = (r.get("syns_fr_source") or "").strip()
        if idx and sf:
            normalized = to_pubmed_urls(
                src, src_file=syn_fr_path, row=row_number, col="syns_fr_source"
            )
            for source in (normalized or [src]):
                syn_fr_map[idx].append((sf, source))

    hpo_by_idx = defaultdict(list)  # dict: {aab_id: [hpo_id, ...]}
    for r in hpo_id:
        idx = (r.get("index") or "").strip()
        if idx not in KEEP:
            continue
        hp_raw = (r.get("hpo_id") or r.get("hpo") or "").strip()
        if not (idx and hp_raw):
            continue
        for hp in split_hpo_values(hp_raw):
            if hp not in hpo_by_idx[idx]:
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
        did_raw = (r.get("related_diseases_id") or "").strip()
        pm = to_pubmed_urls(
            r.get("diseases_pmids"), src_file=dis_path, row=i, col="diseases_pmids"
        )
        if not (idx and did_raw):
            continue
        # A cell may contain one or several disease IDs (e.g., separated by | ; ,)
        for did in [x.strip() for x in re.split(r"[|;,\n\r]+", did_raw) if x.strip()]:
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
    up_names,
    up_urls,
    umls_cn_names,
    umls_cn_urls,
    hpo_cn_names,
    chebi_cn_names,
    chebi_cn_urls,
    umls_synonyms=None,
):
    # External root classes are retained as domain classes. Target is assigned
    # to the selected local individuals, not globally to every protein/molecule.
    uniprot_protein = UP.Protein

    # Declare makaao_core CSV as a Document used by the Activity
    csv_doc = MAKAAO["makaao_core.csv"] # add something with Core version here?
    g.add((csv_doc, RDF.type, MAKAAO.Document))
    g.add((csv_doc, RDFS.seeAlso, URIRef("https://makaao.inria.fr/data/makaao_core.csv"))) # add something with Core version here?
    add_pref(g, csv_doc, os.path.basename(makaao_core_name))
    g.add((GLOBAL_ACTIVITY, PROV.used, csv_doc))

    hpo_local_names = hpo_cn_names or {}
    umls_synonyms = umls_synonyms or {}
    primary_labels, aab_class_uri = {}, {}
    pos_uris_by_idx = {}
    positivity_instances_by_hpo = defaultdict(set)

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

        # Every autoantibody has exactly one local positivity class and one
        # local positivity individual. HPO correspondences are non-logical
        # closeMatch links on that class; no HPO identifier enters a local URI.
        local_pos_label = primary_labels[idx] + " positivity"
        pos_inst = MAKAAO[f"positivity_{idx}_instance"]
        mapped_hpo_codes = []
        if idx == "18":
            structural_pos_cls = MAKAAO.AutoantibodyPositivity
            root_code = "HP:0030057"
            root_label = (
                (hpo_local_names.get(root_code) or "").strip()
                or "Autoimmune antibody positivity"
            )
            add_pref(g, structural_pos_cls, root_label)
            add_pref(g, HPO_AUTOIMMUNE_ANTIBODY_POSITIVITY, root_label)
            positivity_instance_label = root_label
            mapped_hpo_codes.append(root_code)
        else:
            structural_pos_cls = MAKAAO[f"positivity_{idx}"]
            g.add((structural_pos_cls, RDF.type, OWL.Class))
            add_pref(g, structural_pos_cls, local_pos_label)
            positivity_instance_label = local_pos_label

            # Link the one AAb-specific positivity class directly to every
            # corresponding HPO class. Materialize both directions because
            # skos:closeMatch is symmetric and the canonical KG is checked
            # without relying on OWL inference. HPO labels remain on HPO resources.
            seen_hpo_codes = set()
            for hp_code in data["hpo_list"].get(idx, []):
                code_norm = canonical_hpo_code(hp_code)
                hpo_cls = hp_to_obo_uri(code_norm)
                if not hpo_cls or not code_norm or code_norm in seen_hpo_codes:
                    continue
                seen_hpo_codes.add(code_norm)
                hpo_label = (hpo_local_names.get(code_norm) or "").strip()
                g.add((structural_pos_cls, SKOS.closeMatch, hpo_cls))
                g.add((hpo_cls, SKOS.closeMatch, structural_pos_cls))
                add_pref(g, hpo_cls, hpo_label or code_norm)
                mapped_hpo_codes.append(code_norm)

        # The structural class alone participates in the mirrored taxonomy.
        pos_uris_by_idx[idx] = [structural_pos_cls]

        # The individual instantiates only the one AAb-specific local class.
        add_pref(g, pos_inst, positivity_instance_label)
        g.add((pos_inst, RDF.type, structural_pos_cls))
        for hpo_code in mapped_hpo_codes:
            positivity_instances_by_hpo[hpo_code].add(pos_inst)
        # Materialize both directions of the Biolink inverse pair so the
        # canonical KG can be validated without requiring a reasoner.
        g.add((inst, BIOLINK.biomarker_for, pos_inst))
        g.add((pos_inst, BIOLINK.has_biomarker, inst))

    # Mirror only the autoantibody taxonomy in the one-per-AAb local positivity
    # classes. HPO closeMatch mappings do not participate in this hierarchy.
    for idx in sorted(data["indices"], key=lambda value: int(value)):
        if idx == "18":
            continue

        child_uri = pos_uris_by_idx[idx][0]
        parents = data["parents"].get(idx) or []
        valid_parents = sorted(
            {p for p in parents if p != idx and p in data["indices"]},
            key=lambda value: int(value),
        )
        parent_list = (
            [pos_uris_by_idx[parent_idx][0] for parent_idx in valid_parents]
            if valid_parents
            else [MAKAAO.AutoantibodyPositivity]
        )

        for parent_uri in parent_list:
            if child_uri != parent_uri:
                g.add((child_uri, RDFS.subClassOf, parent_uri))


    # UMLS targets
    umls_local_names = umls_cn_names or {}  # umls_names from code_names table
    for idx, items in data["cui"].items():
        inst_aab = MAKAAO[
            f"aab_{idx}_instance"
        ]  # for each aab_id, we get the name of the corresponding instance
        for cui, pmids in items:
            cui_key = canonical_cui_code(cui)
            if not cui_key:
                continue
            preferred = umls_local_names.get(cui_key)
            source_url = (umls_cn_urls or {}).get(cui_key)
            _, cui_uri = ensure_cui_class_and_instance(
                g,
                cui_key,
                preferred_label=preferred,
                source_url=source_url,
                synonyms=umls_synonyms.get(cui_key, ()),
            )
            g.add((cui_uri, RDF.type, MAKAAO.Target))
            g.add((inst_aab, BAO_HAS_TARGET, cui_uri))
            g.add((cui_uri, BAO_IS_TARGET_FOR, inst_aab))
            for p in pmids or [""]:
                add_reified_relation(g, inst_aab, BAO_HAS_TARGET, cui_uri, p)

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
            prot_cls = MAKAAO[f"UP_{make_valid(base)}"]
            prot_ind = MAKAAO[f"UP_{make_valid(base)}_instance"]
            g.add((prot_cls, RDF.type, OWL.Class))
            g.add((prot_cls, RDFS.subClassOf, uniprot_protein))
            g.add((prot_ind, RDF.type, prot_cls))
            g.add((prot_ind, RDF.type, MAKAAO.Target))
            up_name = up_names.get(base) or up_names.get(norm)
            if up_name:
                add_pref(g, prot_cls, up_name)
                add_pref(g, prot_ind, up_name)
                up_names_added += 1
            else:
                add_pref(g, prot_cls, base)
                add_pref(g, prot_ind, base)
            g.add((prot_cls, SKOS.notation, Literal(base)))
            up_url = up_urls.get(base) or up_urls.get(norm)
            external_up = normalize_http_iri(up_url, f"UniProt URL for {base}") if up_url else UNIPROT[base]
            g.add((prot_cls, RDFS.seeAlso, external_up))
            g.add((prot_ind, RDFS.seeAlso, external_up))
            g.add((inst_aab, BAO_HAS_TARGET, prot_ind))
            g.add((prot_ind, BAO_IS_TARGET_FOR, inst_aab))
            for p in pmids or [""]:
                add_reified_relation(g, inst_aab, BAO_HAS_TARGET, prot_ind, p)
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
            # Expose the used external ChEBI class in the canonical MAKAAO KG.
            # Its genuine ChEBI hierarchy remains supplied by the imported module;
            # ``Target`` remains a contextual role on the local individual.
            g.add((chebi_cls, RDF.type, OWL.Class))
            g.add((chebi_ind, RDF.type, chebi_cls))
            g.add((chebi_ind, RDF.type, MAKAAO.Target))

            chebi_name = (chebi_cn_names or {}).get(code_colon)
            if chebi_name:
                add_pref(g, chebi_cls, chebi_name)
                add_pref(g, chebi_ind, chebi_name)
            else:
                add_pref(g, chebi_cls, code_colon)
                add_pref(g, chebi_ind, code_colon)
            g.add((chebi_cls, SKOS.notation, Literal(code_colon)))
            g.add((chebi_ind, SKOS.notation, Literal(code_colon)))

            chebi_url = (chebi_cn_urls or {}).get(code_colon)
            if chebi_url:
                external_chebi = normalize_http_iri(
                    chebi_url, f"ChEBI URL for {code_colon}"
                )
                g.add((chebi_cls, RDFS.seeAlso, external_chebi))
                g.add((chebi_ind, RDFS.seeAlso, external_chebi))
            g.add((inst_aab, BAO_HAS_TARGET, chebi_ind))
            g.add((chebi_ind, BAO_IS_TARGET_FOR, inst_aab))
            for p in pmids or [""]:
                add_reified_relation(
                    g, inst_aab, BAO_HAS_TARGET, chebi_ind, p
                )  # for each source we have, we add a reified relation to carry provenance

    return {
        hpo_code: tuple(sorted(instances, key=str))
        for hpo_code, instances in positivity_instances_by_hpo.items()
    }

def process_diseases(
    g,
    data,
    orpha_hpo_links,
    positivity_instances_by_hpo,
    umls_cn_names,
    umls_cn_urls=None,
    hpo_cn_names=None,
    umls_synonyms=None,
):
    hpo_cn_names = hpo_cn_names or {}
    umls_synonyms = umls_synonyms or {}

    orpha_cn_names, _ = (
        read_code_names_orpha(CODE_NAMES_CSV)
        if os.path.exists(CODE_NAMES_CSV)
        else ({}, {})
    )

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
            cui_norm = canonical_cui_code(code_upper)

            # ORPHANET
            orpha_num = canonical_orpha_code(code_raw)
            if orpha_num:
                d_cls = URIRef(f"http://www.orpha.net/ORDO/Orphanet_{orpha_num}")
                # The ORDO declaration and source hierarchy come from the
                # imported ORDO module. This is the MAKAAO alignment axiom.
                g.add((d_cls, RDFS.subClassOf, MAKAAO.AutoimmuneDisease))
                inst = MAKAAO[
                    f"orpha_{orpha_num}_instance"
                ]  # also add instance of Orphanet class
                g.add((inst, RDF.type, d_cls))
                orpha_name = (orpha_cn_names.get(orpha_num) or "").strip()
                if orpha_name:
                    add_pref(g, inst, orpha_name)

                # Use the ORDO URI string (same as keys in orpha_hpo_links)
                for link in orpha_hpo_links.get(str(d_cls), []):
                    # HPO identifier (may be HP:..., HP_..., or full URI)
                    hpo_id_raw = (link.get("HPOId") or link.get("hpoid") or "").strip()
                    if not hpo_id_raw:
                        continue

                    # Normalize HPO code and build URI
                    code_norm = canonical_hpo_code(hpo_id_raw)
                    pos_uri = hp_to_obo_uri(code_norm)

                    if pos_uri is None or not code_norm:
                        continue

                    # Label priority: code_names.csv HPO first, then Orphadata HPOTerm
                    term = (
                        (hpo_cn_names.get(code_norm) or "").strip()
                        or (link.get("HPOTerm") or link.get("hpoterm") or "").strip()
                    )
                    # Reuse local positivity individuals whenever this HPO
                    # class corresponds to a generated positivity class. Only
                    # otherwise instantiate the external HPO class directly.
                    phenotype_instances = positivity_instances_by_hpo.get(
                        code_norm, ()
                    )
                    if not phenotype_instances:
                        pos_inst = MAKAAO[f"hpo_{make_valid(code_norm)}_instance"]
                        g.add((pos_inst, RDF.type, pos_uri))
                        if term:
                            add_pref(g, pos_inst, term)
                        class_label = g.value(pos_uri, RDFS.label) or g.value(
                            pos_uri, SKOS.prefLabel
                        )
                        if class_label:
                            g.add((pos_inst, RDFS.label, class_label))
                        phenotype_instances = (pos_inst,)

                    for pos_inst in phenotype_instances:
                        g.add((inst, SIO["SIO_001279"], pos_inst))  # has_phenotype
                        g.add((pos_inst, SIO["SIO_001280"], inst))  # is_phenotype_of
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

            # UMLS CUI disease
            elif cui_norm:  # if current disease code contains a CUI:
                code_norm = cui_norm
                preferred = (umls_cn_names or {}).get(code_norm)
                source_url = (umls_cn_urls or {}).get(code_norm)
                cui_cls, inst = ensure_cui_class_and_instance(
                    g,
                    code_norm,
                    preferred_label=preferred,
                    source_url=source_url,
                    synonyms=umls_synonyms.get(code_norm, ()),
                )
                # Harmonize with ORDO disease modelling: classify the concept
                # class, not the individual, as an autoimmune disease.
                g.add((cui_cls, RDFS.subClassOf, MAKAAO.AutoimmuneDisease))

            # Last case: if the disease value is neither an ORPHA identifier nor
            # a CUI, create only a generic local individual. Since no dedicated
            # disease class exists for this fallback, direct typing is retained.
            else:
                code_norm = code_upper.replace("CUI:", "").strip()
                fragment = safe_local_fragment(code_norm, prefix="disease")
                inst = MAKAAO[f"{fragment}_instance"]
                add_pref(g, inst, code_upper)
                g.add((inst, RDF.type, MAKAAO.AutoimmuneDisease))

            if inst is None:
                continue
            g.add(
                (aab_inst, SIO["SIO_001403"], inst)
            )  # we link that disease instance to the current aab instance
            g.add((inst, SIO["SIO_001403"], aab_inst))
            for p in prov or [""]:
                add_reified_relation(
                    g, aab_inst, SIO["SIO_001403"], inst, p
                )  # for each provenance entry, we add a reified relation to the graph, for the current instance

def read_required_json_object(path, description):
    """Read a required JSON object and fail on absence or malformed content."""
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(f"Required {description} not found: {path}")
    try:
        with open(path, "r", encoding="utf-8") as handle:
            value = json.load(handle)
    except OSError as exc:
        raise OSError(f"Could not read required {description}: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"Malformed JSON in required {description} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Required {description} {path} must contain a JSON object")
    return value


def process_loinc_mappings(
    g,
    loinc_index_csv,
    keep_indices,
    part_labels,
    term_labels,
    part_test_json=LOINC_PART_TEST_JSON,
):
    """Create relevant LOINC resources using labels from ``code_names.csv``."""
    map_rows = read_csv_rows(loinc_index_csv)
    if not map_rows:
        return

    raw_part_tests = read_required_json_object(
        part_test_json, "LOINC Part-to-Term dictionary"
    )
    part_test_sets = defaultdict(set)
    for raw_part, raw_tests in raw_part_tests.items():
        part_code = canonical_loinc_code(raw_part, "part")
        if not part_code:
            raise ValueError(
                f"{part_test_json} contains an invalid LOINC Part key: {raw_part!r}"
            )
        if not isinstance(raw_tests, list):
            raise ValueError(
                f"{part_test_json} value for {raw_part!r} must be a JSON array"
            )

        # Preserve valid Parts whose Primary-link list is empty. Accessing the
        # defaultdict here creates the canonical key even when the loop below
        # has no Terms to add.
        part_test_sets[part_code]

        for raw_test in raw_tests:
            test_code = canonical_loinc_code(raw_test, "term")
            if not test_code:
                raise ValueError(
                    f"{part_test_json} contains an invalid LOINC Term for "
                    f"{part_code}: {raw_test!r}"
                )
            part_test_sets[part_code].add(test_code)
    part_tests = {
        part_code: sorted(test_codes, key=lambda value: (value.casefold(), value))
        for part_code, test_codes in part_test_sets.items()
    }

    relevant_rows = [
        row
        for row in map_rows
        if (row.get("aab_id") or "").strip() in keep_indices
    ]
    used_parts = {
        code
        for row in relevant_rows
        if (code := canonical_loinc_code(row.get("loinc_id"), "part"))
    }
    missing_part_mappings = sorted(used_parts - set(part_tests))
    missing_part_labels = sorted(used_parts - set(part_labels))
    used_terms = {
        term
        for part in used_parts
        for term in part_tests.get(part, [])
    }
    missing_term_labels = sorted(used_terms - set(term_labels))
    if missing_part_mappings or missing_part_labels or missing_term_labels:
        details = []
        if missing_part_mappings:
            details.append(
                "Parts absent from Part-to-Term dictionary: "
                + ", ".join(missing_part_mappings[:20])
            )
        if missing_part_labels:
            details.append(
                "Parts absent from label dictionary: "
                + ", ".join(missing_part_labels[:20])
            )
        if missing_term_labels:
            details.append(
                "Terms absent from label dictionary: "
                + ", ".join(missing_term_labels[:20])
            )
        raise ValueError(
            "LOINC enrichment files are incomplete for the processed LOINC input:\n  "
            + "\n  ".join(details)
        )

    # Small local schema used to organize the imported LOINC classes. Reuse
    # CompLOINC's COMPONENT predicate instead of minting a MAKAAO predicate.
    g.add((MAKAAO.LoincPart, RDF.type, OWL.Class))
    g.add((MAKAAO.LoincPart, RDFS.label, Literal("LOINC Part", lang="en")))
    g.add((MAKAAO.LoincTerm, RDF.type, OWL.Class))
    g.add((MAKAAO.LoincTerm, RDFS.label, Literal("LOINC Term", lang="en")))
    g.add((LOINC_COMPONENT, RDF.type, OWL.ObjectProperty))
    g.add((LOINC_COMPONENT, RDFS.label, Literal("has LOINC component", lang="en")))
    g.add((LOINC_COMPONENT, RDFS.domain, MAKAAO.LoincTerm))
    g.add((LOINC_COMPONENT, RDFS.range, MAKAAO.LoincPart))

    seen_parts = set()
    seen_terms = set()

    for r in relevant_rows:
        idx = (r.get("aab_id") or "").strip()
        part_code = canonical_loinc_code(r.get("loinc_id"), "part")
        if not part_code:
            continue

        aab_inst = MAKAAO[f"aab_{idx}_instance"]
        part_cls = MAKAAO_LOINC[part_code]
        part_inst = MAKAAO[f"loinc_{make_valid(part_code)}_instance"]
        part_label = part_labels.get(part_code) or part_code

        if part_code not in seen_parts:
            seen_parts.add(part_code)
            g.add((part_cls, RDF.type, OWL.Class))
            g.add((part_cls, RDFS.subClassOf, MAKAAO.LoincPart))
            add_pref(g, part_cls, part_label)
            g.add((part_cls, SKOS.notation, Literal(part_code)))
            g.add((part_cls, RDFS.seeAlso, LOINC[part_code]))

            # Instantiate the imported LOINC Part class locally.
            g.add((part_inst, RDF.type, part_cls))
            g.add((part_inst, RDF.type, SKOS.Concept))
            add_pref(g, part_inst, part_label)

        # Map the existing MAKAAO Autoantibody concept individual to the Part concept individual.
        g.add((aab_inst, RDF.type, SKOS.Concept))
        g.add((aab_inst, SKOS.closeMatch, part_inst))
        g.add((part_inst, SKOS.closeMatch, aab_inst))

        # Import only terms associated with a Part used by MAKAAO.
        for test_code in part_tests.get(part_code, []):
            term_cls = MAKAAO_LOINC[test_code]
            term_inst = MAKAAO[f"loinc_{make_valid(test_code)}_instance"]
            term_label = term_labels.get(test_code) or test_code

            if test_code not in seen_terms:
                seen_terms.add(test_code)
                g.add((term_cls, RDF.type, OWL.Class))
                g.add((term_cls, RDFS.subClassOf, MAKAAO.LoincTerm))
                add_pref(g, term_cls, term_label)
                g.add((term_cls, SKOS.notation, Literal(test_code)))
                g.add((term_cls, RDFS.seeAlso, LOINC[test_code]))

                # Instantiate the imported LOINC Term class locally.
                g.add((term_inst, RDF.type, term_cls))
                g.add((term_inst, RDF.type, SKOS.Concept))
                add_pref(g, term_inst, term_label)

            # Relate the LOINC Term individual to the LOINC Part individual.
            g.add((term_inst, LOINC_COMPONENT, part_inst))

    print(
        f"Loaded LOINC classes and instances — Parts:{len(seen_parts)} Terms:{len(seen_terms)}"
    )

def set_output_file_metadata(graph: Graph, output_name: str) -> URIRef:
    """Set file-specific DCAT distribution metadata and an exact VoID count.

    Derived graphs must not retain the full KG's download URL or triple count.
    Existing distributions linked from the MAKAAO dataset are removed before a
    distribution for ``output_name`` is created. The VoID count includes its own
    ``void:triples`` assertion.
    """
    dataset = URIRef("http://makaao.inria.fr/kg/")
    output_filename = Path(output_name).name
    if not output_filename:
        raise ValueError("Output filename is empty")

    old_distributions = set(graph.objects(dataset, DCAT.distribution))
    graph.remove((dataset, DCAT.distribution, None))
    graph.remove((dataset, VOID.triples, None))
    for distribution in old_distributions:
        graph.remove((distribution, None, None))
        graph.remove((None, None, distribution))

    distribution_token = safe_local_fragment(
        Path(output_filename).stem, prefix="distribution", max_slug_length=80
    )
    distribution = URIRef(
        f"http://makaao.inria.fr/kg/distribution/{distribution_token}"
    )
    graph.add((dataset, DCAT.distribution, distribution))
    graph.add((distribution, RDF.type, DCAT.Distribution))
    graph.add(
        (
            distribution,
            DCAT.downloadURL,
            URIRef(f"https://makaao.inria.fr/data/{quote(output_filename, safe='._-')}"),
        )
    )
    graph.add(
        (
            distribution,
            DCAT.mediaType,
            URIRef("https://www.iana.org/assignments/media-types/application/rdf+xml"),
        )
    )
    graph.add(
        (
            distribution,
            DCTERMS.license,
            URIRef("https://creativecommons.org/licenses/by/4.0/"),
        )
    )

    graph.add(
        (
            dataset,
            VOID.triples,
            Literal(len(graph) + 1, datatype=XSD.integer),
        )
    )
    recorded_counts = list(graph.objects(dataset, VOID.triples))
    if len(recorded_counts) != 1 or recorded_counts[0].toPython() != len(graph):
        raise RuntimeError(
            f"Could not set exact void:triples metadata for {output_filename}"
        )
    return distribution


# ===================== FAIR / DCAT METADATA =====================
# function to add some hard-coded triples to the KG
def append_fair_metadata(kg: Graph):
    ONT = URIRef("http://makaao.inria.fr/kg/")
    kg.add((ONT, RDF.type, OWL.Ontology))
    kg.add((ONT, RDF.type, DCAT.Dataset))
    kg.add((ONT, RDF.type, VOID.Dataset))

    kg.add((ONT, DCTERMS.identifier, Literal(f"MAKAAO-{version}")))
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

    OPEN_ACCESS = URIRef("http://makaao.inria.fr/kg/access-right/open")
    kg.add((OPEN_ACCESS, RDF.type, DCTERMS.RightsStatement))
    kg.add(
        (
            OPEN_ACCESS,
            RDFS.label,
            Literal("Open access", lang="en"),
        )
    )
    kg.add((ONT, DCTERMS.accessRights, OPEN_ACCESS))

    POLICY_PAGE = URIRef("https://makaao.inria.fr/usage_policies.html")
    kg.add((POLICY_PAGE, RDF.type, DCTERMS.RightsStatement))
    kg.add((ONT, DCTERMS.rights, POLICY_PAGE))


    SERVICE = URIRef("http://makaao.inria.fr/kg/service")
    kg.add((SERVICE, RDF.type, DCAT.DataService))
    kg.add(
        (
            SERVICE,
            DCAT.endpointURL,
            URIRef("http://makaao.inria.fr/kg/"),
        )
    )

    
    kg.add((SERVICE, DCAT.servesDataset, ONT))
    kg.add((ONT, DCAT.landingPage, URIRef("https://makaao.inria.fr/")))

    ACT = URIRef("http://makaao.inria.fr/kg/activity_makaao_core")
    kg.add((ONT, PROV.wasGeneratedBy, ACT))
    kg.add((ACT, RDF.type, PROV.Activity))
    #kg.add((ONT, PROV.wasDerivedFrom, ACT))

    author_uri = URIRef("https://heka.gitlabpages.inria.fr/team/members/maury.html")
    kg.add((ONT, DCTERMS.creator, author_uri))
    kg.add((author_uri, RDF.type, FOAF.Person))
    add_pref(kg, author_uri, "Fabien Maury")

    team_alpha = URIRef("https://heka.gitlabpages.inria.fr/")
    team_beta = URIRef("https://www.institutimagine.org/en")
    for team_uri, team_label in [
        (team_alpha, "Team HeKA"),
        (team_beta, "Institut Imagine, Inserm"),
    ]:
        kg.add((ONT, DCTERMS.contributor, team_uri))
        kg.add((team_uri, RDF.type, FOAF.Organization))
        add_pref(kg, team_uri, team_label)

    kg.add(
        (
            ONT,
            DCTERMS.created,
            Literal("2024-01-15", datatype=XSD.date),
        )
    )

    kg.add(
        (
            ONT,
            DCTERMS.modified,
            Literal(date.today().isoformat(), datatype=XSD.date),
        )
    )
    kg.add((ONT, OWL.versionInfo, Literal(version)))
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

    from rdflib.namespace import Namespace

    VCARD = Namespace("http://www.w3.org/2006/vcard/ns#")

    CONTACT = URIRef("http://makaao.inria.fr/kg/contact.html")

    kg.add((ONT, DCAT.contactPoint, CONTACT))
    kg.add((CONTACT, RDF.type, VCARD.Organization))
    kg.add(
        (
            CONTACT,
            VCARD.hasEmail,
            URIRef("mailto:contact.makaao@inria.fr"),
        )
    )
# ===================== T-BOX EXPORT =====================

_TBOX_ONTOLOGY_HEADER_PREDICATES = {
    RDF.type,
    RDFS.label,
    RDFS.comment,
    SKOS.prefLabel,
    OWL.versionIRI,
    OWL.versionInfo,
    OWL.imports,
}


def validate_tbox_export(tbox: Graph) -> None:
    """Ensure the standalone ontology export contains no dataset-service ABox."""
    forbidden_types = {
        DCAT.Dataset,
        DCAT.Distribution,
        DCAT.DataService,
        VOID.Dataset,
    }
    forbidden_predicates = {
        DCAT.distribution,
        DCAT.downloadURL,
        DCAT.endpointURL,
        DCAT.servesDataset,
        DCAT.contactPoint,
        VOID.triples,
    }
    typed_metadata_individuals = {
        subject
        for subject, obj in tbox.subject_objects(RDF.type)
        if obj in forbidden_types
    }
    forbidden_assertions = [
        (subject, predicate, obj)
        for predicate in forbidden_predicates
        for subject, obj in tbox.subject_objects(predicate)
    ]
    if typed_metadata_individuals or forbidden_assertions:
        raise RuntimeError(
            "Standalone MAKAAO TBox contains dataset/distribution/service ABox metadata"
        )

def extract_tbox(source: Graph, local_ns: str) -> Graph:
    """
    Extract the ontology T-box while retaining the generated LOINC integration.

    Included:
      - schema terms defined in the MAKAAO namespace;
      - generated LOINC Part and Term classes that are direct subclasses of
        makaao:LoincPart or makaao:LoincTerm;
      - the reused loinc_property:COMPONENT object-property schema;
      - local positivity classes and their direct HPO close-match mappings;
      - recursively reachable blank-node structures, such as OWL restrictions.

    Excluded:
      - all individuals and individual-level A-box assertions. Class-level
        skos:closeMatch alignments and generated LOINC schema links are retained.
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

    # Retain the original local T-box.
    for s, o in source.subject_objects(RDF.type):
        if o in schema_types and isinstance(s, URIRef) and str(s).startswith(local_ns):
            seeds.add(s)

    # Retain only the LOINC classes generated by process_loinc_mappings().
    seeds.update(source.subjects(RDFS.subClassOf, MAKAAO.LoincPart))
    seeds.update(source.subjects(RDFS.subClassOf, MAKAAO.LoincTerm))

    # Retain project alignment axioms whose subject is an external class
    # and whose superclass is a MAKAAO class (for example ORDO disease classes).
    for external_cls, local_parent in source.subject_objects(RDFS.subClassOf):
        if (
            isinstance(external_cls, URIRef)
            and not str(external_cls).startswith(local_ns)
            and isinstance(local_parent, URIRef)
            and str(local_parent).startswith(local_ns)
        ):
            seeds.add(external_cls)

    # Retain used ChEBI classes explicitly exposed by build_core(). Their source
    # hierarchy is still supplied by the imported ChEBI reasoning module.
    for external_cls in source.subjects(RDF.type, OWL.Class):
        if _is_chebi_class(external_cls):
            seeds.add(external_cls)

    # Retain the property used to relate LOINC Term and Part individuals.
    if (LOINC_COMPONENT, RDF.type, OWL.ObjectProperty) in source:
        seeds.add(LOINC_COMPONENT)

    ontology_headers = {
        subject
        for subject in source.subjects(RDF.type, OWL.Ontology)
        if isinstance(subject, URIRef) and str(subject).startswith(local_ns)
    }

    # Copy all triples about schema seeds and recursively copy blank nodes.
    # For the ontology header, retain only ontology metadata and imports. The
    # DCAT/VoID publication records belong to the canonical KG, not the TBox.
    stack = list(seeds)
    seen = set(seeds)
    while stack:
        s = stack.pop()
        for _, p, o in source.triples((s, None, None)):
            if s in ontology_headers:
                if p not in _TBOX_ONTOLOGY_HEADER_PREDICATES:
                    continue
                if p == RDF.type and o != OWL.Ontology:
                    continue
            tbox.add((s, p, o))
            if isinstance(o, BNode) and o not in seen:
                seen.add(o)
                stack.append(o)

    validate_tbox_export(tbox)
    return tbox

def build_strict_reasoning_release(
    kg_graph: Graph,
    tbox_graph: Graph,
    non_reified_graph: Graph,
    output_dir: Path,
    canonical_kg_path: Path,
    canonical_tbox_path: Path,
    subprocess_env: dict[str, str],
):
    """Build strict reasoning modules and reasoned outputs from pinned sources."""
    if not BUILD_STRICT_REASONING_RELEASE:
        return {"status": "DISABLED"}

    builder_path = Path(REASONING_RELEASE_BUILDER)
    if not builder_path.is_file():
        raise FileNotFoundError(
            f"Integrated reasoning-release helper not found: {builder_path}"
        )

    module_name = f"makaao_reasoning_release_builder_{SCRIPT_VERSION.replace('.', '_')}"
    spec = importlib.util.spec_from_file_location(module_name, builder_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load reasoning helper: {builder_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        helper_version = getattr(module, "SCRIPT_VERSION", None)
        if helper_version != SCRIPT_VERSION:
            raise RuntimeError(
                "Graph builder and reasoning helper versions differ: "
                f"graph={SCRIPT_VERSION!r}, helper={helper_version!r}"
            )
        result = module.build_reasoning_release(
            kg_graph=kg_graph,
            tbox_graph=tbox_graph,
            non_reified_graph=non_reified_graph,
            output_dir=output_dir,
            source_files=dict(EXTERNAL_ONTOLOGY_FILES),
            skos_source_file=SKOS_SOURCE_FILE,
            kg_version=version,
            robot_command=_robot_command(),
            reasoner_name="HermiT",
            timeout_seconds=REASONER_TIMEOUT_SECONDS,
            require_robot=REQUIRE_STRICT_REASONING_RELEASE,
            canonical_kg_path=canonical_kg_path,
            canonical_tbox_path=canonical_tbox_path,
            subprocess_env=subprocess_env,
        )
    finally:
        sys.modules.pop(spec.name, None)

    print(f"Strict reasoning release: {result['status']}")
    return result


def _serialize_rdfxml_verified(graph: Graph, path: Path, context: str) -> None:
    """Serialize a canonical root graph and verify it can be parsed unchanged."""
    path.parent.mkdir(parents=True, exist_ok=True)
    graph.serialize(destination=str(path), format="xml")
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"RDF/XML serialization did not create {context}: {path}")
    verified = Graph()
    try:
        verified.parse(str(path), format="xml")
        if len(verified) != len(graph):
            raise RuntimeError(
                f"RDF/XML round-trip changed triple count for {context}: "
                f"{len(graph)} != {len(verified)}"
            )
    finally:
        verified.close()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_root_reasoning_catalog(stage_root: Path, reasoning_result: dict) -> Path:
    """Write and validate the KG-level XML catalog for the strict module set."""
    from xml.etree import ElementTree as ET

    if REQUIRE_STRICT_REASONING_RELEASE and reasoning_result.get("status") != "PASSED_OWL2_DL_AND_HERMIT":
        raise RuntimeError(
            "Strict reasoning release did not complete successfully: "
            f"{reasoning_result.get('status')!r}"
        )

    module_iris = {
        "biolink": BIOLINK_MODULE_IRI,
        "hpo": HPO_MODULE_IRI,
        "chebi": CHEBI_MODULE_IRI,
        "ordo": ORDO_MODULE_IRI,
        "sio": SIO_MODULE_IRI,
        "prov": PROV_MODULE_IRI,
        "bao": BAO_MODULE_IRI,
        "uniprot": UNIPROT_MODULE_IRI,
        "skos": SKOS_OWL_DL,
    }
    modules = list(reasoning_result.get("modules", []))
    module_names = [module.get("name") for module in modules]
    if len(module_names) != len(set(module_names)) or set(module_names) != set(module_iris):
        raise RuntimeError(
            "Reasoning result does not contain exactly one module for each expected vocabulary: "
            f"expected={sorted(module_iris)}, actual={sorted(map(str, module_names))}"
        )

    catalog_path = stage_root / "catalog-imports.xml"
    root = ET.Element(
        "catalog",
        {
            "xmlns": "urn:oasis:names:tc:entity:xmlns:xml:catalog",
            "prefer": "public",
        },
    )
    for module in sorted(modules, key=lambda item: item["name"]):
        name = module["name"]
        filename = Path(module["path"]).name
        module_path = stage_root / REASONING_OUTPUT_DIRNAME / "modules" / filename
        if not module_path.is_file() or module_path.stat().st_size == 0:
            raise RuntimeError(f"Reasoning module missing before catalog publication: {module_path}")
        ET.SubElement(
            root,
            "uri",
            {
                "name": str(module_iris[name]),
                "uri": f"{REASONING_OUTPUT_DIRNAME}/modules/{filename}",
            },
        )
    ET.ElementTree(root).write(catalog_path, encoding="UTF-8", xml_declaration=True)

    # Parse the generated XML and verify every relative target while still in staging.
    parsed_root = ET.parse(catalog_path).getroot()
    entries = [element for element in parsed_root if element.tag.endswith("uri")]
    if len(entries) != len(module_iris):
        raise RuntimeError(f"Root catalog contains {len(entries)} entries; expected {len(module_iris)}")
    for entry in entries:
        target = stage_root / entry.attrib["uri"]
        if not target.is_file():
            raise RuntimeError(f"Root catalog target does not exist: {target}")
    return catalog_path


def finalize_release_manifest_and_checksums(
    stage_root: Path,
    stage_reasoning: Path,
    stage_kg: Path,
    stage_tbox: Path,
    root_catalog: Path,
) -> None:
    """Bind every permanent root and reasoning artifact into one checksum set."""
    manifest_path = stage_reasoning / "reasoning-manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"Reasoning manifest missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    canonical_inputs = manifest.setdefault("canonical_inputs", {})

    expected_root_inputs = {
        "kg": stage_kg,
        "tbox": stage_tbox,
    }
    for key, path in expected_root_inputs.items():
        recorded = canonical_inputs.get(key)
        actual_hash = _file_sha256(path)
        if not isinstance(recorded, dict) or recorded.get("sha256") != actual_hash:
            raise RuntimeError(
                f"Manifest hash for canonical {key} does not match staged file {path}"
            )

    canonical_inputs["root_catalog"] = {
        "path": f"../{root_catalog.name}",
        "sha256": _file_sha256(root_catalog),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    checksum_path = stage_reasoning / "SHA256SUMS"
    entries: list[tuple[str, Path]] = []
    for path in stage_reasoning.rglob("*"):
        if path.is_file() and path != checksum_path:
            entries.append((path.relative_to(stage_reasoning).as_posix(), path))
    for path in (root_catalog, stage_kg, stage_tbox):
        entries.append((f"../{path.name}", path))

    names = [name for name, _ in entries]
    if len(names) != len(set(names)):
        raise RuntimeError("Duplicate paths detected while generating SHA256SUMS")
    with checksum_path.open("w", encoding="utf-8") as stream:
        for relative_name, path in sorted(entries):
            stream.write(f"{_file_sha256(path)}  {relative_name}\n")

    # Verify the concrete staged release immediately, using the same relative
    # paths that users will verify after publication.
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        expected_hash, relative_name = line.split("  ", 1)
        target = stage_reasoning / relative_name
        if not target.is_file() or _file_sha256(target) != expected_hash:
            raise RuntimeError(f"Checksum verification failed before publication: {relative_name}")


def _robot_command():
    """Return an explicit ROBOT command, preferring a pinned local JAR."""
    if not re.fullmatch(r"[1-9][0-9]*[mMgG]", JAVA_MAX_HEAP):
        raise ValueError(
            "MAKAAO_JAVA_MAX_HEAP must look like 8G or 4096M, "
            f"got {JAVA_MAX_HEAP!r}"
        )

    candidates = [
        Path(ROBOT_JAR),
        PROJECT_DIR / "kg" / "robot.jar",
        KG_DIR / "robot.jar",
        PROJECT_DIR / "robot.jar",
    ]
    seen = set()
    for jar in candidates:
        resolved = jar.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if jar.is_file():
            java = shutil.which("java")
            if java is None:
                raise RuntimeError(
                    f"Found ROBOT JAR at {jar}, but Java is not available on PATH"
                )
            return [java, f"-Xmx{JAVA_MAX_HEAP}", "-jar", str(jar)]

    executable = shutil.which(ROBOT_EXECUTABLE)
    if executable:
        return [executable]
    return None


def validate_build_prerequisites() -> None:
    """Fail before graph construction when required build inputs are unavailable."""
    filesystem_root = Path(KG_DIR.anchor).resolve()
    project_root = PROJECT_DIR.resolve()
    data_root = DATA_DIR.resolve()
    script_root = SCRIPT_DIR.resolve()

    protected_trees = {data_root}
    if script_root != project_root:
        protected_trees.add(script_root)

    unsafe_location = (
        KG_DIR in {filesystem_root, project_root, data_root, script_root}
        or any(root in KG_DIR.parents for root in protected_trees)
    )
    if unsafe_location or KG_DIR.name.lower() != "kg":
        raise ValueError(
            "The generated-output directory must be the dedicated project 'kg' "
            f"directory outside the data and script trees, got {KG_DIR}"
        )

    required_files = [
        Path(REASONING_RELEASE_BUILDER),
        Path(ORPHA_UMLS_MAPPINGS_JSON),
        Path(LOINC_PART_TEST_JSON),
        Path(LOINC_INDEX_CSV),
        *(Path(path) for path in EXTERNAL_ONTOLOGY_FILES.values()),
        Path(SKOS_SOURCE_FILE),
    ]
    missing = [str(path) for path in required_files if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Required MAKAAO build files are missing:\n  " + "\n  ".join(missing)
        )
    if REASONER_TIMEOUT_SECONDS <= 0:
        raise ValueError("MAKAAO_REASONER_TIMEOUT_SECONDS must be greater than zero")
    if REQUIRE_STRICT_REASONING_RELEASE and _robot_command() is None:
        raise RuntimeError(
            "ROBOT is required. Install `robot` or place robot.jar in "
            "project/tools, project/kg, the configured output directory, or the project root."
        )


def _replace_with_retries(source: Path, destination: Path, attempts: int = 12) -> None:
    """Move a path atomically, retrying short-lived Windows/DrvFS locks.

    Antivirus scanners, file indexers, and GUI applications can briefly hold a
    handle after a file is created or read on a Windows-mounted WSL drive. A
    bounded retry makes publication resilient without concealing persistent
    locks or permission errors.
    """
    delay = 0.10
    for attempt in range(attempts):
        try:
            os.replace(source, destination)
            return
        except PermissionError:
            if attempt + 1 >= attempts:
                raise
            time.sleep(delay)
            delay = min(delay * 1.6, 1.0)


def _iter_tree_files(root: Path) -> list[Path]:
    """Return files and symlinks below root without following directory links."""
    paths: list[Path] = []
    if not root.exists() and not root.is_symlink():
        return paths
    for current, dirnames, filenames in os.walk(root, topdown=True, followlinks=False):
        current_path = Path(current)
        retained_dirs: list[str] = []
        for dirname in dirnames:
            child = current_path / dirname
            if child.is_symlink():
                paths.append(child)
            else:
                retained_dirs.append(dirname)
        dirnames[:] = retained_dirs
        paths.extend(current_path / filename for filename in filenames)
    return sorted(paths, key=lambda path: path.relative_to(root).as_posix())


def _prune_empty_directories(root: Path) -> None:
    """Best-effort removal of stale empty subdirectories, never root itself."""
    if not root.is_dir() or root.is_symlink():
        return
    directories = [path for path in root.rglob("*") if path.is_dir() and not path.is_symlink()]
    for directory in sorted(directories, key=lambda path: len(path.parts), reverse=True):
        try:
            directory.rmdir()
        except OSError:
            pass


def commit_staged_release(
    stage_root: Path,
    final_dir: Path,
    remove_if_absent: tuple[str, ...] = (),
) -> None:
    """Publish owned outputs with rollback and Windows-safe directory handling.

    Regular top-level files are replaced atomically. Directory outputs are
    published file-by-file into the existing directory instead of renaming the
    whole non-empty directory. This avoids the common ``PermissionError`` on
    Windows-mounted WSL drives when Explorer, an editor, a terminal, or an
    indexer holds a directory handle.

    Publication remains atomic per file rather than across the complete release.
    Every destructive action is registered before it occurs, and any handled
    interruption restores the previous files from a private backup directory.
    """
    final_dir.mkdir(parents=True, exist_ok=True)
    backup_dir = Path(tempfile.mkdtemp(prefix=".makaao-backup-", dir=final_dir))
    moved_old: list[tuple[Path, Path]] = []
    installed: list[Path] = []
    created_directory_roots: list[Path] = []
    commit_succeeded = False
    rollback_succeeded = False

    def backup_file(target: Path, backup: Path) -> None:
        backup.parent.mkdir(parents=True, exist_ok=True)
        moved_old.append((backup, target))
        _replace_with_retries(target, backup)

    def backup_tree(target_root: Path, backup_root: Path) -> None:
        for target_file in _iter_tree_files(target_root):
            relative = target_file.relative_to(target_root)
            backup_file(target_file, backup_root / relative)

    def install_file(staged_file: Path, target_file: Path) -> None:
        target_file.parent.mkdir(parents=True, exist_ok=True)
        if target_file.exists() or target_file.is_symlink():
            raise RuntimeError(
                f"Internal publication error: target still exists after backup: {target_file}"
            )
        installed.append(target_file)
        _replace_with_retries(staged_file, target_file)

    def install_tree(staged_root: Path, target_root: Path) -> None:
        if target_root.exists() and not target_root.is_dir():
            raise RuntimeError(
                f"Cannot publish directory {staged_root.name}: existing target is not a directory: "
                f"{target_root}"
            )
        target_preexisted = target_root.exists()
        target_root.mkdir(parents=True, exist_ok=True)
        if not target_preexisted:
            created_directory_roots.append(target_root)
        backup_tree(target_root, backup_dir / staged_root.name)
        for staged_file in _iter_tree_files(staged_root):
            install_file(staged_file, target_root / staged_file.relative_to(staged_root))
        _prune_empty_directories(target_root)

    try:
        # Optional outputs owned by this build must not survive from an older
        # release when absent from the new stage.
        for name in remove_if_absent:
            staged_optional = stage_root / name
            target_optional = final_dir / name
            if staged_optional.exists() or staged_optional.is_symlink():
                continue
            if target_optional.is_dir() and not target_optional.is_symlink():
                backup_tree(target_optional, backup_dir / name)
                _prune_empty_directories(target_optional)
                try:
                    target_optional.rmdir()
                except OSError:
                    # An empty directory may remain when another Windows process
                    # holds a directory handle. It contains no release files.
                    pass
            elif target_optional.exists() or target_optional.is_symlink():
                backup_file(target_optional, backup_dir / name)

        for staged in sorted(stage_root.iterdir(), key=lambda item: item.name):
            target = final_dir / staged.name
            if staged.is_dir() and not staged.is_symlink():
                install_tree(staged, target)
            else:
                if target.is_dir() and not target.is_symlink():
                    raise RuntimeError(
                        f"Cannot publish file {staged.name}: existing target is a directory: {target}"
                    )
                if target.exists() or target.is_symlink():
                    backup_file(target, backup_dir / staged.name)
                install_file(staged, target)
        commit_succeeded = True
    except BaseException as original_error:
        rollback_errors: list[str] = []
        for target in reversed(installed):
            try:
                if target.is_symlink() or target.is_file():
                    target.unlink()
                elif target.exists():
                    raise RuntimeError(f"installed target unexpectedly became a directory: {target}")
            except Exception as exc:
                rollback_errors.append(f"could not remove new output {target}: {exc}")

        for backup, target in reversed(moved_old):
            try:
                if backup.exists() or backup.is_symlink():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    _replace_with_retries(backup, target)
            except Exception as exc:
                rollback_errors.append(f"could not restore {target} from {backup}: {exc}")

        for staged in stage_root.iterdir():
            target = final_dir / staged.name
            if staged.is_dir() and target.is_dir():
                _prune_empty_directories(target)
        for directory in reversed(created_directory_roots):
            try:
                directory.rmdir()
            except OSError:
                pass

        if rollback_errors:
            raise RuntimeError(
                "Release commit failed and rollback was incomplete. "
                f"Preserved backup directory: {backup_dir}\n  "
                + "\n  ".join(rollback_errors)
            ) from original_error
        rollback_succeeded = True
        raise
    finally:
        if commit_succeeded or rollback_succeeded:
            shutil.rmtree(backup_dir, ignore_errors=True)

def validate_staged_root_contents(stage_root: Path, expected_names: set[str]) -> None:
    """Prevent temporary or unexpected files from becoming permanent outputs."""
    actual_names = {path.name for path in stage_root.iterdir()}
    missing = sorted(expected_names - actual_names)
    unexpected = sorted(actual_names - expected_names)
    if missing or unexpected:
        details = []
        if missing:
            details.append("missing staged outputs: " + ", ".join(missing))
        if unexpected:
            details.append("unexpected staged outputs: " + ", ".join(unexpected))
        raise RuntimeError("Invalid staged release contents: " + "; ".join(details))

def _strip_controlled_java_options(value: str, variable_name: str) -> str:
    """Remove heap and Java-temp settings that this build controls explicitly."""
    if not value.strip():
        return ""
    try:
        tokens = shlex.split(value, posix=True)
    except ValueError as exc:
        raise ValueError(
            f"{variable_name} contains invalid shell-style quoting: {exc}"
        ) from exc
    retained = [
        token
        for token in tokens
        if not token.startswith("-Xmx")
        and not token.startswith("-Djava.io.tmpdir=")
    ]
    return shlex.join(retained)


def build_reasoner_environment(stage_java_tmp: Path) -> dict[str, str]:
    """Return a controlled Java environment for ROBOT/HermiT subprocesses.

    Java installations can process ``JDK_JAVA_OPTIONS``, ``JAVA_TOOL_OPTIONS``,
    command-line VM arguments, and ``_JAVA_OPTIONS`` in implementation-specific
    precedence orders. Remove inherited heap and ``java.io.tmpdir`` settings
    from all three environment variables, preserve unrelated options, and add
    the build-controlled values once through ``JAVA_TOOL_OPTIONS``. A local JAR
    also receives the same explicit ``-Xmx`` command-line argument. The caller's
    process environment is not modified.
    """
    stage_java_tmp = stage_java_tmp.resolve()
    stage_java_tmp.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()

    for variable_name in (
        "JDK_JAVA_OPTIONS",
        "JAVA_TOOL_OPTIONS",
        "_JAVA_OPTIONS",
    ):
        cleaned = _strip_controlled_java_options(
            environment.get(variable_name, ""), variable_name
        )
        if cleaned:
            environment[variable_name] = cleaned
        else:
            environment.pop(variable_name, None)

    configured_java_options = shlex.join(
        [f"-Xmx{JAVA_MAX_HEAP}", f"-Djava.io.tmpdir={stage_java_tmp}"]
    )
    existing_tool_options = environment.get("JAVA_TOOL_OPTIONS", "").strip()
    environment["JAVA_TOOL_OPTIONS"] = (
        f"{existing_tool_options} {configured_java_options}"
    ).strip()
    environment["TMPDIR"] = str(stage_java_tmp)
    environment["TMP"] = str(stage_java_tmp)
    environment["TEMP"] = str(stage_java_tmp)
    return environment


# ===================== MAIN =====================
def main():
    validate_build_prerequisites()
    KG_DIR.mkdir(parents=True, exist_ok=True)
    stage_root = Path(tempfile.mkdtemp(prefix=".makaao-build-", dir=KG_DIR))
    stage_kg = stage_root / Path(OUTPUT_OWL_ENRICHED).name
    stage_tbox = stage_root / Path(OUTPUT_OWL_TBOX).name
    stage_reasoning = stage_root / REASONING_OUTPUT_DIRNAME
    stage_java_tmp = stage_root / ".java-tmp"
    stage_java_tmp.mkdir(parents=True, exist_ok=True)

    reasoner_env = build_reasoner_environment(stage_java_tmp)

    try:
        # create empty KG
        g = init_graph()
    
        # Optional enrichment datasets
    
        # orpha_hpo_links is a dict where keys are Orphanet codes and values are [hpo, freq]; from the csv where we can find the HPO terms linked to the Orphanet diseases (from en_product4.xml an Orphanet file listing HPO terms associated to Orphanet dieases with their frequency)
        orpha_hpo_links = (
            grouped(read_csv_rows(ORPHANET_HPO_LINKS), "orpha_code")
            if os.path.exists(ORPHANET_HPO_LINKS)
            else {}
        )
        #print(orpha_hpo_links)
    
        # Names read from code_names.csv, including LOINC Parts and Terms.
        if os.path.exists(CODE_NAMES_CSV):
            up_names, up_urls = read_code_names_uniprot(
                CODE_NAMES_CSV
            )  # 2 dictionnaries: {uniprot_id: english_name} and {uniprot_id: uniprot_url}
    
            umls_cn_names, umls_cn_urls = read_code_names_umls(
                CODE_NAMES_CSV
            )  # dictinoary: {CUI: english_name} (CUIs of differnt concept sof interest: targets, diseases...)
            umls_synonyms = read_code_names_umls_synonyms(CODE_NAMES_CSV)
            hpo_cn_names, _ = read_code_names_hpo(
                CODE_NAMES_CSV
            )  # dict: {HP:nnnnnnn: english_name}
            orpha_cn_names, orpha_cn_urls = read_code_names_orpha(
                CODE_NAMES_CSV
            )
            chebi_cn_names, chebi_cn_urls = read_code_names_chebi(
                CODE_NAMES_CSV
            )  # 2 dictionaries: {chebi_id: english_name} and {chebi_id: chebi_url}
            loinc_part_names, loinc_term_names = read_code_names_loinc(
                CODE_NAMES_CSV
            )
            print(
                "Loaded names — "
                f"UniProt:{len(up_names)} "
                f"UMLS:{len(umls_cn_names)} "
                f"HPO:{len(hpo_cn_names)} "
                f"ORPHA:{len(orpha_cn_names)} "
                f"ChEBI:{len(chebi_cn_names)} "
                f"LOINC Parts:{len(loinc_part_names)} "
                f"LOINC Terms:{len(loinc_term_names)}"
            )
        else:
            up_names, up_urls = {}, {}
            umls_cn_names, umls_cn_urls = {}, {}
            umls_synonyms = {}
            hpo_cn_names = {}
            orpha_cn_names, orpha_cn_urls = {}, {}
            chebi_cn_names, chebi_cn_urls = {}, {}
            loinc_part_names, loinc_term_names = {}, {}
            print(
                f"WARN: {CODE_NAMES_CSV} not found; no enrichment labels were loaded."
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

        orpha_umls_mapping_pairs = read_orpha_umls_mappings(
            ORPHA_UMLS_MAPPINGS_JSON
        )
        print(
            "Loaded ORPHA/UMLS enrichment mappings: "
            f"pairs={len(orpha_umls_mapping_pairs)}"
        )
    
        # this function call add AAb with their class hierarchy, and their positivity phenotypes. It also add their targets (UMLS, CHebi, or Uniprot), instantiate all these things and link them together as needed
        positivity_instances_by_hpo = build_core(
            g,
            data,
            up_names,
            up_urls,
            umls_cn_names,
            umls_cn_urls,
            hpo_cn_names,
            chebi_cn_names,
            chebi_cn_urls,
            umls_synonyms=umls_synonyms,
        )
    
        # this function call add diseases linked to AAbs, instantiate them and link them to AAb instances. Diseases can be Orphanet diseases (preferred) or UMLS CUIs (with possible SNOMED mappings)
        process_diseases(
            g,
            data,
            orpha_hpo_links,
            positivity_instances_by_hpo,
            umls_cn_names,
            umls_cn_urls,
            hpo_cn_names,
            umls_synonyms=umls_synonyms,
        )
        keep = data["indices"]
        process_loinc_mappings(
            g,
            LOINC_INDEX_CSV,
            keep,
            loinc_part_names,
            loinc_term_names,
        )

        orpha_umls_triples_added = add_orpha_umls_close_matches(
            g,
            orpha_umls_mapping_pairs,
            umls_names=umls_cn_names,
            umls_urls=umls_cn_urls,
            orpha_names=orpha_cn_names,
            umls_synonyms=umls_synonyms,
        )
        print(
            "ORPHA/UMLS closeMatch enrichment: "
            f"pairs={len(orpha_umls_mapping_pairs)} "
            f"closeMatch_triples_added={orpha_umls_triples_added}"
        )

        append_fair_metadata(g)

        # Run the mutating collision audit at the latest possible point: after
        # every input-derived and metadata enrichment step, immediately before
        # final validation, TBox extraction, serialization, and reasoning.
        label_collision_rows = add_label_collision_close_matches(g)
        linked_label_pairs = sum(
            1 for row in label_collision_rows if row["decision"] == "linked"
        )
        label_close_match_triples_added = sum(
            int(row["close_match_triples_added"]) for row in label_collision_rows
        )
        print(
            "Class-label collision audit: "
            f"candidate_pairs={len(label_collision_rows)} "
            f"linked_pairs={linked_label_pairs} "
            f"closeMatch_triples_added={label_close_match_triples_added}"
        )

        # Set file-level metadata only after the last graph-mutating audit so
        # void:triples records the final provenance-bearing KG size.
        set_output_file_metadata(g, Path(OUTPUT_OWL_ENRICHED).name)

        labelled_cui_resources = validate_local_cui_labels(g)
        print(f"Local UMLS resources with labels: {labelled_cui_resources}")
        validate_graph_iris(g, "assembled MAKAAO KG")
    
        # Export ontology (T-box) separately (no instances)
        tbox = extract_tbox(g, str(MAKAAO))
        validate_graph_iris(tbox, "MAKAAO TBox")
        _serialize_rdfxml_verified(tbox, stage_tbox, "standalone MAKAAO TBox")
        print(f"Staged {stage_tbox}  triples={len(tbox)}")
    
        _serialize_rdfxml_verified(
            g, stage_kg, "canonical provenance-bearing MAKAAO KG"
        )
        print(f"Staged {stage_kg}  triples={len(g)}")

        # Build the reasoner-oriented projection in memory. It is consumed by
        # the strict reasoning helper but is not published as a duplicate root file.
        no_reification = build_non_reified_graph(g)
        print(
            f"Built in-memory non-reified projection  triples={len(no_reification)} "
            "(Statement/Document resources removed)"
        )
    
        # Build the only external module set directly from the pinned upstream
        # ontology files, then validate/classify the strict OWL 2 DL release.
        reasoning_result = build_strict_reasoning_release(
            g,
            tbox,
            no_reification,
            stage_reasoning,
            stage_kg,
            stage_tbox,
            reasoner_env,
        )
        label_collision_report_path = (
            stage_reasoning / "reports" / LABEL_COLLISION_REPORT_FILENAME
        )
        write_label_collision_report(
            label_collision_rows,
            label_collision_report_path,
        )
        print(f"Class-label close-match report: {label_collision_report_path}")
        root_catalog = write_root_reasoning_catalog(stage_root, reasoning_result)
        finalize_release_manifest_and_checksums(
            stage_root,
            stage_reasoning,
            stage_kg,
            stage_tbox,
            root_catalog,
        )
        # Java temporary files are build-time artifacts and must not be moved
        # into the permanent release by the directory-level commit.
        shutil.rmtree(stage_java_tmp, ignore_errors=True)
        validate_staged_root_contents(
            stage_root,
            {
                stage_kg.name,
                stage_tbox.name,
                "catalog-imports.xml",
                REASONING_OUTPUT_DIRNAME,
            },
        )
        commit_staged_release(
            stage_root,
            KG_DIR,
            remove_if_absent=(
                ".java-tmp",
                "imports",
                f"makaao_kg_{version}_no-reification.rdf",
                f"makaao_kg_{version}_import-closure.rdf",
                f"makaao_kg_{version}_no-reification_import-closure.rdf",
                f"makaao_kg_{version}_no-reification_import-closure.owl.xml",
                f"makaao_kg_{version}_no-reification_import-closure_reasoned.owl",
                REASONING_OUTPUT_DIRNAME,
            ),
        )
        print(
            f"MAKAAO build complete: script={SCRIPT_VERSION} kg={version} "
            f"triples={len(g)} reasoning={reasoning_result['status']}"
        )
    finally:
        if stage_root.exists():
            shutil.rmtree(stage_root, ignore_errors=True)
if __name__ == "__main__":  
    main()