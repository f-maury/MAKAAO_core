from __future__ import annotations

from types import ModuleType
from pathlib import Path
import re

import pytest
from rdflib import Graph, Namespace, URIRef, RDF, RDFS, OWL, Literal
from rdflib.namespace import SKOS


# ----------------------------
# Script runner with IN_PATH / OUT_PATH overrides
# ----------------------------
def _find_script_path() -> Path:
    candidates = [
        Path("scripts/04_make_lite_graph_from_makaao-kg.py"),
        Path("04_make_lite_graph_from_makaao-kg.py"),
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        "Cannot find 04_make_lite_graph_from_makaao-kg.py in: "
        + ", ".join(str(c) for c in candidates)
    )


def exec_script_with_overrides(script_path: Path, *, in_path: Path, out_path: Path) -> ModuleType:
    """
    Execute a script that runs on import, but first rewrite IN_PATH and OUT_PATH
    in its source so it uses our test files.
    """
    code = script_path.read_text(encoding="utf-8", errors="ignore")

    code = re.sub(
        r'^\s*IN_PATH\s*=\s*.*$',
        f'IN_PATH = r"{str(in_path)}"',
        code,
        flags=re.MULTILINE,
    )
    code = re.sub(
        r'^\s*OUT_PATH\s*=\s*.*$',
        f'OUT_PATH = r"{str(out_path)}"',
        code,
        flags=re.MULTILINE,
    )

    m = ModuleType("kg_lite_exec")
    exec(compile(code, str(script_path), "exec"), m.__dict__)
    return m


# ----------------------------
# Invariants for the lite KG
# ----------------------------
MAK = Namespace("http://makaao.inria.fr/kg/")
OBO = Namespace("http://purl.obolibrary.org/obo/")
UNI = Namespace("http://purl.uniprot.org/core/")
SIO = Namespace("http://semanticscience.org/resource/")
BIOL = Namespace("https://w3id.org/biolink/vocab/")
BAO = Namespace("http://www.bioassayontology.org/bao#")
PROV = Namespace("http://www.w3.org/ns/prov#")
ORDO = Namespace("http://www.orpha.net/ORDO/")

KEEP_PRED = {RDFS.label, SKOS.prefLabel, RDFS.subClassOf}
REL_PRED = {
    SIO["SIO_001403"],         # is_associated_with
    BIOL["biomarker_for"],
    BIOL["has_biomarker"],
    BAO["BAO_0000598"],        # is_target_for
    BAO["BAO_0000211"],        # has_target
    SIO["SIO_001279"],         # has_phenotype
    SIO["SIO_001280"],         # is_phenotype_of
}
ALLOWED_PRED = KEEP_PRED | REL_PRED | {RDF.type}

EXCLUDE_CLASSES = {URIRef(MAK["Document"]), URIRef(MAK["Relation"])}

CHEBI_23367 = URIRef(OBO["CHEBI_23367"])
PROTEIN = URIRef(UNI["Protein"])
TARGET = URIRef(MAK["Target"])
AUTOIMMUNE_DZ = URIRef(MAK["AutoimmuneDisease"])
ORPHANET_C001 = URIRef(ORDO["Orphanet_C001"])


def _tail(u: URIRef) -> str:
    return str(u).rsplit("/", 1)[-1]


def _is_instance_uri(u: URIRef) -> bool:
    s = str(u)
    if s.endswith("_instance"):
        return True
    if s.startswith(str(MAK)):
        t = _tail(u)
        if t.startswith("document_"):
            return True
        if t.startswith("r") and t[1:].isdigit():
            return True
    return False


def _involves_prov(t) -> bool:
    s, p, o = t
    return (
        (isinstance(s, URIRef) and str(s).startswith(str(PROV)))
        or (isinstance(p, URIRef) and str(p).startswith(str(PROV)))
        or (isinstance(o, URIRef) and str(o).startswith(str(PROV)))
    )


def assert_lite_graph_invariants(g: Graph) -> None:
    # Non-empty
    assert len(g) > 0

    # No PROV anywhere
    prov_triples = [t for t in g if _involves_prov(t)]
    assert not prov_triples, f"Found PROV-involving triples, e.g. {prov_triples[:3]!r}"

    # Only whitelisted predicates (+ rdf:type)
    bad_preds = sorted({p for _, p, _ in g if p not in ALLOWED_PRED}, key=str)
    assert not bad_preds, f"Found non-whitelisted predicates, e.g. {bad_preds[:10]!r}"

    # No instances remain
    bad_nodes = sorted(
        {n for (s, _, o) in g for n in (s, o) if isinstance(n, URIRef) and _is_instance_uri(n)},
        key=str,
    )
    assert not bad_nodes, f"Found instance URIs in output, e.g. {bad_nodes[:10]!r}"

    # KEY REQUIREMENT: the only rdf:type assertions are owl:Class
    bad_types = sorted(
        {(s, o) for s, o in g.subject_objects(RDF.type) if o != OWL.Class},
        key=lambda x: (str(x[0]), str(x[1])),
    )
    assert not bad_types, (
        "Output contains rdf:type assertions other than owl:Class. "
        f"Examples: {bad_types[:10]!r}"
    )

    # Excluded classes never appear as subject or object
    banned_hits = []
    for banned in EXCLUDE_CLASSES:
        if (banned, None, None) in g:
            banned_hits.append(("as_subject", banned))
        if (None, None, banned) in g:
            banned_hits.append(("as_object", banned))
    assert not banned_hits, f"Excluded classes appear in output: {banned_hits!r}"

    # Literals only on label predicates
    bad_literal_triples = [
        (s, p, o)
        for (s, p, o) in g
        if isinstance(o, Literal) and p not in {RDFS.label, SKOS.prefLabel}
    ]
    assert not bad_literal_triples, (
        "Found literals on non-label predicates. "
        f"Examples: {bad_literal_triples[:10]!r}"
    )

    # Sanity checks on enforced roots (only if present)
    if (PROTEIN, RDFS.subClassOf, None) in g:
        supers = {o for o in g.objects(PROTEIN, RDFS.subClassOf)}
        assert supers == {TARGET}, f"Protein superclasses must be exactly {{Target}}, got {supers!r}"

    if (CHEBI_23367, RDFS.subClassOf, None) in g:
        supers = {o for o in g.objects(CHEBI_23367, RDFS.subClassOf)}
        assert supers == {TARGET}, f"CHEBI_23367 superclasses must be exactly {{Target}}, got {supers!r}"

    if (ORPHANET_C001, RDFS.subClassOf, None) in g:
        supers = {o for o in g.objects(ORPHANET_C001, RDFS.subClassOf)}
        assert AUTOIMMUNE_DZ in supers, f"Orphanet_C001 must be subclassOf AutoimmuneDisease, got {supers!r}"


# ----------------------------
# Fixtures (your repo files)
# ----------------------------
@pytest.fixture(scope="session")
def full_input_path() -> Path:
    p = Path("kg/makaao_kg_sample.rdf")
    assert p.exists(), "Missing input fixture: kg/makaao_kg_sample.rdf"
    return p


@pytest.fixture(scope="session")
def lite_golden_path() -> Path:
    p = Path("kg/makaao_kg_lite_sample.rdf")
    assert p.exists(), "Missing golden fixture: kg/makaao_kg_lite_sample.rdf"
    return p


# ----------------------------
# Tests
# ----------------------------
def test_04_script_produces_lite_graph_from_full(full_input_path: Path, tmp_path: Path):
    out_path = tmp_path / "lite_out.rdf"
    script_path = _find_script_path()

    exec_script_with_overrides(script_path, in_path=full_input_path, out_path=out_path)

    assert out_path.exists(), "Script did not write OUT_PATH"

    g = Graph()
    g.parse(str(out_path), format="xml")
    assert_lite_graph_invariants(g)


def test_04_expected_lite_graph_file_is_valid(lite_golden_path: Path):
    g = Graph()
    g.parse(str(lite_golden_path), format="xml")
    assert_lite_graph_invariants(g)


def test_04_script_is_idempotent_on_lite_graph(lite_golden_path: Path, tmp_path: Path):
    out_path = tmp_path / "lite_out_idempotent.rdf"
    script_path = _find_script_path()

    exec_script_with_overrides(script_path, in_path=lite_golden_path, out_path=out_path)

    g = Graph()
    g.parse(str(out_path), format="xml")
    assert_lite_graph_invariants(g)
