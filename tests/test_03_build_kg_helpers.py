import importlib.util
from pathlib import Path
from rdflib import Graph
from rdflib.namespace import RDFS

def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod

def test_make_valid_and_hp_uri_and_id_canon():
    mod = load_module(Path("scripts/03_build_kg_from_tables.py"), "build_kg")

    assert mod.make_valid("a b:c/d|e") == "a_b_c_d_e"
    assert str(mod.hp_to_obo_uri("HP_0003731")) == "http://purl.obolibrary.org/obo/HP_0003731"
    assert str(mod.hp_to_obo_uri("HP:0003731")) == "http://purl.obolibrary.org/obo/HP_0003731"

    assert mod.canon_uniprot_id("P12345") == "P12345"
    assert mod.canon_uniprot_id("uniprot:P12345") == "UP:P12345"
    assert mod.canon_chebi_id("CHEBI:23367") == "CHEBI:23367"
    assert mod.canon_chebi_id("obo:CHEBI_23367") == "CHEBI:23367"

def test_add_label_deduplication_rules():
    mod = load_module(Path("scripts/03_build_kg_from_tables.py"), "build_kg")
    g = Graph()
    n = mod.MAKAAO["X"]

    # first label added
    assert mod.add_label(g, n, RDFS.label, "Quadriceps muscle weakness") is True
    # duplicate same prop rejected
    assert mod.add_label(g, n, RDFS.label, "Quadriceps muscle weakness") is False
    # cross-duplicate between rdfs:label and skos:prefLabel allowed once
    assert mod.add_pref(g, n, "Quadriceps muscle weakness") is True
    # but further duplicates rejected
    assert mod.add_pref(g, n, "Quadriceps muscle weakness") is False

    # duplicates across unrelated props are rejected
    assert mod.add_pref(g, n, "X") is True
    assert mod.add_label(g, n, RDFS.label, "X") is False
