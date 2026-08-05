import sys
import importlib.util
import json
from pathlib import Path

import pytest
from rdflib import Graph, Literal, OWL, RDF, RDFS, URIRef
from rdflib.namespace import SKOS


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return load_module(Path("scripts/03_build_kg_from_tables.py"), "build_kg")


def fresh_graph(mod):
    return mod.init_graph()


def add_class(graph, iri, label, parent=None):
    graph.add((iri, RDF.type, OWL.Class))
    if parent is not None:
        graph.add((iri, RDFS.subClassOf, parent))
    graph.add((iri, RDFS.label, Literal(label)))


def test_current_script_version_and_reasoning_helper_are_synchronized(mod):
    helper = load_module(
        Path("scripts/03-1_build_reasoning_release.py"), "reasoning_helper_version"
    )
    assert mod.SCRIPT_VERSION == "1.2.28"
    assert helper.SCRIPT_VERSION == mod.SCRIPT_VERSION
    assert helper.SCRIPT_ITERATION == mod.SCRIPT_ITERATION


def test_identifier_normalizers(mod):
    assert mod.make_valid(" A:B/C|D ") == "_A_B_C_D_"
    assert mod.canonical_loinc_code("https://loinc.org/LP123-4", "part") == "LP123-4"
    assert mod.canonical_loinc_code("LOINC:1234-5", "term") == "1234-5"
    assert mod.canonical_hpo_code("HP_0000001") == "HP:0000001"
    assert mod.canonical_orpha_code("ORPHA:123") == "123"
    assert mod.canonical_cui_code("CUI:C1234567") == "C1234567"


def test_label_deduplication_and_cui_fallbacks(mod):
    graph = fresh_graph(mod)
    node = mod.MAKAAO["Example"]
    assert mod.add_label(graph, node, RDFS.label, "Example") is True
    assert mod.add_label(graph, node, RDFS.label, "Example") is False
    assert mod.add_pref(graph, node, "Example") is True
    assert mod.add_pref(graph, node, "Example") is False

    cui_class, cui_instance = mod.ensure_cui_class_and_instance(graph, "C1234567")
    for resource in (cui_class, cui_instance):
        assert Literal("C1234567") in set(graph.objects(resource, RDFS.label))
        assert Literal("C1234567") in set(graph.objects(resource, SKOS.prefLabel))
    assert mod.validate_local_cui_labels(graph) == 2


def test_exact_label_audit_searches_only_the_three_permitted_kind_pairs(mod, tmp_path):
    graph = fresh_graph(mod)

    cui_protein = mod.MAKAAO["CUI_C0000001"]
    uniprot = mod.MAKAAO["UP_P12345"]
    add_class(graph, cui_protein, "Protein X", mod.MAKAAO.CUI)
    add_class(graph, uniprot, "  protein   x  ")

    antibody = mod.MAKAAO["aab_1"]
    loinc_part = mod.MAKAAO_LOINC["LP12345-6"]
    add_class(graph, antibody, "Marker Y", mod.MAKAAO.Autoantibody)
    add_class(graph, loinc_part, "marker y", mod.MAKAAO.LoincPart)

    ordo = URIRef("http://www.orpha.net/ORDO/Orphanet_123")
    cui_disease = mod.MAKAAO["CUI_C0000002"]
    add_class(graph, ordo, "Disease Z", mod.MAKAAO.AutoimmuneDisease)
    add_class(graph, cui_disease, "disease z", mod.MAKAAO.CUI)

    # Exact collisions in prohibited combinations must never become candidates.
    hpo = URIRef("http://purl.obolibrary.org/obo/HP_0000001")
    chebi = URIRef("http://purl.obolibrary.org/obo/CHEBI_23367")
    second_uniprot = mod.MAKAAO["UP_Q99999"]
    third_uniprot = mod.MAKAAO["UP_Q99998"]
    add_class(graph, hpo, "Forbidden collision")
    add_class(graph, chebi, "Forbidden collision")
    add_class(graph, second_uniprot, "Same UniProt label")
    add_class(graph, third_uniprot, "Same UniProt label")

    rows = mod.add_label_collision_close_matches(graph)
    assert len(rows) == 3
    assert {row["decision"] for row in rows} == {"linked"}
    assert sum(int(row["close_match_triples_added"]) for row in rows) == 6
    expected_kind_pairs = frozenset(
        {
            frozenset({"umls_cui", "uniprot"}),
            frozenset({"autoantibody", "loinc_part"}),
            frozenset({"ordo_disease", "umls_cui"}),
        }
    )
    assert mod.LABEL_MATCH_KIND_PAIRS == expected_kind_pairs
    assert {
        frozenset((row["class_1_kind"], row["class_2_kind"])) for row in rows
    } == expected_kind_pairs

    for left, right in (
        (cui_protein, uniprot),
        (antibody, loinc_part),
        (ordo, cui_disease),
    ):
        assert (left, SKOS.closeMatch, right) in graph
        assert (right, SKOS.closeMatch, left) in graph

    assert (hpo, SKOS.closeMatch, chebi) not in graph
    assert (second_uniprot, SKOS.closeMatch, third_uniprot) not in graph

    report = tmp_path / "class-label-close-match-report.tsv"
    mod.write_label_collision_report(rows, report)
    text = report.read_text(encoding="utf-8")
    assert "review_required" not in text
    assert text.count("\n") == 4



def test_exact_label_audit_excludes_all_other_kind_combinations(mod):
    graph = fresh_graph(mod)
    shared_label = "One shared label"

    resources = {
        "autoantibody": mod.MAKAAO["aab_900"],
        "loinc_part": mod.MAKAAO_LOINC["LP900-1"],
        "umls_cui": mod.MAKAAO["CUI_C0000900"],
        "uniprot": mod.MAKAAO["UP_P00900"],
        "ordo_disease": URIRef("http://www.orpha.net/ORDO/Orphanet_900"),
        "hpo": URIRef("http://purl.obolibrary.org/obo/HP_0000900"),
        "chebi": URIRef("http://purl.obolibrary.org/obo/CHEBI_900"),
    }

    add_class(graph, resources["autoantibody"], shared_label, mod.MAKAAO.Autoantibody)
    add_class(graph, resources["loinc_part"], shared_label, mod.MAKAAO.LoincPart)
    add_class(graph, resources["umls_cui"], shared_label, mod.MAKAAO.CUI)
    add_class(graph, resources["uniprot"], shared_label)
    add_class(
        graph,
        resources["ordo_disease"],
        shared_label,
        mod.MAKAAO.AutoimmuneDisease,
    )
    add_class(graph, resources["hpo"], shared_label)
    add_class(graph, resources["chebi"], shared_label)

    rows = mod.add_label_collision_close_matches(graph)
    actual_pairs = {
        frozenset((row["class_1_kind"], row["class_2_kind"])) for row in rows
    }
    expected_pairs = {
        frozenset({"umls_cui", "uniprot"}),
        frozenset({"autoantibody", "loinc_part"}),
        frozenset({"ordo_disease", "umls_cui"}),
    }
    assert actual_pairs == expected_pairs
    assert len(rows) == 3

    linked_resources = {
        frozenset((URIRef(row["class_1"]), URIRef(row["class_2"]))) for row in rows
    }
    expected_resources = {
        frozenset((resources["umls_cui"], resources["uniprot"])),
        frozenset((resources["autoantibody"], resources["loinc_part"])),
        frozenset((resources["ordo_disease"], resources["umls_cui"])),
    }
    assert linked_resources == expected_resources

    for left_kind, left in resources.items():
        for right_kind, right in resources.items():
            if str(left) >= str(right):
                continue
            pair = frozenset((left, right))
            if pair in expected_resources:
                assert (left, SKOS.closeMatch, right) in graph
                assert (right, SKOS.closeMatch, left) in graph
            else:
                assert (left, SKOS.closeMatch, right) not in graph, (
                    f"Unexpected mapping for {left_kind} and {right_kind}"
                )
                assert (right, SKOS.closeMatch, left) not in graph, (
                    f"Unexpected reverse mapping for {left_kind} and {right_kind}"
                )


def test_exact_label_audit_uses_pref_and_alt_labels(mod):
    graph = fresh_graph(mod)
    cui = mod.MAKAAO["CUI_C0000100"]
    uniprot = mod.MAKAAO["UP_P00001"]
    graph.add((cui, RDF.type, OWL.Class))
    graph.add((cui, RDFS.subClassOf, mod.MAKAAO.CUI))
    graph.add((cui, SKOS.altLabel, Literal("  Alias   Protein ")))
    graph.add((uniprot, RDF.type, OWL.Class))
    graph.add((uniprot, SKOS.prefLabel, Literal("alias protein")))

    rows = mod.add_label_collision_close_matches(graph)
    assert len(rows) == 1
    assert rows[0]["shared_normalized_labels"] == "alias protein"
    assert rows[0]["close_match_triples_added"] == 2
    assert (cui, SKOS.closeMatch, uniprot) in graph
    assert (uniprot, SKOS.closeMatch, cui) in graph


def test_exact_label_audit_is_idempotent_for_existing_mappings(mod):
    graph = fresh_graph(mod)
    ordo = URIRef("http://www.orpha.net/ORDO/Orphanet_999")
    cui = mod.MAKAAO["CUI_C0000999"]
    add_class(graph, ordo, "Example disease", mod.MAKAAO.AutoimmuneDisease)
    add_class(graph, cui, "example disease", mod.MAKAAO.CUI)
    graph.add((ordo, SKOS.closeMatch, cui))
    graph.add((cui, SKOS.closeMatch, ordo))

    rows = mod.add_label_collision_close_matches(graph)
    assert len(rows) == 1
    assert rows[0]["decision"] == "linked"
    assert rows[0]["close_match_triples_added"] == 0

def test_orpha_umls_dictionary_must_be_symmetric(mod, tmp_path):
    path = tmp_path / "orpha_umls_mappings.json"
    payload = {
        "orpha_to_umls": {"123": ["C1234567"]},
        "umls_to_orpha": {"C1234567": ["123"]},
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    assert mod.read_orpha_umls_mappings(path) == [("123", "C1234567")]

    payload["umls_to_orpha"] = {}
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="not symmetric"):
        mod.read_orpha_umls_mappings(path)


def test_empty_loinc_primary_array_creates_part_without_fabricated_term(mod, tmp_path):
    graph = fresh_graph(mod)
    index_csv = tmp_path / "index_loinc.csv"
    index_csv.write_text("aab_id,loinc_id\n1,LP100-1\n", encoding="utf-8")
    mappings = tmp_path / "loinc_part_test_dict.json"
    mappings.write_text(json.dumps({"LP100-1": []}), encoding="utf-8")
    labels = tmp_path / "loinc_labels.json"
    labels.write_text(
        json.dumps({"parts": {"LP100-1": "Example Part"}, "tests": {}}),
        encoding="utf-8",
    )

    mod.process_loinc_mappings(
        graph,
        str(index_csv),
        {"1"},
        part_test_json=str(mappings),
        labels_json=str(labels),
    )

    part_class = mod.MAKAAO_LOINC["LP100-1"]
    part_instance = mod.MAKAAO["loinc_LP100-1_instance"]
    assert (part_class, RDF.type, OWL.Class) in graph
    assert (part_instance, RDF.type, part_class) in graph
    assert not list(graph.triples((None, mod.MAKAAO.hasLoincComponent, None)))


def test_non_reified_projection_removes_records_but_preserves_assertion(mod):
    graph = fresh_graph(mod)
    subject = mod.MAKAAO["subject"]
    predicate = mod.MAKAAO["predicate"]
    obj = mod.MAKAAO["object"]
    relation = mod.MAKAAO["r1"]
    document = mod.MAKAAO["document_1"]

    graph.add((subject, predicate, obj))
    graph.add((relation, RDF.type, RDF.Statement))
    graph.add((relation, RDF.type, mod.MAKAAO.Relation))
    graph.add((relation, RDF.subject, subject))
    graph.add((relation, RDF.predicate, predicate))
    graph.add((relation, RDF.object, obj))
    graph.add((relation, mod.PROV.wasDerivedFrom, document))
    graph.add((document, RDF.type, mod.MAKAAO.Document))

    projected = mod.build_non_reified_graph(graph)
    assert (subject, predicate, obj) in projected
    assert not list(projected.triples((relation, None, None)))
    assert not list(projected.triples((document, None, None)))
    assert not list(projected.triples((None, RDF.subject, None)))


def test_invalid_iri_is_rejected(mod):
    graph = Graph()
    graph.add((URIRef("http://example.org/has space"), RDF.type, OWL.Class))
    with pytest.raises(ValueError, match="Invalid IRI"):
        mod.validate_graph_iris(graph)
