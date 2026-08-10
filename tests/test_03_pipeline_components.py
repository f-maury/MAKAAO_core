from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest
from rdflib import BNode, Graph, Literal, OWL, RDF, RDFS, URIRef
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
    return load_module(Path("scripts/03_build_kg_from_tables.py"), "build_kg_components")


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_processed_tables(base: Path) -> None:
    write_csv(
        base / "index_name_en.csv",
        ["index", "name_en"],
        [
            {"index": "18", "name_en": "Autoantibody"},
            {"index": "1", "name_en": "Protein X antibody"},
            {"index": "2", "name_en": "Second antibody"},
        ],
    )
    write_csv(
        base / "index_parent_index.csv",
        ["index", "parent_index"],
        [
            {"index": "1", "parent_index": "18"},
            {"index": "2", "parent_index": "1"},
            {"index": "2", "parent_index": "2"},
        ],
    )
    write_csv(
        base / "index_syn_source_en.csv",
        ["index", "syns_en", "syns_en_source"],
        [
            {"index": "1", "syns_en": "Protein X autoantibody", "syns_en_source": "PMID:123; https://example.org/paper"},
            {"index": "99", "syns_en": "ignored", "syns_en_source": "PMID:999"},
        ],
    )
    write_csv(
        base / "index_syn_source_fr.csv",
        ["index", "syns_fr", "syns_fr_source"],
        [{"index": "1", "syns_fr": "autoanticorps protéine X", "syns_fr_source": "PMID:124"}],
    )
    write_csv(
        base / "index_hpo_id.csv",
        ["index", "hpo_id"],
        [
            {"index": "1", "hpo_id": "HP:0000001 | HP_0000002"},
            {"index": "2", "hpo_id": ""},
        ],
    )
    write_csv(
        base / "index_cui_source.csv",
        ["index", "umls_target_cui", "umls_pmids"],
        [{"index": "1", "umls_target_cui": "CUI:C0000001", "umls_pmids": "PMID:200"}],
    )
    write_csv(
        base / "index_uniprot_source.csv",
        ["index", "uniprot_target_id", "uniprot_pmids"],
        [
            {"index": "1", "uniprot_target_id": "https://rest.uniprot.org/uniprotkb/P12345.fasta", "uniprot_pmids": "PMID:201"},
            {"index": "2", "uniprot_target_id": "Q99999", "uniprot_pmids": ""},
        ],
    )
    write_csv(
        base / "index_chebi_source.csv",
        ["index", "chebi_target_id", "chebi_pmids"],
        [{"index": "1", "chebi_target_id": "CHEBI:23367", "chebi_pmids": "PMID:202"}],
    )
    write_csv(
        base / "index_disease_source.csv",
        ["index", "related_diseases_id", "diseases_pmids"],
        [
            {"index": "1", "related_diseases_id": "ORPHA:123|C0000009|free disease", "diseases_pmids": "PMID:300"},
            {"index": "99", "related_diseases_id": "ORPHA:999", "diseases_pmids": "PMID:999"},
        ],
    )


def test_synthetic_processed_table_components_metadata_and_tbox(
    mod, tmp_path, monkeypatch
):
    processed = tmp_path / "processed"
    processed.mkdir()
    make_processed_tables(processed)

    code_names = tmp_path / "code_names.csv"
    write_csv(
        code_names,
        ["source", "id", "name", "url"],
        [
            {"source": "ORPHA", "id": "123", "name": "Disease Z", "url": "https://www.orpha.net/en/disease/detail/123"},
            {"source": "UMLS", "id": "C0000001", "name": "Protein X", "url": "https://uts.nlm.nih.gov/uts/umls/concept/C0000001"},
            {"source": "UMLS", "id": "C0000009", "name": "Disease CUI", "url": "https://uts.nlm.nih.gov/uts/umls/concept/C0000009"},
            {"source": "HPO", "id": "HP:0000001", "name": "All", "url": "http://purl.obolibrary.org/obo/HP_0000001"},
            {"source": "HPO", "id": "HP:0000002", "name": "Abnormality", "url": "http://purl.obolibrary.org/obo/HP_0000002"},
            {"source": "UniProt", "id": "P12345", "name": "Protein X", "url": "https://www.uniprot.org/uniprotkb/P12345"},
            {"source": "ChEBI", "id": "CHEBI:23367", "name": "molecular entity", "url": "http://purl.obolibrary.org/obo/CHEBI_23367"},
        ],
    )
    monkeypatch.setattr(mod, "CODE_NAMES_CSV", str(code_names))
    # This is deliberately a synthetic component fixture. The separate smoke
    # test exercises the committed makaao_sample.csv input.
    monkeypatch.setattr(
        mod, "makaao_core_name", str(tmp_path / "synthetic_core_fixture.csv")
    )
    monkeypatch.setattr(mod, "OUTPUT_OWL_ENRICHED", str(tmp_path / "makaao_kg_test.rdf"))

    data = mod.load_processed_tables(str(processed))
    assert data["indices"] == {"1", "2", "18"}
    assert data["parents"]["2"] == {"1"}
    assert data["hpo_list"]["1"] == ["HP:0000001", "HP:0000002"]
    assert data["syn_en"]["1"][0][1] == "https://pubmed.ncbi.nlm.nih.gov/123"
    assert data["diseases"]["1"] == [
        ("ORPHA:123", ["https://pubmed.ncbi.nlm.nih.gov/300"]),
        ("C0000009", ["https://pubmed.ncbi.nlm.nih.gov/300"]),
        ("free disease", ["https://pubmed.ncbi.nlm.nih.gov/300"]),
    ]

    graph = mod.init_graph()
    positivity_instances_by_hpo = mod.build_core(
        graph,
        data,
        {"P12345": "Protein X"},
        {"P12345": "https://www.uniprot.org/uniprotkb/P12345"},
        {"C0000001": "Protein X", "C0000009": "Disease CUI"},
        {
            "C0000001": "https://uts.nlm.nih.gov/uts/umls/concept/C0000001",
            "C0000009": "https://uts.nlm.nih.gov/uts/umls/concept/C0000009",
        },
        {"HP:0000001": "All", "HP:0000002": "Abnormality"},
        {"CHEBI:23367": "molecular entity"},
        {"CHEBI:23367": "http://purl.obolibrary.org/obo/CHEBI_23367"},
    )

    aab_instance = mod.MAKAAO["aab_1_instance"]
    positivity_class = mod.MAKAAO["positivity_1"]
    positivity_instance = mod.MAKAAO["positivity_1_instance"]
    hpo_class = URIRef("http://purl.obolibrary.org/obo/HP_0000001")
    assert positivity_instances_by_hpo["HP:0000001"] == (positivity_instance,)
    assert (positivity_instance, RDF.type, positivity_class) in graph
    assert (positivity_instance, RDF.type, hpo_class) not in graph
    assert (positivity_class, SKOS.closeMatch, hpo_class) in graph
    assert (hpo_class, SKOS.closeMatch, positivity_class) in graph
    assert (aab_instance, mod.BIOLINK.biomarker_for, positivity_instance) in graph
    assert (positivity_instance, mod.BIOLINK.has_biomarker, aab_instance) in graph

    ordo = URIRef("http://www.orpha.net/ORDO/Orphanet_123")
    orpha_links = {
        str(ordo): [
            {"HPOId": "HP:0000001", "HPOTerm": "All"},
            {"HPOId": "bad", "HPOTerm": "ignored"},
            {"HPOId": "", "HPOTerm": "ignored"},
        ]
    }
    mod.process_diseases(
        graph,
        data,
        orpha_links,
        positivity_instances_by_hpo,
        {"C0000009": "Disease CUI"},
        {"C0000009": "https://uts.nlm.nih.gov/uts/umls/concept/C0000009"},
        {"HP:0000001": "All"},
    )

    orpha_instance = mod.MAKAAO["orpha_123_instance"]
    assert (orpha_instance, mod.SIO["SIO_001279"], positivity_instance) in graph
    assert (positivity_instance, mod.SIO["SIO_001280"], orpha_instance) in graph
    assert not list(
        graph.triples((mod.MAKAAO["hpo_HP_0000001_instance"], None, None))
    )

    loinc_index = tmp_path / "index_loinc.csv"
    write_csv(
        loinc_index,
        ["aab_id", "loinc_id"],
        [
            {"aab_id": "1", "loinc_id": "LP100-1"},
            {"aab_id": "2", "loinc_id": "LP100-1"},
            {"aab_id": "99", "loinc_id": "LP999-9"},
        ],
    )
    part_tests = tmp_path / "part_tests.json"
    part_tests.write_text(json.dumps({"LP100-1": ["1234-5", "1234-5"]}), encoding="utf-8")
    mod.process_loinc_mappings(
        graph,
        str(loinc_index),
        data["indices"],
        {"LP100-1": "Protein X antibody"},
        {"1234-5": "Protein X test"},
        part_test_json=str(part_tests),
    )
    assert (
        mod.MAKAAO["loinc_1234-5_instance"],
        mod.LOINC_COMPONENT,
        mod.MAKAAO["loinc_LP100-1_instance"],
    ) in graph

    added = mod.add_orpha_umls_close_matches(
        graph,
        [("123", "C0000009")],
        umls_names={"C0000009": "Disease Z"},
        orpha_names={"123": "Disease Z"},
    )
    assert added == 2
    cui_disease_class = mod.MAKAAO["CUI_C0000009"]
    assert (ordo, SKOS.closeMatch, cui_disease_class) in graph
    assert (cui_disease_class, SKOS.closeMatch, ordo) in graph
    # The label-collision audit may report additional equal-label candidates,
    # but it must not remove the explicit ORPHA/UMLS enrichment mapping.
    mod.add_label_collision_close_matches(graph)
    assert (ordo, SKOS.closeMatch, cui_disease_class) in graph
    assert (cui_disease_class, SKOS.closeMatch, ordo) in graph

    mod.append_fair_metadata(graph)
    distribution = mod.set_output_file_metadata(
        graph, Path(mod.OUTPUT_OWL_ENRICHED).name
    )
    dataset = URIRef("http://makaao.inria.fr/kg/")
    distributions = list(graph.objects(dataset, mod.DCAT.distribution))
    assert distributions == [distribution]
    assert list(graph.objects(dataset, mod.VOID.triples))[0].toPython() == len(graph)

    # Add a blank-node restriction to a local class to exercise recursive TBox copying.
    restriction = BNode()
    graph.add((mod.MAKAAO["aab_1"], RDFS.subClassOf, restriction))
    graph.add((restriction, RDF.type, OWL.Restriction))
    graph.add((restriction, OWL.onProperty, mod.BAO_HAS_TARGET))
    graph.add((restriction, OWL.someValuesFrom, mod.MAKAAO.Target))

    tbox = mod.extract_tbox(graph, str(mod.MAKAAO))
    assert (mod.MAKAAO["aab_1"], RDF.type, OWL.Class) in tbox
    assert (ordo, RDFS.subClassOf, mod.MAKAAO.AutoimmuneDisease) in tbox
    assert (restriction, RDF.type, OWL.Restriction) in tbox
    assert not list(tbox.triples((None, mod.DCAT.downloadURL, None)))
    mod.validate_tbox_export(tbox)
    mod.validate_graph_iris(graph, "synthetic assembled KG")
    assert mod.validate_local_cui_labels(graph) >= 4


def test_csv_and_identifier_error_paths(mod, tmp_path):
    duplicate = tmp_path / "duplicate.csv"
    duplicate.write_text(" Index ,index\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate CSV headers"):
        mod.read_csv_rows(duplicate)

    with pytest.raises(ValueError, match="Unexpected separator"):
        mod.to_pubmed_urls("PMID:1\nPMID:2", src_file="x.csv", row=2, col="pmids")

    assert mod.canonical_loinc_code("not-loinc") is None
    assert mod.canonical_hpo_code("not-hpo") is None
    assert mod.canonical_orpha_code("https://example.org/123") is None
    assert mod.hp_to_obo_uri("not-hpo") is None
    with pytest.raises(ValueError, match="Invalid UniProtKB accession"):
        mod.canon_uniprot_id("not accession")
    with pytest.raises(ValueError, match="Invalid ChEBI"):
        mod.canon_chebi_id("not-chebi")

    with pytest.raises(FileNotFoundError):
        mod.read_required_json_object(tmp_path / "missing.json", "test dictionary")
    malformed = tmp_path / "malformed.json"
    malformed.write_text("{", encoding="utf-8")
    with pytest.raises(ValueError, match="Malformed JSON"):
        mod.read_required_json_object(malformed, "test dictionary")
    array = tmp_path / "array.json"
    array.write_text("[]", encoding="utf-8")
    with pytest.raises(ValueError, match="must contain a JSON object"):
        mod.read_required_json_object(array, "test dictionary")


def test_reification_reuses_relation_and_documents(mod):
    graph = mod.init_graph()
    subject = mod.MAKAAO["s"]
    predicate = mod.BAO_HAS_TARGET
    obj = mod.MAKAAO["o"]
    graph.add((subject, predicate, obj))

    mod.add_reified_relation(graph, subject, predicate, obj, "PMID:123")
    mod.add_reified_relation(graph, subject, predicate, obj, "PMID:124")
    mod.add_reified_relation(graph, subject, predicate, obj, "PMID:123")
    mod.add_reified_relation(graph, subject, predicate, obj, "")

    relations = list(graph.subjects(RDF.subject, subject))
    assert len(relations) == 1
    relation = relations[0]
    documents = set(graph.objects(relation, mod.PROV.wasDerivedFrom))
    assert len(documents) == 2
    assert all((document, RDF.type, mod.MAKAAO.Document) in graph for document in documents)
    # PMID provenance is normalized to an external PubMed URL.
    assert all(list(graph.objects(document, RDFS.seeAlso)) for document in documents)


def test_output_metadata_replaces_old_distribution_and_rejects_metadata_tbox(mod):
    graph = Graph()
    dataset = URIRef("http://makaao.inria.fr/kg/")
    old = URIRef("http://example.org/old")
    graph.add((dataset, mod.DCAT.distribution, old))
    graph.add((old, RDF.type, mod.DCAT.Distribution))
    graph.add((old, mod.DCAT.downloadURL, URIRef("http://example.org/old.rdf")))
    graph.add((dataset, mod.VOID.triples, Literal(999)))

    distribution = mod.set_output_file_metadata(graph, "file name.rdf")
    assert distribution != old
    assert not list(graph.triples((old, None, None)))
    assert str(next(graph.objects(distribution, mod.DCAT.downloadURL))).endswith("file%20name.rdf")
    assert next(graph.objects(dataset, mod.VOID.triples)).toPython() == len(graph)

    with pytest.raises(ValueError, match="empty"):
        mod.set_output_file_metadata(graph, "")

    bad_tbox = Graph()
    bad_tbox.add((dataset, RDF.type, mod.DCAT.Dataset))
    with pytest.raises(RuntimeError, match="dataset/distribution/service"):
        mod.validate_tbox_export(bad_tbox)