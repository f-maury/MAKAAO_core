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
            {"index": "1", "related_diseases_id": "ORPHA:123", "diseases_pmids": "PMID:301"},
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
        ("ORPHA:123", ["https://pubmed.ncbi.nlm.nih.gov/301"]),
    ]

    graph = mod.init_graph()
    assert (
        mod.MAKAAO.Autoantibody,
        RDFS.subClassOf,
        mod.BIOLINK.ChemicalEntityOrGeneOrGeneProduct,
    ) in graph
    assert not list(graph.triples((mod.MAKAAO.BiomolecularEntity, None, None)))
    assert not list(graph.triples((None, None, mod.MAKAAO.BiomolecularEntity)))
    assert (
        mod.MAKAAO.AutoimmuneDisease,
        RDFS.subClassOf,
        mod.MAKAAO.AutoimmunityRelatedDisease,
    ) in graph

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
    cui_disease_instance = mod.MAKAAO["CUI_C0000009_instance"]
    free_disease_instance = mod.MAKAAO[f"{mod.safe_local_fragment('FREE DISEASE', prefix='disease')}_instance"]
    assert (ordo, RDF.type, OWL.Class) in graph
    assert (ordo, RDFS.subClassOf, mod.MAKAAO.AutoimmunityRelatedDisease) in graph
    assert (
        mod.MAKAAO["CUI_C0000009"],
        RDFS.subClassOf,
        mod.MAKAAO.AutoimmunityRelatedDisease,
    ) in graph
    assert (orpha_instance, RDF.type, mod.MAKAAO.AutoimmunityRelatedDisease) not in graph
    assert (cui_disease_instance, RDF.type, mod.MAKAAO.AutoimmunityRelatedDisease) not in graph
    assert (free_disease_instance, RDF.type, mod.MAKAAO.AutoimmunityRelatedDisease) in graph
    assert (orpha_instance, mod.SIO["SIO_001279"], positivity_instance) in graph
    assert (positivity_instance, mod.SIO["SIO_001280"], orpha_instance) in graph
    assert not list(
        graph.triples((mod.MAKAAO["hpo_HP_0000001_instance"], None, None))
    )

    # Two MAKAAO source rows assert the exact same AAb -> ORDO triple, so they
    # become two Relation occurrences. The ORDO -> HPO source assertion itself
    # still comes from one Orphadata record and must not be duplicated merely
    # because the ORDO disease was encountered twice.
    aab_orpha_relations = {
        relation
        for relation in graph.subjects(RDF.subject, aab_instance)
        if (relation, RDF.predicate, mod.SIO["SIO_001403"]) in graph
        and (relation, RDF.object, orpha_instance) in graph
    }
    assert aab_orpha_relations == {
        mod.deterministic_relation_uri(
            aab_instance, mod.SIO["SIO_001403"], orpha_instance, 1
        ),
        mod.deterministic_relation_uri(
            aab_instance, mod.SIO["SIO_001403"], orpha_instance, 2
        ),
    }

    ordo_hpo_relations = {
        relation
        for relation in graph.subjects(RDF.subject, orpha_instance)
        if (relation, RDF.predicate, mod.SIO["SIO_001279"]) in graph
        and (relation, RDF.object, positivity_instance) in graph
    }
    assert ordo_hpo_relations == {
        mod.deterministic_relation_uri(
            orpha_instance, mod.SIO["SIO_001279"], positivity_instance, 1
        )
    }

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
    loinc_term_instance = mod.MAKAAO["loinc_1234-5_instance"]
    loinc_part_instance = mod.MAKAAO["loinc_LP100-1_instance"]
    assert (
        loinc_term_instance,
        mod.LOINC_COMPONENT,
        loinc_part_instance,
    ) in graph
    loinc_component_relations = {
        relation
        for relation in graph.subjects(RDF.subject, loinc_term_instance)
        if (relation, RDF.predicate, mod.LOINC_COMPONENT) in graph
        and (relation, RDF.object, loinc_part_instance) in graph
    }
    assert loinc_component_relations == {
        mod.deterministic_relation_uri(
            loinc_term_instance, mod.LOINC_COMPONENT, loinc_part_instance, 1
        )
    }

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
    assert (ordo, RDF.type, OWL.Class) in tbox
    assert (ordo, RDFS.subClassOf, mod.MAKAAO.AutoimmunityRelatedDisease) in tbox
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


def test_reification_preserves_occurrences_and_reuses_documents(mod):
    graph = mod.init_graph()
    subject = mod.MAKAAO["s"]
    predicate = mod.BAO_HAS_TARGET
    obj = mod.MAKAAO["o"]
    graph.add((subject, predicate, obj))

    # One source occurrence may have several provenance documents; this remains
    # one Relation occurrence. A second source occurrence of the exact same S/P/O
    # gets the next local occurrence number.
    relation_1 = mod.add_reified_relation(
        graph, subject, predicate, obj, ["PMID:123", "PMID:124"]
    )
    relation_2 = mod.add_reified_relation(
        graph, subject, predicate, obj, "PMID:123"
    )
    relation_3 = mod.add_reified_relation(
        graph, subject, predicate, obj, "PMID:125"
    )
    missing = mod.add_reified_relation(graph, subject, predicate, obj, "")

    assert missing is None
    assert relation_1 == mod.MAKAAO["rel_s_BAO~u0000211_o_1"]
    assert relation_2 == mod.MAKAAO["rel_s_BAO~u0000211_o_2"]
    assert relation_3 == mod.MAKAAO["rel_s_BAO~u0000211_o_3"]

    relations = set(graph.subjects(RDF.subject, subject))
    assert relations == {relation_1, relation_2, relation_3}
    documents_1 = set(graph.objects(relation_1, mod.PROV.wasDerivedFrom))
    documents_2 = set(graph.objects(relation_2, mod.PROV.wasDerivedFrom))
    documents_3 = set(graph.objects(relation_3, mod.PROV.wasDerivedFrom))
    assert len(documents_1) == 2
    assert len(documents_2) == len(documents_3) == 1
    assert documents_2 < documents_1
    documents = documents_1 | documents_2 | documents_3
    assert len(documents) == 3
    assert all(
        (document, RDF.type, mod.MAKAAO.Document) in graph
        for document in documents
    )
    # PMID provenance is normalized to an external PubMed URL.
    assert all(list(graph.objects(document, RDFS.seeAlso)) for document in documents)

    # Formerly ambiguous component boundaries must not produce the same URI.
    collision_a = mod.deterministic_relation_uri(
        mod.MAKAAO["a"],
        predicate,
        mod.MAKAAO["b_BAO_0000211_c"],
        1,
    )
    collision_b = mod.deterministic_relation_uri(
        mod.MAKAAO["a_BAO_0000211_b"],
        predicate,
        mod.MAKAAO["c"],
        1,
    )
    assert collision_a != collision_b


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


def test_strict_autoimmune_dictionary_and_cui_ordo_propagation(mod, tmp_path):
    graph = mod.init_graph()

    ordo_from_cui = URIRef("http://www.orpha.net/ORDO/Orphanet_123")
    cui_direct = mod.MAKAAO["CUI_C0000001"]
    graph.add((ordo_from_cui, RDF.type, OWL.Class))
    graph.add((ordo_from_cui, RDFS.subClassOf, mod.MAKAAO.AutoimmunityRelatedDisease))
    graph.add((cui_direct, RDF.type, OWL.Class))
    graph.add((cui_direct, RDFS.subClassOf, mod.MAKAAO.CUI))
    graph.add((cui_direct, RDFS.subClassOf, mod.MAKAAO.AutoimmunityRelatedDisease))
    graph.add((ordo_from_cui, SKOS.closeMatch, cui_direct))
    graph.add((cui_direct, SKOS.closeMatch, ordo_from_cui))

    ordo_direct = URIRef("http://www.orpha.net/ORDO/Orphanet_456")
    cui_from_orpha = mod.MAKAAO["CUI_C0000002"]
    graph.add((ordo_direct, RDF.type, OWL.Class))
    graph.add((ordo_direct, RDFS.subClassOf, mod.MAKAAO.AutoimmunityRelatedDisease))
    graph.add((cui_from_orpha, RDF.type, OWL.Class))
    graph.add((cui_from_orpha, RDFS.subClassOf, mod.MAKAAO.CUI))
    graph.add((cui_from_orpha, RDFS.subClassOf, mod.MAKAAO.AutoimmunityRelatedDisease))
    graph.add((ordo_direct, SKOS.closeMatch, cui_from_orpha))
    graph.add((cui_from_orpha, SKOS.closeMatch, ordo_direct))

    # A closeMatch to a Target-only CUI is deliberately not an eligible disease edge.
    target_cui = mod.MAKAAO["CUI_C0000003"]
    target_instance = mod.MAKAAO["CUI_C0000003_instance"]
    graph.add((target_cui, RDF.type, OWL.Class))
    graph.add((target_cui, RDFS.subClassOf, mod.MAKAAO.CUI))
    graph.add((target_instance, RDF.type, target_cui))
    graph.add((target_instance, RDF.type, mod.MAKAAO.Target))
    graph.add((ordo_direct, SKOS.closeMatch, target_cui))
    graph.add((target_cui, SKOS.closeMatch, ordo_direct))

    dictionary = tmp_path / "dico_MAI_strict.json"
    dictionary.write_text(
        json.dumps(
            {
                "CUI:C0000001": ["direct CUI disease"],
                "ORPHA:456": ["direct ORPHA disease"],
                "ORPHA:999": ["not present in this synthetic graph"],
            }
        ),
        encoding="utf-8",
    )
    strict_cuis, strict_orphas = mod.read_strict_autoimmune_disease_ids(dictionary)
    assert strict_cuis == {"C0000001"}
    assert strict_orphas == {"456", "999"}

    result = mod.apply_strict_autoimmune_classification(
        graph, strict_cuis, strict_orphas
    )

    assert (cui_direct, RDFS.subClassOf, mod.MAKAAO.AutoimmuneDisease) in graph
    assert (ordo_from_cui, RDFS.subClassOf, mod.MAKAAO.AutoimmuneDisease) in graph
    assert (ordo_direct, RDFS.subClassOf, mod.MAKAAO.AutoimmuneDisease) in graph
    assert (cui_from_orpha, RDFS.subClassOf, mod.MAKAAO.AutoimmuneDisease) in graph
    assert (target_cui, RDFS.subClassOf, mod.MAKAAO.AutoimmuneDisease) not in graph

    assert len(result["direct_classes"]) == 2
    assert len(result["propagated_classes"]) == 2
    assert len(result["strict_classes"]) == 4
    assert result["unmatched_orphas"] == {"999"}
    assert not result["unmatched_cuis"]

    malformed = tmp_path / "bad_strict_dictionary.json"
    malformed.write_text(json.dumps({"SNOMED:123": ["bad key"]}), encoding="utf-8")
    with pytest.raises(ValueError, match="expected CUI:... or ORPHA:..."):
        mod.read_strict_autoimmune_disease_ids(malformed)



def test_release_manifest_records_strict_autoimmune_curation(mod, tmp_path):
    stage_root = tmp_path / "stage"
    stage_reasoning = stage_root / "reasoning"
    stage_reasoning.mkdir(parents=True)
    stage_kg = stage_root / "makaao_kg_1.0.5.rdf"
    stage_tbox = stage_root / "makaao_kg_1.0.5_ontology.owl"
    root_catalog = stage_root / "catalog-imports.xml"
    stage_kg.write_text("kg", encoding="utf-8")
    stage_tbox.write_text("tbox", encoding="utf-8")
    root_catalog.write_text("<catalog/>", encoding="utf-8")

    manifest_path = stage_reasoning / "reasoning-manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "canonical_inputs": {
                    "kg": {"sha256": mod._file_sha256(stage_kg)},
                    "tbox": {"sha256": mod._file_sha256(stage_tbox)},
                }
            }
        ),
        encoding="utf-8",
    )
    curation = {
        "source": "data/dico_MAI_strict.json",
        "sha256": "abc123",
        "entries": 177,
        "direct_classes": 104,
        "propagated_classes": 82,
        "total_classes": 186,
        "unmatched_entries": 73,
    }

    mod.finalize_release_manifest_and_checksums(
        stage_root,
        stage_reasoning,
        stage_kg,
        stage_tbox,
        root_catalog,
        curation,
    )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["strict_autoimmune_curation"] == curation
    assert (stage_reasoning / "SHA256SUMS").is_file()