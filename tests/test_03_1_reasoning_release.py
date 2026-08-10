import sys
import importlib.util
from pathlib import Path

import pytest
from rdflib import BNode, Graph, Literal, OWL, RDF, RDFS, URIRef


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return load_module(
        Path("scripts/03-1_build_reasoning_release.py"), "reasoning_release"
    )


def test_signature_collects_terms_and_predicates(mod):
    graph = Graph()
    subject = URIRef("http://example.org/vocab/A")
    predicate = URIRef("http://example.org/vocab/p")
    obj = URIRef("http://example.org/other/B")
    graph.add((subject, predicate, obj))

    terms, predicates = mod.collect_signature(
        graph, ("http://example.org/vocab/",)
    )
    assert terms == {subject, predicate}
    assert predicates == {predicate}


def test_expected_import_validation_accepts_exact_nine_module_set(mod):
    graph = Graph()
    ontology = URIRef("http://example.org/ontology")
    graph.add((ontology, RDF.type, OWL.Ontology))
    for spec in mod.VOCABULARIES:
        graph.add((ontology, OWL.imports, spec.module_iri))

    expected_imports = {
        "http://makaao.inria.fr/imports/bao",
        "http://makaao.inria.fr/imports/biolink",
        "http://makaao.inria.fr/imports/chebi",
        "http://makaao.inria.fr/imports/hpo",
        "http://makaao.inria.fr/imports/ordo",
        "http://makaao.inria.fr/imports/prov",
        "http://makaao.inria.fr/imports/sio",
        "http://makaao.inria.fr/imports/uniprot-core",
        "http://www.w3.org/TR/skos-reference/skos-owl1-dl.rdf",
    }
    result = mod.validate_expected_imports(graph)
    assert set(result["expected"]) == expected_imports
    assert set(result["actual"]) == expected_imports

    graph.add((ontology, OWL.imports, URIRef("http://example.org/unexpected")))
    with pytest.raises(RuntimeError, match="unsupported imports"):
        mod.validate_expected_imports(graph)


def test_remove_generic_reification_preserves_ordinary_assertion(mod):
    graph = Graph()
    subject = URIRef("http://example.org/s")
    predicate = URIRef("http://example.org/p")
    obj = URIRef("http://example.org/o")
    relation = URIRef("http://makaao.inria.fr/kg/r1")
    document = URIRef("http://makaao.inria.fr/kg/document_1")

    graph.add((subject, predicate, obj))
    graph.add((relation, RDF.type, RDF.Statement))
    graph.add((relation, RDF.type, mod.MAKAAO.Relation))
    graph.add((relation, RDF.subject, subject))
    graph.add((relation, RDF.predicate, predicate))
    graph.add((relation, RDF.object, obj))
    graph.add((document, RDF.type, mod.MAKAAO.Document))

    cleaned = mod.remove_generic_reification(graph)
    assert (subject, predicate, obj) in cleaned
    assert not list(cleaned.triples((relation, None, None)))
    assert not list(cleaned.triples((document, None, None)))
    assert mod.graph_postconditions(cleaned, require_no_imports=True) == {
        "rdf_statement_instances": 0,
        "makaao_relation_instances": 0,
        "makaao_document_instances": 0,
        "rdf_reification_predicate_triples": 0,
        "owl_imports": 0,
    }


def test_graph_postconditions_reject_reification_and_imports(mod):
    graph = Graph()
    graph.add((URIRef("http://example.org/r"), RDF.type, RDF.Statement))
    with pytest.raises(RuntimeError, match="Reification/provenance"):
        mod.graph_postconditions(graph, require_no_imports=True)

    graph = Graph()
    ontology = URIRef("http://example.org/ontology")
    graph.add((ontology, RDF.type, OWL.Ontology))
    graph.add((ontology, OWL.imports, URIRef("http://example.org/import")))
    with pytest.raises(RuntimeError, match="owl:imports"):
        mod.graph_postconditions(graph, require_no_imports=True)


def test_bao_validation_requires_both_object_properties(mod):
    source = Graph()
    module = Graph()
    for prop in (mod.BAO_HAS_TARGET, mod.BAO_IS_TARGET_FOR):
        source.add((prop, RDF.type, OWL.ObjectProperty))
        module.add((prop, RDF.type, OWL.ObjectProperty))
    mod.validate_bao_target_properties(source, module)

    module.remove((mod.BAO_IS_TARGET_FOR, RDF.type, OWL.ObjectProperty))
    with pytest.raises(RuntimeError, match="BAO module validation failed"):
        mod.validate_bao_target_properties(source, module)


def test_biolink_biomarker_axioms_are_preserved_and_validated(mod):
    source = Graph()
    module = Graph()
    autoantibody = URIRef("http://makaao.inria.fr/kg/Autoantibody")
    phenotype = mod.BIOLINK.PhenotypicFeature
    axioms = {
        mod.BIOLINK.biomarker_for: (autoantibody, phenotype),
        mod.BIOLINK.has_biomarker: (phenotype, autoantibody),
    }
    for prop, (domain, range_) in axioms.items():
        for graph in (source, module):
            graph.add((prop, RDF.type, OWL.ObjectProperty))
            graph.add((prop, RDFS.domain, domain))
            graph.add((prop, RDFS.range, range_))
    for graph in (source, module):
        graph.add(
            (
                mod.BIOLINK.biomarker_for,
                OWL.inverseOf,
                mod.BIOLINK.has_biomarker,
            )
        )

    mod.validate_biolink_biomarker_axioms(source, module)

    module.remove(
        (
            mod.BIOLINK.biomarker_for,
            OWL.inverseOf,
            mod.BIOLINK.has_biomarker,
        )
    )
    with pytest.raises(RuntimeError, match="inverse axioms differ"):
        mod.validate_biolink_biomarker_axioms(source, module)

    anonymous_source = Graph()
    anonymous_source += source
    anonymous_source.set(
        (mod.BIOLINK.biomarker_for, RDFS.range, BNode())
    )
    with pytest.raises(RuntimeError, match="anonymous"):
        mod.validate_biolink_biomarker_axioms(anonymous_source, source)


def test_bao_module_extraction_preserves_inverse_axiom_and_audit_counts(mod, tmp_path):
    bao_spec = next(spec for spec in mod.VOCABULARIES if spec.name == "bao")
    source = Graph()
    source_ontology = URIRef("http://example.org/bao-source")
    source.add((source_ontology, RDF.type, OWL.Ontology))
    source.add((mod.BAO_HAS_TARGET, RDF.type, OWL.ObjectProperty))
    source.add((mod.BAO_IS_TARGET_FOR, RDF.type, OWL.ObjectProperty))
    source.add((mod.BAO_HAS_TARGET, RDFS.label, Literal("has target")))
    source.add((mod.BAO_IS_TARGET_FOR, RDFS.label, Literal("is target of")))
    source.add((mod.BAO_HAS_TARGET, OWL.inverseOf, mod.BAO_IS_TARGET_FOR))
    source_path = tmp_path / "bao-source.owl"
    source.serialize(source_path, format="xml")

    signature = Graph()
    antibody = URIRef("http://example.org/antibody")
    target = URIRef("http://example.org/target")
    signature.add((antibody, mod.BAO_HAS_TARGET, target))
    signature.add((target, mod.BAO_IS_TARGET_FOR, antibody))

    output_path = tmp_path / "bao-module.owl"
    result = mod.build_reasoning_module(
        bao_spec, str(source_path), signature, output_path, "test"
    )
    module = Graph().parse(output_path, format="xml")

    assert len(module) == result.triples == 12
    assert result.used_terms == 2
    assert result.used_predicates == 2
    assert result.selected_classes == 0
    assert result.selected_object_properties == 2
    assert result.retained_axiom_counts["object_property_declarations"] == 2
    assert result.retained_axiom_counts["inverse_property_axioms"] == 1
    assert (mod.BAO_HAS_TARGET, OWL.inverseOf, mod.BAO_IS_TARGET_FOR) in module


def _minimal_source_files(mod, tmp_path):
    sources = {}
    for spec in mod.VOCABULARIES:
        graph = Graph()
        ontology = (
            mod.SKOS_ONTOLOGY_IRI
            if spec.name == "skos"
            else URIRef(f"http://example.org/source/{spec.name}")
        )
        graph.add((ontology, RDF.type, OWL.Ontology))
        path = tmp_path / f"{spec.name}.owl"
        graph.serialize(path, format="xml")
        sources[spec.source_key] = str(path)
    return sources


def test_complete_reasoning_release_orchestration_without_external_robot(
    mod, tmp_path, monkeypatch
):
    sources = _minimal_source_files(mod, tmp_path)
    skos_path = sources.pop("skos")

    kg = Graph()
    tbox = Graph()
    ontology = URIRef("http://makaao.inria.fr/kg/")
    tbox.add((ontology, RDF.type, OWL.Ontology))
    for spec in mod.VOCABULARIES:
        tbox.add((ontology, OWL.imports, spec.module_iri))
    local_class = URIRef("http://makaao.inria.fr/kg/LocalClass")
    local_instance = URIRef("http://makaao.inria.fr/kg/local_instance")
    tbox.add((local_class, RDF.type, OWL.Class))
    kg.add((local_class, RDF.type, OWL.Class))
    kg.add((local_instance, RDF.type, local_class))

    canonical_kg = tmp_path / "canonical.rdf"
    canonical_tbox = tmp_path / "canonical.owl"
    kg.serialize(canonical_kg, format="xml")
    tbox.serialize(canonical_tbox, format="xml")

    def fake_validate_profile(
        robot_command,
        input_path,
        report_path,
        work_dir,
        timeout_seconds,
        subprocess_env=None,
    ):
        assert robot_command == ["robot"]
        assert timeout_seconds == 30
        report_path.write_text(
            "OWL 2 DL Profile Report: [Ontology and imports closure in profile]\n",
            encoding="utf-8",
        )

    def fake_run_reasoner(
        robot_command,
        input_path,
        output_path,
        work_dir,
        timeout_seconds,
        reasoner_name,
        subprocess_env=None,
    ):
        reasoned = Graph().parse(input_path, format="xml")
        reasoned.add(
            (
                URIRef("http://example.org/inferred"),
                RDF.type,
                OWL.Class,
            )
        )
        reasoned.serialize(output_path, format="xml")

    monkeypatch.setattr(mod, "validate_profile", fake_validate_profile)
    monkeypatch.setattr(mod, "run_reasoner", fake_run_reasoner)

    output = tmp_path / "reasoning"
    manifest = mod.build_reasoning_release(
        kg_graph=kg,
        tbox_graph=tbox,
        non_reified_graph=kg,
        output_dir=output,
        source_files=sources,
        skos_source_file=skos_path,
        kg_version="test",
        robot_command=["robot"],
        reasoner_name="HermiT",
        timeout_seconds=30,
        require_robot=True,
        canonical_kg_path=canonical_kg,
        canonical_tbox_path=canonical_tbox,
        subprocess_env={"TEST": "1"},
    )

    assert manifest["status"] == "PASSED_OWL2_DL_AND_HERMIT"
    assert len(manifest["modules"]) == 9
    assert manifest["design"]["asserted_dl_triples_preserved_in_reasoned_graph"] is True
    assert manifest["canonical_inputs"]["kg"]["sha256"] == mod.file_sha256(canonical_kg)
    assert manifest["outputs"]["reasoned_dl_kg"].endswith("-reasoned.owl")
    assert manifest["kg_reasoning"]["assertion_preservation"][
        "asserted_triples_missing_after_restore"
    ] == 0

    assert (output / "catalog-v001.xml").is_file()
    assert (output / "reasoning-manifest.json").is_file()
    assert (output / "README.txt").is_file()
    assert (output / "EXTRACTION_POLICY.txt").is_file()
    assert (output / "SHA256SUMS").is_file()
    assert len(list((output / "modules").glob("*.owl"))) == 9
    assert len(list((output / "reports").glob("*-selection.json"))) == 9
    assert len(list((output / "reports").glob("*-profile.txt"))) == 12


def test_reasoning_release_rejects_bad_timeout_and_missing_robot(mod, tmp_path):
    sources = _minimal_source_files(mod, tmp_path)
    skos_path = sources.pop("skos")
    tbox = Graph()
    ontology = URIRef("http://makaao.inria.fr/kg/")
    tbox.add((ontology, RDF.type, OWL.Ontology))
    for spec in mod.VOCABULARIES:
        tbox.add((ontology, OWL.imports, spec.module_iri))

    with pytest.raises(ValueError, match="greater than zero"):
        mod.build_reasoning_release(
            kg_graph=Graph(),
            tbox_graph=tbox,
            non_reified_graph=Graph(),
            output_dir=tmp_path / "bad-timeout",
            source_files=sources,
            skos_source_file=skos_path,
            kg_version="test",
            robot_command=None,
            timeout_seconds=0,
            require_robot=False,
        )

    with pytest.raises(RuntimeError, match="ROBOT is required"):
        mod.build_reasoning_release(
            kg_graph=Graph(),
            tbox_graph=tbox,
            non_reified_graph=Graph(),
            output_dir=tmp_path / "missing-robot",
            source_files=sources,
            skos_source_file=skos_path,
            kg_version="test",
            robot_command=None,
            timeout_seconds=30,
            require_robot=True,
        )