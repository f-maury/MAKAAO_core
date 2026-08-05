from __future__ import annotations

import importlib.util
import subprocess
import sys
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
    return load_module(
        Path("scripts/04_make_lite_graph_from_makaao-kg.py"), "lite_builder"
    )


def synthetic_graphs(mod):
    data = Graph()
    schema = Graph()

    aab_class = mod.MAKAAO["aab_1"]
    aab_instance = mod.MAKAAO["aab_1_instance"]
    disease_class = mod.ORDO["Orphanet_123"]
    disease_instance = mod.MAKAAO["orpha_123_instance"]
    target_class = mod.MAKAAO["UP_P12345"]
    target_instance = mod.MAKAAO["UP_P12345_instance"]
    phenotype_class = mod.OBO["HP_0030057"]
    phenotype_instance = mod.MAKAAO["hpo_HP_0030057_instance"]
    loinc_part_class = mod.MAKAAO_LOINC["LP12345-6"]
    loinc_part_instance = mod.MAKAAO["loinc_LP12345-6_instance"]
    loinc_term_class = mod.MAKAAO_LOINC["1234-5"]
    loinc_term_instance = mod.MAKAAO["loinc_1234-5_instance"]

    typed = [
        (aab_instance, aab_class, mod.MAKAAO.AutoantibodyPositivity),
        (disease_instance, disease_class, mod.MAKAAO.AutoimmuneDisease),
        (target_instance, target_class, mod.MAKAAO.Target),
        (phenotype_instance, phenotype_class, mod.MAKAAO.AutoantibodyPositivity),
        (loinc_part_instance, loinc_part_class, SKOS.Concept),
        (loinc_term_instance, loinc_term_class, SKOS.Concept),
    ]
    for instance, domain_class, role in typed:
        data.add((instance, RDF.type, domain_class))
        data.add((instance, RDF.type, role))
        data.add((domain_class, RDF.type, OWL.Class))
        schema.add((domain_class, RDF.type, OWL.Class))
        schema.add((domain_class, RDFS.label, Literal(str(domain_class).rsplit("/", 1)[-1])))

    data.add((aab_instance, mod.SIO["SIO_001403"], disease_instance))
    data.add((disease_instance, mod.SIO["SIO_001403"], aab_instance))
    data.add((aab_instance, mod.BIOLINK.biomarker_for, phenotype_instance))
    data.add((disease_instance, mod.SIO["SIO_001279"], phenotype_instance))
    data.add((aab_instance, mod.BAO["BAO_0000211"], target_instance))
    data.add((loinc_term_instance, mod.MAKAAO.hasLoincComponent, loinc_part_instance))
    data.add((aab_instance, SKOS.closeMatch, loinc_part_instance))
    data.add((loinc_part_instance, SKOS.closeMatch, aab_instance))

    # Reification and provenance must not leak into the lite graph.
    relation = mod.MAKAAO["r1"]
    data.add((relation, RDF.type, RDF.Statement))
    data.add((relation, RDF.subject, aab_instance))
    data.add((relation, RDF.predicate, mod.SIO["SIO_001403"]))
    data.add((relation, RDF.object, disease_instance))
    data.add((relation, mod.PROV.wasDerivedFrom, mod.MAKAAO["document_1"]))

    return data, schema


def test_current_version_and_explicit_path_derivation(mod, tmp_path):
    assert mod.SCRIPT_VERSION == "2.0.8"
    input_path = tmp_path / "makaao_kg_1.0.4test.rdf"
    output_path, schema_path = mod.derive_paths_from_explicit_input(input_path)
    assert output_path.name == "makaao_kg_1.0.4test_lite.rdf"
    assert schema_path == input_path

    ontology = tmp_path / "makaao_kg_1.0.4test_ontology.owl"
    ontology.write_text("placeholder", encoding="utf-8")
    _, schema_path = mod.derive_paths_from_explicit_input(input_path)
    assert schema_path == ontology

    reasoning = tmp_path / "reasoning"
    reasoning.mkdir()
    merged = reasoning / "makaao_kg_1.0.4test_curated-tbox-merged.owl"
    merged.write_text("placeholder", encoding="utf-8")
    _, schema_path = mod.derive_paths_from_explicit_input(input_path)
    assert schema_path == merged


def test_build_lite_graph_projects_current_relationships(mod):
    data, schema = synthetic_graphs(mod)
    lite, report = mod.build_lite_graph(data, schema)
    mod.validate_lite_graph(lite)

    assert report["skipped_relationships"] == 0
    assert report["projected_instances"] == 6
    assert not any(mod.is_instance_uri(node) for triple in lite for node in (triple[0], triple[2]))
    assert not list(lite.triples((None, RDF.subject, None)))
    assert not list(lite.triples((None, mod.PROV.wasDerivedFrom, None)))

    aab_class = mod.MAKAAO["aab_1"]
    disease_class = mod.ORDO["Orphanet_123"]
    target_class = mod.MAKAAO["UP_P12345"]
    phenotype_class = mod.OBO["HP_0030057"]
    loinc_part_class = mod.MAKAAO_LOINC["LP12345-6"]
    loinc_term_class = mod.MAKAAO_LOINC["1234-5"]

    assert (aab_class, mod.SIO["SIO_001403"], disease_class) in lite
    assert (aab_class, mod.BIOLINK.biomarker_for, phenotype_class) in lite
    assert (phenotype_class, mod.BIOLINK.has_biomarker, aab_class) in lite
    assert (disease_class, mod.SIO["SIO_001279"], phenotype_class) in lite
    assert (phenotype_class, mod.SIO["SIO_001280"], disease_class) in lite
    assert (aab_class, mod.BAO["BAO_0000211"], target_class) in lite
    assert (target_class, mod.BAO["BAO_0000598"], aab_class) in lite
    assert (loinc_term_class, mod.MAKAAO.hasLoincComponent, loinc_part_class) in lite
    assert (aab_class, SKOS.closeMatch, loinc_part_class) in lite
    assert report["relationship_output_counts"] == {
        str(mod.MAKAAO.hasLoincComponent): 1,
        str(mod.SIO["SIO_001279"]): 1,
        str(mod.SIO["SIO_001280"]): 1,
        str(mod.SIO["SIO_001403"]): 2,
        str(mod.BAO["BAO_0000211"]): 1,
        str(mod.BAO["BAO_0000598"]): 1,
        str(SKOS.closeMatch): 2,
        str(mod.BIOLINK.biomarker_for): 1,
        str(mod.BIOLINK.has_biomarker): 1,
    }


def test_explicit_cli_paths_work_outside_project_checkout(mod, tmp_path):
    data, schema = synthetic_graphs(mod)
    input_path = tmp_path / "input.rdf"
    schema_path = tmp_path / "schema.owl"
    output_path = tmp_path / "output.rdf"
    data.serialize(input_path, format="xml")
    schema.serialize(schema_path, format="xml")

    script = Path("scripts/04_make_lite_graph_from_makaao-kg.py").resolve()
    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--input",
            str(input_path),
            "--schema",
            str(schema_path),
            "--output",
            str(output_path),
        ],
        cwd=tmp_path,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout
    assert "script=2.0.8" in completed.stdout
    output = Graph().parse(output_path)
    mod.validate_lite_graph(output)


def test_lite_validation_rejects_non_class_type(mod):
    graph = Graph()
    graph.add((mod.MAKAAO["x"], RDF.type, mod.MAKAAO.Target))
    with pytest.raises(RuntimeError, match="non-class rdf:type"):
        mod.validate_lite_graph(graph)


def test_path_discovery_and_default_paths(mod, tmp_path, monkeypatch):
    project = tmp_path / "project"
    scripts = project / "scripts"
    data_dir = project / "data"
    kg_dir = project / "kg"
    reasoning = kg_dir / "reasoning"
    scripts.mkdir(parents=True)
    data_dir.mkdir()
    reasoning.mkdir(parents=True)

    builder = scripts / "03_build_kg_from_tables.py"
    builder.write_text('KG_VERSION = "9.9.9"\n', encoding="utf-8")
    script_path = scripts / "04_make_lite_graph_from_makaao-kg.py"
    script_path.write_text("# marker\n", encoding="utf-8")
    canonical = kg_dir / "makaao_kg_9.9.9.rdf"
    canonical.write_text("placeholder", encoding="utf-8")
    merged = reasoning / "makaao_kg_9.9.9_curated-tbox-merged.owl"
    merged.write_text("placeholder", encoding="utf-8")

    assert mod.find_project_dir(scripts) == project
    assert mod.read_builder_kg_version(project, scripts) == "9.9.9"
    monkeypatch.setattr(mod, "__file__", str(script_path))
    input_path, output_path, schema_path = mod.default_paths()
    assert input_path == canonical
    assert output_path == kg_dir / "makaao_kg_9.9.9_lite.rdf"
    assert schema_path == merged

    with pytest.raises(RuntimeError, match="Could not locate"):
        mod.find_project_dir(tmp_path / "unrelated")


def test_canonical_discovery_filters_derived_files_and_rejects_ambiguity(mod, tmp_path):
    kg_dir = tmp_path / "kg"
    kg_dir.mkdir()
    (kg_dir / "makaao_kg_1_lite.rdf").write_text("x", encoding="utf-8")
    (kg_dir / "makaao_kg_1_curated-full.rdf").write_text("x", encoding="utf-8")
    with pytest.raises(FileNotFoundError, match="No canonical"):
        mod.discover_canonical_kg(kg_dir)

    first = kg_dir / "makaao_kg_1.rdf"
    first.write_text("x", encoding="utf-8")
    assert mod.discover_canonical_kg(kg_dir) == first
    second = kg_dir / "makaao_kg_2.rdf"
    second.write_text("x", encoding="utf-8")
    with pytest.raises(RuntimeError, match="More than one canonical"):
        mod.discover_canonical_kg(kg_dir)


def test_projection_fallbacks_and_ambiguity(mod):
    assert mod.fallback_class_from_instance(mod.MAKAAO["aab_18_instance"]) == mod.MAKAAO.Autoantibody
    assert mod.fallback_class_from_instance(mod.MAKAAO["orpha_123_instance"]) == mod.ORDO["Orphanet_123"]
    assert mod.fallback_class_from_instance(mod.MAKAAO["hpo_HP_0000001_instance"]) == mod.OBO["HP_0000001"]
    assert mod.fallback_class_from_instance(mod.MAKAAO["CHEBI_23367_instance"]) == mod.OBO["CHEBI_23367"]
    assert mod.fallback_class_from_instance(mod.MAKAAO["loinc_LP1-1_instance"]) == mod.MAKAAO_LOINC["LP1-1"]
    assert mod.fallback_class_from_instance(mod.MAKAAO["x_instance"]) == mod.MAKAAO["x"]
    assert mod.fallback_class_from_instance(mod.MAKAAO["x"]) is None

    node = mod.MAKAAO["positivity_1_instance"]
    local = mod.MAKAAO["positivity_1"]
    hpo = mod.OBO["HP_0000001"]
    assert mod.choose_projection_class(node, {node: {local, hpo}}) == hpo

    ambiguous = mod.MAKAAO["unknown_instance"]
    with pytest.raises(RuntimeError, match="Ambiguous"):
        mod.choose_projection_class(
            ambiguous,
            {ambiguous: {mod.MAKAAO["A"], mod.MAKAAO["B"]}},
        )
    assert mod.choose_projection_class(Literal("literal"), {}) is None


def test_direct_main_and_atomic_serialization_are_covered(mod, tmp_path, monkeypatch, capsys):
    data, schema = synthetic_graphs(mod)
    input_path = tmp_path / "input.rdf"
    schema_path = tmp_path / "schema.owl"
    output_path = tmp_path / "output.rdf"
    data.serialize(input_path, format="xml")
    schema.serialize(schema_path, format="xml")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "04_make_lite_graph_from_makaao-kg.py",
            "--input",
            str(input_path),
            "--schema",
            str(schema_path),
            "--output",
            str(output_path),
        ],
    )
    mod.main()
    output = Graph().parse(output_path, format="xml")
    mod.validate_lite_graph(output)
    stdout = capsys.readouterr().out
    assert "MAKAAO lite build complete: script=2.0.8" in stdout
    assert "Skipped relationships: 0" in stdout

    with pytest.raises(ValueError, match="different"):
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "04_make_lite_graph_from_makaao-kg.py",
                "--input",
                str(input_path),
                "--schema",
                str(schema_path),
                "--output",
                str(input_path),
            ],
        )
        mod.main()


def test_validation_reports_multiple_invalid_constructs(mod):
    graph = Graph()
    instance = mod.MAKAAO["x_instance"]
    graph.add((instance, RDF.type, mod.MAKAAO.Target))
    graph.add((mod.MAKAAO["x"], mod.PROV.wasDerivedFrom, mod.MAKAAO["document_1"]))
    graph.add((mod.MAKAAO["x"], URIRef("http://example.org/unexpected"), Literal("bad")))
    with pytest.raises(RuntimeError, match="Lite graph validation failed") as exc:
        mod.validate_lite_graph(graph)
    message = str(exc.value)
    assert "unexpected predicate" in message
    assert "instance URI retained" in message
    assert "PROV term retained" in message
