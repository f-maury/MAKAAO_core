from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

from rdflib import Graph, Literal, OWL, RDF, RDFS, URIRef


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _write_processed_tables(processor, frame, output_dir: Path) -> None:
    processor.write_index_name_en(frame, output_dir)
    processor.write_index_hpo_id(frame, output_dir)
    processor.write_index_parent_index(frame, output_dir)
    processor.write_index_syn_en(frame, output_dir)
    processor.write_index_syn_fr(frame, output_dir)
    processor.write_index_cui_source(frame, output_dir)
    processor.write_index_disease_source(frame, output_dir)
    processor.write_index_uniprot_source(frame, output_dir)
    processor.write_index_chebi_source(frame, output_dir)
    processor.write_index_loinc(frame, output_dir)


def test_committed_sample_generates_01_03_and_04_graphs(
    sample_core_path, tmp_path, monkeypatch
):
    """Run scripts 01, core 03 graph construction, and 04 on the public sample.

    Script 02 enrichment and script 03-1 release construction have dedicated
    tests because their complete runs require large or licensed source files
    and an external ROBOT/Java installation.
    """
    processor = load_module(
        Path("scripts/01_process_makaao_core_to_tables.py"),
        "sample_pipeline_processor",
    )
    build = load_module(
        Path("scripts/03_build_kg_from_tables.py"),
        "sample_pipeline_build",
    )
    lite_builder = load_module(
        Path("scripts/04_make_lite_graph_from_makaao-kg.py"),
        "sample_pipeline_lite",
    )

    frame = processor.load_core(sample_core_path)
    assert not frame.empty

    processed = tmp_path / "processed_tables"
    processed.mkdir()
    _write_processed_tables(processor, frame, processed)

    expected = {
        "index_name_en.csv",
        "index_hpo_id.csv",
        "index_parent_index.csv",
        "index_syn_source_en.csv",
        "index_syn_source_fr.csv",
        "index_cui_source.csv",
        "index_disease_source.csv",
        "index_uniprot_source.csv",
        "index_chebi_source.csv",
        "index_loinc.csv",
    }
    missing = sorted(name for name in expected if not (processed / name).is_file())
    assert not missing, f"Script 01 did not create expected tables: {missing}"

    data = build.load_processed_tables(str(processed))
    assert data["indices"], "The committed sample produced no retained indices"

    # Keep this sample test hermetic: labels are deliberately allowed to fall
    # back to identifiers instead of reading a developer's enrichment folder.
    monkeypatch.setattr(build, "CODE_NAMES_CSV", str(tmp_path / "missing-code-names.csv"))
    graph = build.init_graph()
    positivity_instances_by_hpo = build.build_core(
        graph,
        data,
        {},  # UniProt labels: code fallbacks are sufficient for this smoke test.
        {},
        {},  # UMLS labels and URLs.
        {},
        {},  # HPO labels.
        {},  # ChEBI labels and URLs.
        {},
    )
    build.process_diseases(
        graph,
        data,
        {},
        positivity_instances_by_hpo,
        {},
        {},
        {},
    )

    # Exercise the current LOINC interface with Parts taken from the generated
    # sample table. A single synthetic linked Term is sufficient to verify the
    # official loinc:COMPONENT predicate without requiring licensed LOINC data.
    loinc_rows = build.read_csv_rows(processed / "index_loinc.csv")
    sample_parts = sorted(
        {
            part
            for row in loinc_rows
            if (row.get("aab_id") or "").strip() in data["indices"]
            if (part := build.canonical_loinc_code(row.get("loinc_id"), "part"))
        }
    )
    assert sample_parts, "The committed sample no longer exercises LOINC mappings"
    linked_term = "1234-5"
    part_tests_path = tmp_path / "loinc_part_test_dict.json"
    part_tests = {
        part: [linked_term] if part == sample_parts[0] else []
        for part in sample_parts
    }
    part_tests_path.write_text(
        json.dumps(part_tests, indent=2) + "\n", encoding="utf-8"
    )
    build.process_loinc_mappings(
        graph,
        str(processed / "index_loinc.csv"),
        data["indices"],
        {part: f"Sample LOINC Part {part}" for part in sample_parts},
        {linked_term: "Synthetic linked LOINC Term"},
        part_test_json=str(part_tests_path),
    )

    build.append_fair_metadata(graph)
    build.add_label_collision_close_matches(graph)
    canonical_path = tmp_path / "makaao_kg_sample.rdf"
    build.set_output_file_metadata(graph, canonical_path.name)
    build.validate_graph_iris(graph, "sample-driven smoke KG")
    assert build.validate_local_cui_labels(graph) > 0

    graph.serialize(canonical_path, format="xml")
    generated = Graph().parse(canonical_path, format="xml")
    tbox = build.extract_tbox(generated, str(build.MAKAAO))
    projected = build.build_non_reified_graph(generated)
    lite, stats = lite_builder.build_lite_graph(projected, tbox)
    lite_builder.validate_lite_graph(lite)

    biomarker_links = list(
        generated.triples((None, build.BIOLINK.biomarker_for, None))
    )
    assert biomarker_links
    for autoantibody, _, positivity in biomarker_links:
        assert (positivity, build.BIOLINK.has_biomarker, autoantibody) in generated
        assert str(positivity).startswith(str(build.MAKAAO) + "positivity_")
        assert not any(
            str(rdf_type).startswith("http://purl.obolibrary.org/obo/HP_")
            for rdf_type in generated.objects(positivity, RDF.type)
        )

    assert list(generated.triples((None, build.LOINC_COMPONENT, None)))

    assert len(generated) > 0
    assert len(tbox) > 0
    assert len(projected) > 0
    assert len(lite) > 0
    assert stats["projected_instances"] > 0
    assert stats["relationship_output_counts"][str(build.BIOLINK.biomarker_for)] > 0
    assert stats["relationship_output_counts"][str(build.BIOLINK.has_biomarker)] > 0


def test_03_non_reified_projection_is_accepted_by_04_lite_builder():
    build = load_module(Path("scripts/03_build_kg_from_tables.py"), "pipeline_build")
    lite_builder = load_module(
        Path("scripts/04_make_lite_graph_from_makaao-kg.py"), "pipeline_lite"
    )

    graph = build.init_graph()
    left_class = build.MAKAAO["aab_1"]
    right_class = URIRef("http://www.orpha.net/ORDO/Orphanet_123")
    left_instance = build.MAKAAO["aab_1_instance"]
    right_instance = build.MAKAAO["orpha_123_instance"]
    graph.add((left_class, RDF.type, OWL.Class))
    graph.add((right_class, RDF.type, OWL.Class))
    graph.add((left_class, RDFS.label, Literal("Example antibody")))
    graph.add((right_class, RDFS.label, Literal("Example disease")))
    graph.add((left_instance, RDF.type, left_class))
    graph.add((right_instance, RDF.type, right_class))
    graph.add((left_instance, build.SIO["SIO_001403"], right_instance))
    graph.add((right_instance, build.SIO["SIO_001403"], left_instance))

    projected = build.build_non_reified_graph(graph)
    lite, _ = lite_builder.build_lite_graph(projected, projected)
    lite_builder.validate_lite_graph(lite)
    assert (left_class, build.SIO["SIO_001403"], right_class) in lite