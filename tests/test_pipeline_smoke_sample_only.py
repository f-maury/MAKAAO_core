from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path

import pytest
from rdflib import Graph, Literal, Namespace, OWL, RDF, RDFS, URIRef


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


def test_committed_sample_runs_through_01_03_and_04(sample_core_path, tmp_path):
    """Exercise the public sample without private data or ontology downloads."""
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

    graph = build.init_graph()
    build.build_core(
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
    build.process_diseases(graph, data, {}, {}, {}, {})
    build.validate_graph_iris(graph, "sample-driven smoke KG")
    assert build.validate_local_cui_labels(graph) >= 0

    projected = build.build_non_reified_graph(graph)
    lite, stats = lite_builder.build_lite_graph(projected, projected)
    lite_builder.validate_lite_graph(lite)

    assert len(graph) > 0
    assert len(projected) > 0
    assert len(lite) > 0
    assert stats["projected_instances"] >= 0


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


def test_checked_in_sample_graph_is_parseable_when_present():
    path = Path("kg/makaao_kg_sample.rdf")
    if not path.is_file():
        pytest.skip("kg/makaao_kg_sample.rdf is not present")
    graph = Graph().parse(path)
    assert len(graph) > 0
    dcterms = Namespace("http://purl.org/dc/terms/")
    versions = {str(value) for value in graph.objects(None, dcterms.hasVersion)}

    # Legacy sample fixtures may predate dataset-version metadata. When metadata
    # is present, it must agree with the current builder rather than a hard-coded
    # version string.
    if versions:
        script_path = Path("scripts/03_build_kg_from_tables.py")
        tree = ast.parse(
            script_path.read_text(encoding="utf-8"),
            filename=str(script_path),
        )
        expected = None
        for node in tree.body:
            if not isinstance(node, ast.Assign):
                continue
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "SCRIPT_VERSION":
                    expected = ast.literal_eval(node.value)
                    break
            if expected is not None:
                break

        assert expected is not None, "SCRIPT_VERSION is missing from script 03"
        assert expected in versions, (
            f"{path} is stale: expected script version {expected!r}; "
            f"found {sorted(versions)!r}"
        )
