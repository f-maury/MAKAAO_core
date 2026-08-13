from __future__ import annotations

import importlib.util
from pathlib import Path

from rdflib import Graph, Namespace, RDF
from rdflib.plugins.sparql.parser import parseQuery
from pyshacl import validate


SH = Namespace("http://www.w3.org/ns/shacl#")
MAKAAO = Namespace("http://makaao.inria.fr/kg/")
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SHAPES_PATH = PROJECT_ROOT / "tests" / "shacl_shapes.ttl"
FIXTURE_SCRIPT = PROJECT_ROOT / "tests" / "build_shacl_fixture.py"


def load_fixture_module():
    spec = importlib.util.spec_from_file_location("build_shacl_fixture", FIXTURE_SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_all_embedded_shacl_sparql_queries_parse():
    shapes = Graph().parse(SHAPES_PATH, format="turtle")
    queries = [str(query) for query in shapes.objects(None, SH.select)]
    assert queries, "No sh:select queries were found in the SHACL shapes"

    for query in queries:
        parseQuery(query)

    required_shapes = {
        MAKAAO.BiomarkerInverseShape,
        MAKAAO.BiomarkerReverseInverseShape,
        MAKAAO.AutoantibodyBiomarkerShape,
        MAKAAO.PositivityHasBiomarkerShape,
        MAKAAO.AutoimmunityRelatedDiseaseShape,
        MAKAAO.TargetShape,
    }
    assert required_shapes <= set(shapes.subjects(RDF.type, SH.NodeShape))


def test_synthetic_shacl_fixture_is_nonempty_and_serializable(tmp_path):
    module = load_fixture_module()
    fixture = module.build_fixture()
    assert len(fixture) > 0

    output = tmp_path / "makaao-shacl-fixture.ttl"
    fixture.serialize(output, format="turtle")
    reparsed = Graph().parse(output, format="turtle")
    assert len(reparsed) == len(fixture)


def validate_fixture(fixture: Graph) -> tuple[bool, str]:
    shapes = Graph().parse(SHAPES_PATH, format="turtle")
    conforms, _, report = validate(
        data_graph=fixture,
        shacl_graph=shapes,
        advanced=True,
        inference="none",
        allow_infos=True,
        allow_warnings=True,
    )
    return bool(conforms), str(report)


def test_synthetic_fixture_conforms_to_the_actual_shapes():
    fixture = load_fixture_module().build_fixture()
    conforms, report = validate_fixture(fixture)
    assert conforms, report


def test_missing_has_biomarker_inverse_is_rejected():
    module = load_fixture_module()
    fixture = module.build_fixture()
    phenotype = module.MAKAAO["test_phenotype_instance"]
    autoantibody = module.MAKAAO["test_autoantibody_instance"]
    fixture.remove((phenotype, module.BIOLINK.has_biomarker, autoantibody))

    conforms, report = validate_fixture(fixture)
    assert not conforms
    assert "has_biomarker" in report