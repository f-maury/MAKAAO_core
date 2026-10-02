from __future__ import annotations

import ast
import importlib.util
from pathlib import Path

from rdflib import BNode, Graph, Namespace, OWL, RDF, RDFS
from rdflib.plugins.sparql.parser import parseQuery
from pyshacl import validate


SH = Namespace("http://www.w3.org/ns/shacl#")
MAKAAO = Namespace("http://makaao.inria.fr/kg/")
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SHAPES_PATH = PROJECT_ROOT / "tests" / "shacl_shapes.ttl"
FIXTURE_SCRIPT = PROJECT_ROOT / "tests" / "build_shacl_fixture.py"
BUILDER_SCRIPT = PROJECT_ROOT / "scripts" / "03_build_kg_from_tables.py"
KG_DIR = PROJECT_ROOT / "kg"


def load_fixture_module():
    spec = importlib.util.spec_from_file_location("build_shacl_fixture", FIXTURE_SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_current_kg_version() -> str:
    """Read KG_VERSION from the current graph builder without importing it."""
    assert BUILDER_SCRIPT.is_file(), f"MAKAAO graph builder is missing: {BUILDER_SCRIPT}"
    tree = ast.parse(
        BUILDER_SCRIPT.read_text(encoding="utf-8"),
        filename=str(BUILDER_SCRIPT),
    )
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if not any(
            isinstance(target, ast.Name) and target.id == "KG_VERSION"
            for target in targets
        ):
            continue
        value = node.value
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            return value.value
    raise AssertionError(f"Could not read KG_VERSION from {BUILDER_SCRIPT}")


def current_canonical_kg_path() -> Path:
    version = read_current_kg_version()
    return KG_DIR / f"makaao_kg_{version}.rdf"


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
        MAKAAO.NoAnonymousClassExpressionShape,
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


def validate_graph(data_graph: Graph) -> tuple[bool, str]:
    shapes = Graph().parse(SHAPES_PATH, format="turtle")
    conforms, _, report = validate(
        data_graph=data_graph,
        shacl_graph=shapes,
        advanced=True,
        inference="none",
        allow_infos=True,
        allow_warnings=True,
    )
    return bool(conforms), str(report)


def test_synthetic_fixture_conforms_to_the_actual_shapes():
    fixture = load_fixture_module().build_fixture()
    conforms, report = validate_graph(fixture)
    assert conforms, report


def test_missing_has_biomarker_inverse_is_rejected():
    module = load_fixture_module()
    fixture = module.build_fixture()
    phenotype = module.MAKAAO["test_phenotype_instance"]
    autoantibody = module.MAKAAO["test_autoantibody_instance"]
    fixture.remove((phenotype, module.BIOLINK.has_biomarker, autoantibody))

    conforms, report = validate_graph(fixture)
    assert not conforms
    assert "has_biomarker" in report


def test_anonymous_class_expressions_are_rejected():
    """Restrictions and unions must never appear in the MAKAAO KG."""
    module = load_fixture_module()
    positivity_class = module.MAKAAO["positivity_test"]

    with_restriction = module.build_fixture()
    restriction = BNode()
    with_restriction.add((module.MAKAAO["TestAutoantibody"], RDFS.subClassOf, restriction))
    with_restriction.add((restriction, RDF.type, OWL.Restriction))
    with_restriction.add((restriction, OWL.onProperty, module.BIOLINK.biomarker_for))
    with_restriction.add((restriction, OWL.someValuesFrom, positivity_class))
    conforms, report = validate_graph(with_restriction)
    assert not conforms
    assert "anonymous OWL class expressions" in report

    with_union = module.build_fixture()
    union, members = BNode(), BNode()
    with_union.add((union, RDF.type, OWL.Class))
    with_union.add((union, OWL.unionOf, members))
    with_union.add((members, RDF.first, positivity_class))
    with_union.add((members, RDF.rest, RDF.nil))
    conforms, report = validate_graph(with_union)
    assert not conforms
    assert "anonymous OWL class expressions" in report


def test_current_canonical_kg_conforms_to_shacl():
    """Run the complete SHACL suite against the current production KG artifact."""
    kg_path = current_canonical_kg_path()

    if not kg_path.is_file():
        pytest.skip(f"Canonical KG not available: {kg_path}")


    graph = Graph()
    try:
        graph.parse(kg_path, format="xml")
        assert len(graph) > 0, f"Canonical KG contains no RDF triples: {kg_path}"
        conforms, report = validate_graph(graph)
    finally:
        graph.close()

    assert conforms, (
        f"Canonical MAKAAO KG {kg_path.name} failed SHACL validation.\n"
        f"{report}"
    )