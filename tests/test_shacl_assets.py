from __future__ import annotations

import importlib.util
from pathlib import Path

from rdflib import Graph, Namespace
from rdflib.plugins.sparql.parser import parseQuery


SH = Namespace("http://www.w3.org/ns/shacl#")
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


def test_synthetic_shacl_fixture_is_nonempty_and_serializable(tmp_path):
    module = load_fixture_module()
    fixture = module.build_fixture()
    assert len(fixture) > 0

    output = tmp_path / "makaao-shacl-fixture.ttl"
    fixture.serialize(output, format="turtle")
    reparsed = Graph().parse(output, format="turtle")
    assert len(reparsed) == len(fixture)
