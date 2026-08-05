from __future__ import annotations

import argparse
from pathlib import Path

from rdflib import Graph, Literal, Namespace, URIRef
from rdflib.namespace import OWL, RDF, RDFS, SKOS

MAKAAO = Namespace("http://makaao.inria.fr/kg/")
SIO = Namespace("http://semanticscience.org/resource/")
BAO = Namespace("http://www.bioassayontology.org/bao#")
BIOLINK = Namespace("https://w3id.org/biolink/vocab/")
PROV = Namespace("http://www.w3.org/ns/prov#")


def build_fixture() -> Graph:
    """Build a small, license-safe graph that exercises every SHACL shape."""
    graph = Graph()
    graph.bind("makaao", MAKAAO)
    graph.bind("sio", SIO)
    graph.bind("bao", BAO)
    graph.bind("biolink", BIOLINK)
    graph.bind("prov", PROV)
    graph.bind("skos", SKOS)

    autoantibody_class = MAKAAO["TestAutoantibody"]
    autoantibody = MAKAAO["test_autoantibody_instance"]
    disease = MAKAAO["test_disease_instance"]
    phenotype = MAKAAO["test_phenotype_instance"]
    target = MAKAAO["test_target_instance"]
    local_umls = MAKAAO["CUI_C0000001"]
    external_match = URIRef("http://purl.uniprot.org/uniprot/P00001")
    document = MAKAAO["test_document"]
    relation = MAKAAO["test_relation"]

    graph.add((autoantibody_class, RDF.type, OWL.Class))
    graph.add((autoantibody_class, RDFS.subClassOf, MAKAAO.Autoantibody))
    graph.add((autoantibody_class, RDFS.label, Literal("Test autoantibody class")))

    graph.add((autoantibody, RDF.type, autoantibody_class))
    graph.add((autoantibody, RDFS.label, Literal("Test autoantibody")))
    graph.add((autoantibody, BIOLINK.biomarker_for, phenotype))
    graph.add((autoantibody, BAO.BAO_0000211, target))

    graph.add((disease, RDF.type, MAKAAO.AutoimmuneDisease))
    graph.add((disease, RDFS.label, Literal("Test autoimmune disease")))
    graph.add((disease, SIO.SIO_001403, autoantibody))
    graph.add((autoantibody, SIO.SIO_001403, disease))

    graph.add((autoantibody, SIO.SIO_001279, phenotype))
    graph.add((phenotype, SIO.SIO_001280, autoantibody))
    graph.add((phenotype, RDFS.label, Literal("Test phenotype")))

    graph.add((target, RDF.type, MAKAAO.Target))
    graph.add((target, RDFS.label, Literal("Test target")))
    graph.add((target, BAO.BAO_0000598, autoantibody))

    graph.add((local_umls, RDF.type, OWL.Class))
    graph.add((local_umls, RDFS.label, Literal("Test UMLS concept")))
    graph.add((external_match, RDF.type, OWL.Class))
    graph.add((external_match, RDFS.label, Literal("Test UMLS concept")))
    graph.add((local_umls, SKOS.closeMatch, external_match))
    graph.add((external_match, SKOS.closeMatch, local_umls))

    graph.add((document, RDF.type, MAKAAO.Document))
    graph.add((document, RDFS.label, Literal("Synthetic test source")))
    graph.add((document, RDFS.seeAlso, URIRef("https://example.org/makaao-ci-fixture")))

    graph.add((relation, RDF.type, MAKAAO.Relation))
    graph.add((relation, RDF.type, RDF.Statement))
    graph.add((relation, RDF.subject, autoantibody))
    graph.add((relation, RDF.predicate, BIOLINK.biomarker_for))
    graph.add((relation, RDF.object, phenotype))
    graph.add((relation, PROV.wasDerivedFrom, document))

    assert_fixture_exercises_shapes(graph)
    return graph


def assert_fixture_exercises_shapes(graph: Graph) -> None:
    """Prevent the SHACL run from becoming a vacuous empty-graph check."""
    required = {
        "relation": (None, RDF.type, MAKAAO.Relation),
        "document": (None, RDF.type, MAKAAO.Document),
        "closeMatch": (None, SKOS.closeMatch, None),
        "association": (None, SIO.SIO_001403, None),
        "has phenotype": (None, SIO.SIO_001279, None),
        "is phenotype of": (None, SIO.SIO_001280, None),
        "has target": (None, BAO.BAO_0000211, None),
        "is target of": (None, BAO.BAO_0000598, None),
        "target instance": (None, RDF.type, MAKAAO.Target),
        "disease instance": (None, RDF.type, MAKAAO.AutoimmuneDisease),
        "autoantibody subclass": (None, RDFS.subClassOf, MAKAAO.Autoantibody),
    }
    missing = [name for name, triple in required.items() if triple not in graph]
    if missing:
        raise AssertionError(f"Synthetic SHACL fixture misses targets: {missing}")

    local_umls_nodes = {
        node
        for node in set(graph.subjects()) | set(graph.objects())
        if isinstance(node, URIRef)
        and str(node).startswith("http://makaao.inria.fr/kg/CUI_C")
    }
    if not local_umls_nodes:
        raise AssertionError("Synthetic SHACL fixture has no local UMLS node")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    graph = build_fixture()
    graph.serialize(args.output, format="turtle")
    print(f"Wrote {len(graph)} triples to {args.output}")


if __name__ == "__main__":
    main()
