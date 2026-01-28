from types import ModuleType
from pathlib import Path

from rdflib import Graph, Namespace, URIRef, RDF, RDFS
from rdflib.namespace import SKOS


def exec_script_with_overrides(script_path: Path, overrides: dict) -> ModuleType:
    code = script_path.read_text(encoding="utf-8", errors="ignore")
    m = ModuleType("kg_simplifier_exec")
    for k, v in overrides.items():
        setattr(m, k, v)
    exec(compile(code, str(script_path), "exec"), m.__dict__)
    return m


def test_04_simplifier_writes_output_graph(tmp_path):
    # Build a tiny input RDF/XML graph
    MAK = Namespace("http://makaao.inria.fr/kg/")
    UNI = Namespace("http://purl.uniprot.org/core/")
    SIO = Namespace("http://semanticscience.org/resource/")

    g = Graph()
    target_class = MAK["Target"]
    protein_class = UNI["Protein"]
    disease_class = MAK["AutoimmuneDisease"]

    target_inst = URIRef(str(target_class) + "_instance")
    disease_inst = URIRef(str(disease_class) + "_instance")

    g.add((target_class, RDFS.subClassOf, protein_class))
    g.add((target_class, RDFS.label, g.literal("Target")))
    g.add((target_inst, RDF.type, target_class))
    g.add((target_inst, RDFS.label, g.literal("Target inst")))

    g.add((disease_class, RDFS.label, g.literal("Autoimmune disease")))
    g.add((disease_inst, RDF.type, disease_class))
    g.add((disease_inst, SKOS.prefLabel, g.literal("Autoimmune disease inst")))

    # add one relation predicate the script is likely to keep
    g.add((target_inst, SIO["SIO_001403"], disease_inst))

    in_path = tmp_path / "in.rdf"
    out_path = tmp_path / "out.rdf"
    g.serialize(destination=str(in_path), format="xml")

    # Execute the script with overridden paths
    exec_script_with_overrides(
        Path("scripts/04_make_lite_graph_from_makaao-kg.py"),
        overrides={"IN_PATH": str(in_path), "OUT_PATH": str(out_path)},
    )

    assert out_path.exists()

    gout = Graph()
    gout.parse(str(out_path), format="xml")
    assert len(gout) > 0
