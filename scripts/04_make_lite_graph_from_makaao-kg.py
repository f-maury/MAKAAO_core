#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a class-only MAKAAO lite graph.

The current MAKAAO build separates:
  * the canonical KG, which contains domain instances, reification, provenance,
    and dataset metadata; and
  * the strict merged TBox, which contains the selected external class schema.

This script combines only the useful parts of those two files:
  * named classes;
  * rdfs:label and skos:prefLabel;
  * named rdfs:subClassOf links; and
  * selected domain relationships projected from instance-to-instance links to
    class-to-class links.

It intentionally excludes reification, PROV, dataset metadata, imports,
anonymous OWL expressions, and all individuals. Every resource occurring as a
relationship endpoint in the output is explicitly declared owl:Class.
"""

from __future__ import annotations

import argparse
import ast
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

from rdflib import BNode, Graph, Literal, Namespace, OWL, RDF, RDFS, URIRef
from rdflib.namespace import SKOS


SCRIPT_VERSION = "2.0.10"

# ----------------------------- Namespaces -----------------------------
MAKAAO = Namespace("http://makaao.inria.fr/kg/")
MAKAAO_LOINC = Namespace("http://makaao.inria.fr/loinc/")
LOINC_PROPERTY = Namespace("http://loinc.org/property/")
OBO = Namespace("http://purl.obolibrary.org/obo/")
UNIPROT_CORE = Namespace("http://purl.uniprot.org/core/")
SIO = Namespace("http://semanticscience.org/resource/")
BIOLINK = Namespace("https://w3id.org/biolink/vocab/")
BAO = Namespace("http://www.bioassayontology.org/bao#")
PROV = Namespace("http://www.w3.org/ns/prov#")
ORDO = Namespace("http://www.orpha.net/ORDO/")
UMLS = Namespace("http://uts.nlm.nih.gov/uts/umls/concept/")

LOINC_COMPONENT = LOINC_PROPERTY["COMPONENT"]

EXCLUDED_CLASSES = {MAKAAO.Document, MAKAAO.Relation} # exclude reified relations
GENERIC_INSTANCE_TYPES = {
    RDF.Statement,
    OWL.NamedIndividual,
    SKOS.Concept,
    MAKAAO.Document,
    MAKAAO.Relation,
    MAKAAO.Target,
    MAKAAO.AutoimmuneDisease,
    MAKAAO.CUI,
}

# Domain predicates retained in the class-only graph. Molecular-target links
# use the two official BAO predicates only.
RELATION_PREDICATES = {
    SIO["SIO_001403"],       # is associated with
    SIO["SIO_001279"],       # has phenotype
    SIO["SIO_001280"],       # is phenotype of
    BIOLINK.biomarker_for,
    BIOLINK.has_biomarker,
    LOINC_COMPONENT,
    SKOS.closeMatch,
    BAO["BAO_0000598"],      # is target for
    BAO["BAO_0000211"],      # has target
}

# Inverses retained or materialized in the lite graph to preserve bidirectional
# access. Duplicate triples are automatically ignored when the canonical KG
# already asserts both directions.
INVERSE_PREDICATES = {
    BIOLINK.biomarker_for: BIOLINK.has_biomarker,
    BIOLINK.has_biomarker: BIOLINK.biomarker_for,
    SIO["SIO_001279"]: SIO["SIO_001280"],
    SIO["SIO_001280"]: SIO["SIO_001279"],
    BAO["BAO_0000598"]: BAO["BAO_0000211"],
    BAO["BAO_0000211"]: BAO["BAO_0000598"],
}

ALLOWED_OUTPUT_PREDICATES = {
    RDF.type,
    RDFS.label,
    RDFS.subClassOf,
    SKOS.prefLabel,
    *RELATION_PREDICATES,
    *INVERSE_PREDICATES.values(),
}

REIFICATION_PREDICATES = {RDF.subject, RDF.predicate, RDF.object}


# ----------------------------- Paths -----------------------------
def find_project_dir(start: Path) -> Path:
    """Find the project root from project/, scripts/, or scripts/04/."""
    for candidate in (start, *start.parents):
        if (candidate / "data").is_dir():
            return candidate
    raise RuntimeError(
        f"Could not locate the MAKAAO project root above {start}; "
        "expected a data/ directory."
    )

# read some parameters from 03 graph building script
def read_builder_kg_version(project_dir: Path, script_dir: Path) -> str | None:
    """Read KG_VERSION from the main builder without importing/executing it."""
    candidates = (
        script_dir / "03_build_kg_from_tables.py",
        project_dir / "scripts" / "03_build_kg_from_tables.py",
        project_dir / "03_build_kg_from_tables.py",
    )
    for path in candidates:
        if not path.is_file():
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (OSError, SyntaxError, UnicodeError):
            continue
        for node in tree.body:
            if not isinstance(node, (ast.Assign, ast.AnnAssign)):
                continue
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if not any(isinstance(t, ast.Name) and t.id == "KG_VERSION" for t in targets):
                continue
            value = node.value
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                return value.value
    return None

# detect original KG in the folder
def discover_canonical_kg(kg_dir: Path) -> Path:
    """Select the only canonical makaao_kg_<version>.rdf file in kg/."""
    candidates = []
    for path in kg_dir.glob("makaao_kg_*.rdf"):
        name = path.name
        if any(
            marker in name
            for marker in (
                "_lite",
                "_no-reification",
                "_import-closure",
                "_curated-",
            )
        ):
            continue
        candidates.append(path)
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"No canonical makaao_kg_*.rdf found in {kg_dir}")
    names = "\n  ".join(sorted(p.name for p in candidates))
    raise RuntimeError(
        "More than one canonical KG candidate was found; specify --input explicitly:\n  "
        + names
    )

# set some paths for KG-related files
def default_paths() -> tuple[Path, Path, Path | None]:
    script_dir = Path(__file__).resolve().parent
    project_dir = find_project_dir(script_dir)
    kg_dir = project_dir / "kg"
    version = read_builder_kg_version(project_dir, script_dir)

    if version:
        input_path = kg_dir / f"makaao_kg_{version}.rdf"
        if not input_path.is_file():
            input_path = discover_canonical_kg(kg_dir)
    else:
        input_path = discover_canonical_kg(kg_dir)
        version = input_path.name[len("makaao_kg_") : -len(".rdf")]

    output_path = kg_dir / f"makaao_kg_{version}_lite.rdf"
    merged_tbox = (
        kg_dir
        / "reasoning"
        / f"makaao_kg_{version}_curated-tbox-merged.owl"
    )
    if not merged_tbox.is_file():
        fallback = kg_dir / f"makaao_kg_{version}_ontology.owl"
        merged_tbox = fallback if fallback.is_file() else None
    return input_path, output_path, merged_tbox


# ----------------------------- RDF helpers -----------------------------
# select URI part
def tail(term: URIRef) -> str:
    return str(term).rsplit("/", 1)[-1]

# detect entities from PORV-O ontology
def is_prov_term(term) -> bool:
    return isinstance(term, URIRef) and str(term).startswith(str(PROV))

# detect URI of entities to remove
def is_instance_uri(term) -> bool:
    if not isinstance(term, URIRef):
        return False
    text = str(term)
    slug = tail(term)
    if text.endswith("_instance"):
        return True
    if text.startswith(str(MAKAAO)):
        lowered = slug.lower()
        if lowered.startswith(("document_", "relation_")):
            return True
        if slug.startswith("r") and slug[1:].isdigit():
            return True
    return False

# detect some classes
def is_named_class_candidate(term) -> bool:
    return (
        isinstance(term, URIRef)
        and not is_instance_uri(term)
        and not is_prov_term(term)
        and term not in EXCLUDED_CLASSES
    )

# more function to detect entities from different types
def is_hpo_class(term: URIRef) -> bool:
    return str(term).startswith(str(OBO)) and tail(term).startswith("HP_")


def is_chebi_class(term: URIRef) -> bool:
    return str(term).startswith(str(OBO)) and tail(term).startswith("CHEBI_")


def is_ordo_class(term: URIRef) -> bool:
    return str(term).startswith(str(ORDO)) and tail(term).startswith("Orphanet_")


def is_local_up_class(term: URIRef) -> bool:
    return str(term).startswith(str(MAKAAO)) and tail(term).startswith("UP_")


def is_cui_class(term: URIRef) -> bool:
    text = str(term)
    slug = tail(term)
    return (
        (text.startswith(str(MAKAAO)) and slug.startswith("CUI_C"))
        or "/concept/C" in text
    )

# alternative method to detect some types of entities
def fallback_class_from_instance(node: URIRef) -> URIRef | None:
    """Fallback for older/current MAKAAO naming when rdf:type is incomplete."""
    slug = tail(node)
    if not slug.endswith("_instance"):
        return None
    stem = slug[: -len("_instance")]

    if stem == "aab_18":
        return MAKAAO.Autoantibody
    if stem == "positivity_18":
        return MAKAAO.AutoantibodyPositivity
    if stem.startswith("orpha_") and stem[len("orpha_") :].isdigit():
        return ORDO[f"Orphanet_{stem[len('orpha_'):]}"]
    if stem.startswith("hpo_HP_"):
        return OBO[stem[len("hpo_") :]]
    if stem.startswith("CHEBI_"):
        return OBO[stem]
    if stem.startswith("loinc_"):
        return MAKAAO_LOINC[stem[len("loinc_") :]]
    return MAKAAO[stem]

# build list of different types in KG
def build_types_index(graph: Graph) -> dict[URIRef, set[URIRef]]:
    result: dict[URIRef, set[URIRef]] = defaultdict(set)
    for subject, obj in graph.subject_objects(RDF.type):
        if isinstance(subject, URIRef) and isinstance(obj, URIRef):
            result[subject].add(obj)
    return result

# replace some insances connected by a relations, by corresponding classes, in the relation head and tail
def choose_projection_class(
    node,
    types_by_node: dict[URIRef, set[URIRef]],
) -> URIRef | None:
    """Map a relationship endpoint from an individual to its domain class."""
    if not isinstance(node, URIRef):
        return None
    if not is_instance_uri(node):
        return node if is_named_class_candidate(node) else None

    candidates = {
        rdf_type
        for rdf_type in types_by_node.get(node, set())
        if is_named_class_candidate(rdf_type)
        and rdf_type not in GENERIC_INSTANCE_TYPES
        and rdf_type
        not in {
            OWL.Class,
            OWL.ObjectProperty,
            OWL.DatatypeProperty,
            OWL.AnnotationProperty,
        }
    }

    # Positivity individuals project only to their local positivity class.
    # HPO classes are mapping targets and must never replace that local class.
    if tail(node).startswith("positivity_"):
        local_candidates = sorted(
            (
                candidate
                for candidate in candidates
                if candidate == MAKAAO.AutoantibodyPositivity
                or (
                    str(candidate).startswith(str(MAKAAO))
                    and tail(candidate).startswith("positivity_")
                )
            ),
            key=str,
        )
        if len(local_candidates) == 1:
            return local_candidates[0]
        if len(local_candidates) > 1:
            fallback = fallback_class_from_instance(node)
            if fallback in local_candidates:
                return fallback
            raise RuntimeError(
                "Ambiguous local positivity classes for "
                f"{node}: {', '.join(map(str, local_candidates))}"
            )
        hpo_candidates = sorted((c for c in candidates if is_hpo_class(c)), key=str)
        if hpo_candidates:
            raise RuntimeError(
                f"Positivity individual {node} is typed only with external HPO "
                "class(es); a local positivity class is required"
            )

    if len(candidates) == 1:
        return next(iter(candidates))

    fallback = fallback_class_from_instance(node)
    if fallback in candidates:
        return fallback

    if not candidates:
        return fallback

    raise RuntimeError(
        "Ambiguous instance-to-class projection for "
        f"{node}: {', '.join(sorted(map(str, candidates)))}"
    )

# apply lbels from instances to corresponding classes
def copy_labels(
    source_graphs: Iterable[Graph],
    source_iri: URIRef,
    output: Graph,
    *,
    target_iri: URIRef | None = None,
) -> None:
    destination = target_iri or source_iri
    for source in source_graphs:
        for predicate in (RDFS.label, SKOS.prefLabel):
            for value in source.objects(source_iri, predicate):
                if isinstance(value, Literal):
                    output.add((destination, predicate, value))

# check if a class has label
def has_label(graph: Graph, class_iri: URIRef) -> bool:
    return any(graph.objects(class_iri, RDFS.label)) or any(
        graph.objects(class_iri, SKOS.prefLabel)
    )

# return all classes from a KG
def all_class_nodes(graph: Graph) -> set[URIRef]:
    nodes: set[URIRef] = set()
    for subject, predicate, obj in graph:
        if predicate == RDF.type and obj == OWL.Class and isinstance(subject, URIRef):
            nodes.add(subject)
        elif predicate == RDFS.subClassOf:
            if isinstance(subject, URIRef):
                nodes.add(subject)
            if isinstance(obj, URIRef):
                nodes.add(obj)
        elif predicate in RELATION_PREDICATES or predicate in INVERSE_PREDICATES.values():
            if isinstance(subject, URIRef):
                nodes.add(subject)
            if isinstance(obj, URIRef):
                nodes.add(obj)
    return nodes


# ----------------------------- Projection -----------------------------
# collect classes in KG
def collect_schema_classes(schema: Graph) -> set[URIRef]:
    classes = {
        subject
        for subject in schema.subjects(RDF.type, OWL.Class)
        if is_named_class_candidate(subject)
    }
    for subject, obj in schema.subject_objects(RDFS.subClassOf):
        if is_named_class_candidate(subject) and is_named_class_candidate(obj):
            classes.add(subject)
            classes.add(obj)
    return classes

# big function to build lite KG from original KG, check result, and write reports
def build_lite_graph(data: Graph, schema: Graph) -> tuple[Graph, dict]:
    output = Graph()
    for prefix, namespace in (
        ("mak", MAKAAO),
        ("mloinc", MAKAAO_LOINC),
        ("loinc_property", LOINC_PROPERTY),
        ("rdfs", RDFS),
        ("owl", OWL),
        ("skos", SKOS),
        ("hp", OBO),
        ("ordo", ORDO),
        ("sio", SIO),
        ("biolink", BIOLINK),
        ("bao", BAO),
        ("up", UNIPROT_CORE),
    ):
        output.bind(prefix, namespace)

    types_by_node = build_types_index(data)
    schema_classes = collect_schema_classes(schema)

    # Named class schema only: declarations, labels and named subclass links.
    for class_iri in sorted(schema_classes, key=str):
        output.add((class_iri, RDF.type, OWL.Class))
        copy_labels((schema, data), class_iri, output)

    for subject, obj in schema.subject_objects(RDFS.subClassOf):
        if (
            subject in schema_classes
            and obj in schema_classes
            and subject != obj
            and not is_prov_term(subject)
            and not is_prov_term(obj)
        ):
            output.add((subject, RDFS.subClassOf, obj))

    relation_counts: Counter[str] = Counter()
    endpoint_projection: dict[URIRef, URIRef] = {}
    skipped_relationships = 0

    for predicate in sorted(RELATION_PREDICATES, key=str):
        for subject, obj in data.subject_objects(predicate):
            subject_class = choose_projection_class(subject, types_by_node)
            object_class = choose_projection_class(obj, types_by_node)
            if subject_class is None or object_class is None:
                skipped_relationships += 1
                continue
            if subject_class in EXCLUDED_CLASSES or object_class in EXCLUDED_CLASSES:
                skipped_relationships += 1
                continue
            if is_prov_term(subject_class) or is_prov_term(object_class):
                skipped_relationships += 1
                continue

            output.add((subject_class, RDF.type, OWL.Class))
            output.add((object_class, RDF.type, OWL.Class))
            output.add((subject_class, predicate, object_class))
            relation_counts[str(predicate)] += 1

            if isinstance(subject, URIRef) and is_instance_uri(subject):
                endpoint_projection[subject] = subject_class
            if isinstance(obj, URIRef) and is_instance_uri(obj):
                endpoint_projection[obj] = object_class

            inverse = INVERSE_PREDICATES.get(predicate)
            if inverse is not None:
                output.add((object_class, inverse, subject_class))

    # Reuse individual labels only when the selected class has no class label.
    for instance, class_iri in endpoint_projection.items():
        if has_label(output, class_iri):
            continue
        copy_labels((data,), instance, output, target_iri=class_iri)

    role_by_class: dict[URIRef, set[URIRef]] = defaultdict(set)
    for instance, class_iri in endpoint_projection.items():
        for rdf_type in types_by_node.get(instance, set()):
            if rdf_type in {
                MAKAAO.Target,
                MAKAAO.AutoimmuneDisease,
                MAKAAO.AutoantibodyPositivity,
            }:
                role_by_class[class_iri].add(rdf_type)

    # A UMLS concept aligned to an ORDO class is disease evidence even when
    # the enrichment-only CUI has no individual participating in a relation.
    for left, right in data.subject_objects(SKOS.closeMatch):
        if not isinstance(left, URIRef) or not isinstance(right, URIRef):
            continue
        if is_cui_class(left) and is_ordo_class(right):
            role_by_class[left].add(MAKAAO.AutoimmuneDisease)
        elif is_ordo_class(left) and is_cui_class(right):
            role_by_class[right].add(MAKAAO.AutoimmuneDisease)

    add_application_role_hierarchy(output, role_by_class)
    purge_excluded_resources(output)
    ensure_every_endpoint_is_a_class(output)

    relationship_predicates = RELATION_PREDICATES | set(INVERSE_PREDICATES.values())
    output_relation_counts = Counter(
        str(predicate)
        for _, predicate, _ in output
        if predicate in relationship_predicates
    )

    report = {
        "input_triples": len(data),
        "schema_triples": len(schema),
        "output_triples": len(output),
        "classes": len(set(output.subjects(RDF.type, OWL.Class))),
        "projected_instances": len(endpoint_projection),
        "skipped_relationships": skipped_relationships,
        "relationship_source_counts": dict(sorted(relation_counts.items())),
        "relationship_output_counts": dict(sorted(output_relation_counts.items())),
    }
    return output, report

# add makaao classes hierarchy
def add_application_role_hierarchy(
    output: Graph,
    role_by_class: dict[URIRef, set[URIRef]],
) -> None:
    """Add evidence-backed MAKAAO roles without rewriting source hierarchies."""
    classes = all_class_nodes(output)

    for class_iri in sorted(classes, key=str):
        roles = role_by_class.get(class_iri, set())
        supported_roles: set[URIRef] = set()

        if (
            MAKAAO.Target in roles
            and (
                is_chebi_class(class_iri)
                or is_local_up_class(class_iri)
                or is_cui_class(class_iri)
            )
        ):
            supported_roles.add(MAKAAO.Target)

        if (
            MAKAAO.AutoimmuneDisease in roles
            and (is_ordo_class(class_iri) or is_cui_class(class_iri))
        ):
            supported_roles.add(MAKAAO.AutoimmuneDisease)

        for parent in supported_roles:
            output.add((class_iri, RDFS.subClassOf, parent))
            output.add((parent, RDF.type, OWL.Class))

        # Keep an otherwise unclassified CUI in the neutral local CUI branch;
        # absence of disease evidence is not evidence that it is a target.
        if is_cui_class(class_iri) and not supported_roles:
            output.add((class_iri, RDFS.subClassOf, MAKAAO.CUI))
            output.add((MAKAAO.CUI, RDF.type, OWL.Class))

        output.add((class_iri, RDF.type, OWL.Class))

# remove entity types that are not allowed in lite KG
def purge_excluded_resources(output: Graph) -> None:
    for banned in EXCLUDED_CLASSES:
        output.remove((banned, None, None))
        output.remove((None, None, banned))

    for triple in list(output):
        subject, predicate, obj = triple
        if predicate not in ALLOWED_OUTPUT_PREDICATES:
            output.remove(triple)
            continue
        if isinstance(subject, BNode) or isinstance(obj, BNode):
            output.remove(triple)
            continue
        if is_prov_term(subject) or is_prov_term(predicate) or is_prov_term(obj):
            output.remove(triple)
            continue
        if is_instance_uri(subject) or is_instance_uri(obj):
            output.remove(triple)
            continue
        if predicate in REIFICATION_PREDICATES:
            output.remove(triple)
            continue
        if predicate == RDF.type and obj != OWL.Class:
            output.remove(triple)
            continue
        if predicate in {RDFS.label, SKOS.prefLabel} and not isinstance(obj, Literal):
            output.remove(triple)
            continue
        if predicate == RDFS.subClassOf and not isinstance(obj, URIRef):
            output.remove(triple)

# make sure we have only classes and no instances
def ensure_every_endpoint_is_a_class(output: Graph) -> None:
    class_nodes = all_class_nodes(output)
    for class_iri in class_nodes:
        if class_iri not in EXCLUDED_CLASSES and not is_prov_term(class_iri):
            output.add((class_iri, RDF.type, OWL.Class))


# ----------------------------- Validation/output -----------------------------
# functions to perform check on the resulting lite KG
def validate_lite_graph(graph: Graph) -> None:
    errors: list[str] = []

    for subject, predicate, obj in graph:
        if predicate not in ALLOWED_OUTPUT_PREDICATES:
            errors.append(f"unexpected predicate: {predicate}")
        if isinstance(subject, BNode) or isinstance(obj, BNode):
            errors.append(f"blank node retained: {(subject, predicate, obj)}")
        if is_instance_uri(subject) or is_instance_uri(obj):
            errors.append(f"instance URI retained: {(subject, predicate, obj)}")
        if is_prov_term(subject) or is_prov_term(predicate) or is_prov_term(obj):
            errors.append(f"PROV term retained: {(subject, predicate, obj)}")
        if subject in EXCLUDED_CLASSES or obj in EXCLUDED_CLASSES:
            errors.append(f"excluded resource retained: {(subject, predicate, obj)}")
        if predicate == RDF.type and obj != OWL.Class:
            errors.append(f"non-class rdf:type retained: {(subject, predicate, obj)}")
        if predicate in REIFICATION_PREDICATES:
            errors.append(f"reification predicate retained: {predicate}")

    declared_classes = set(graph.subjects(RDF.type, OWL.Class))
    for subject, predicate, obj in graph:
        if predicate in RELATION_PREDICATES or predicate in INVERSE_PREDICATES.values():
            if subject not in declared_classes or obj not in declared_classes:
                errors.append(
                    f"relationship endpoint lacks owl:Class declaration: "
                    f"{(subject, predicate, obj)}"
                )
        if predicate == RDFS.subClassOf:
            if subject not in declared_classes or obj not in declared_classes:
                errors.append(
                    f"subclass endpoint lacks owl:Class declaration: "
                    f"{(subject, predicate, obj)}"
                )

    if errors:
        preview = "\n  ".join(errors[:25])
        suffix = "" if len(errors) <= 25 else f"\n  ... {len(errors) - 25} more"
        raise RuntimeError(
            f"Lite graph validation failed with {len(errors)} issue(s):\n  "
            + preview
            + suffix
        )

# export new lite KG
def serialize_atomically(graph: Graph, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent
    )
    os.close(fd)
    temporary_path = Path(temporary_name)
    try:
        graph.serialize(destination=str(temporary_path), format="xml")
        check = Graph()
        try:
            check.parse(temporary_path, format="xml")
            if len(check) != len(graph):
                raise RuntimeError(
                    f"RDF/XML round-trip count mismatch: {len(graph)} != {len(check)}"
                )
            validate_lite_graph(check)
        finally:
            check.close()
        os.replace(temporary_path, output_path)
        try:
            output_path.chmod(0o644)
        except OSError:
            # Some Windows-mounted WSL filesystems ignore POSIX mode changes.
            pass
    finally:
        temporary_path.unlink(missing_ok=True)

# function to perform transofrmation from original KG to lite KG as a command line call (do we really need this? maybe remove it)
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a class-only MAKAAO lite graph from the current full KG."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="canonical full KG; auto-discovered from the project when omitted",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=None,
        help=(
            "merged TBox used for labels and class hierarchy; derived from the "
            "input path or project when omitted"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="lite RDF/XML output; derived from the input path when omitted",
    )
    return parser.parse_args()

# obtain new file names from original file names, to use as output path
def derive_paths_from_explicit_input(input_path: Path) -> tuple[Path, Path]:
    """Derive output and schema paths without requiring a project checkout."""
    name = input_path.name
    match = re.fullmatch(r"makaao_kg_(.+)\.rdf", name)
    if match:
        version = match.group(1)
        output_path = input_path.with_name(f"makaao_kg_{version}_lite.rdf")
        merged_tbox = (
            input_path.parent
            / "reasoning"
            / f"makaao_kg_{version}_curated-tbox-merged.owl"
        )
        ontology_fallback = input_path.with_name(
            f"makaao_kg_{version}_ontology.owl"
        )
        if merged_tbox.is_file():
            schema_path = merged_tbox
        elif ontology_fallback.is_file():
            schema_path = ontology_fallback
        else:
            schema_path = input_path
        return output_path, schema_path

    return (
        input_path.with_name(f"{input_path.stem}_lite{input_path.suffix}"),
        input_path,
    )

# main function where we call defined functions
def main() -> None:
    args = parse_args()

    if args.input is None:
        default_input, default_output, default_schema = default_paths()
        input_choice = default_input
        output_choice = args.output or default_output
        schema_choice = args.schema or default_schema or default_input
    else:
        input_choice = args.input
        derived_output, derived_schema = derive_paths_from_explicit_input(
            input_choice.expanduser().resolve()
        )
        output_choice = args.output or derived_output
        schema_choice = args.schema or derived_schema

    input_path = input_choice.expanduser().resolve()
    output_path = output_choice.expanduser().resolve()
    schema_path = schema_choice.expanduser().resolve()

    if not input_path.is_file():
        raise FileNotFoundError(f"Canonical KG not found: {input_path}")
    if not schema_path.is_file():
        raise FileNotFoundError(f"Schema graph not found: {schema_path}")
    if input_path == output_path:
        raise ValueError("Input and output paths must be different")
    if schema_path == output_path:
        raise ValueError("Schema and output paths must be different")

    data = Graph()
    schema = Graph()
    try:
        data.parse(input_path)
        if schema_path == input_path:
            schema = data
        else:
            schema.parse(schema_path)

        lite, report = build_lite_graph(data, schema)
        validate_lite_graph(lite)
        serialize_atomically(lite, output_path)
    finally:
        if schema is not data:
            schema.close()
        data.close()

    print(f"MAKAAO lite build complete: script={SCRIPT_VERSION}")
    print(f"Input KG:       {input_path}")
    print(f"Schema graph:   {schema_path}")
    print(f"Output:         {output_path}")
    print(f"Input triples:  {report['input_triples']}")
    print(f"Schema triples: {report['schema_triples']}")
    print(f"Output triples: {report['output_triples']}")
    print(f"Classes:        {report['classes']}")
    print(f"Projected instances: {report['projected_instances']}")
    print(f"Skipped relationships: {report['skipped_relationships']}")
    print("Class-level relationship source counts:")
    for predicate, count in report["relationship_source_counts"].items():
        print(f"  {predicate}: {count}")
    print("Final class-level relationship counts:")
    for predicate, count in report["relationship_output_counts"].items():
        print(f"  {predicate}: {count}")


if __name__ == "__main__":
    main()