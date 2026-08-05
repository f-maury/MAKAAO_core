#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the strict OWL 2 DL reasoning release for MAKAAO.

This is an internal helper called by ``03_build_kg_from_tables.py``. It builds
the only external ontology module set directly from the pinned ontology files
in ``data/external_ontologies``. No intermediate or previously extracted
ontology modules are used as source material.

The reasoning modules retain:

* classes and properties referenced by the generated MAKAAO graph;
* named superclass and superproperty closure within each source vocabulary;
* named inverse/equivalent-property relations;
* named domains and ranges;
* safe literal annotations.

Anonymous restrictions, RDF lists, n-ary axioms, and unresolved property-category
overlap are intentionally excluded from this strict reasoning representation.
"""

from __future__ import annotations

import gc
import hashlib
import json
import os
import signal
import subprocess
import tempfile
import threading
from collections import deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence
from xml.etree import ElementTree as ET

from rdflib import BNode, Graph, Literal, Namespace, RDF, RDFS, OWL, URIRef
from rdflib.namespace import DCTERMS, SKOS, XSD


SCRIPT_VERSION = "1.2.28"
SCRIPT_ITERATION = "2026-08-05-add-ordo-umls-exact-label-matches"

PROV = Namespace("http://www.w3.org/ns/prov#")
BIOLINK = Namespace("https://w3id.org/biolink/vocab/")
BAO = Namespace("http://www.bioassayontology.org/bao#")
BAO_HAS_TARGET = BAO["BAO_0000211"]
BAO_IS_TARGET_FOR = BAO["BAO_0000598"]
MAKAAO = Namespace("http://makaao.inria.fr/kg/")
SKOS_ONTOLOGY_IRI = URIRef("http://www.w3.org/2004/02/skos/core")
DCAT = Namespace("http://www.w3.org/ns/dcat#")
VCARD = Namespace("http://www.w3.org/2006/vcard/ns#")
FOAF = Namespace("http://xmlns.com/foaf/0.1/")
VOID = Namespace("http://rdfs.org/ns/void#")
SCHEMA = Namespace("http://schema.org/")
ODRL = Namespace("http://www.w3.org/ns/odrl/2/")

MODULE_BASE = "http://makaao.inria.fr/imports/"


@dataclass(frozen=True)
class VocabularySpec:
    name: str
    prefixes: tuple[str, ...]
    module_iri: URIRef
    filename: str
    source_key: str


VOCABULARIES = (
    VocabularySpec(
        "biolink",
        ("https://w3id.org/biolink/vocab/",),
        URIRef(MODULE_BASE + "biolink"),
        "biolink-makaao-reasoning-module.owl",
        "biolink",
    ),
    VocabularySpec(
        "hpo",
        ("http://purl.obolibrary.org/obo/HP_",),
        URIRef(MODULE_BASE + "hpo"),
        "hpo-makaao-reasoning-module.owl",
        "hpo",
    ),
    VocabularySpec(
        "chebi",
        ("http://purl.obolibrary.org/obo/CHEBI_",),
        URIRef(MODULE_BASE + "chebi"),
        "chebi-makaao-reasoning-module.owl",
        "chebi",
    ),
    VocabularySpec(
        "ordo",
        ("http://www.orpha.net/ORDO/",),
        URIRef(MODULE_BASE + "ordo"),
        "ordo-makaao-reasoning-module.owl",
        "ordo",
    ),
    VocabularySpec(
        "sio",
        ("http://semanticscience.org/resource/",),
        URIRef(MODULE_BASE + "sio"),
        "sio-makaao-reasoning-module.owl",
        "sio",
    ),
    VocabularySpec(
        "prov",
        ("http://www.w3.org/ns/prov#",),
        URIRef(MODULE_BASE + "prov"),
        "prov-makaao-reasoning-module.owl",
        "prov",
    ),
    VocabularySpec(
        "bao",
        ("http://www.bioassayontology.org/bao#",),
        URIRef(MODULE_BASE + "bao"),
        "bao-makaao-reasoning-module.owl",
        "bao",
    ),
    VocabularySpec(
        "uniprot",
        ("http://purl.uniprot.org/core/",),
        URIRef(MODULE_BASE + "uniprot-core"),
        "uniprot-core-makaao-reasoning-module.owl",
        "uniprot",
    ),
    VocabularySpec(
        "skos",
        ("http://www.w3.org/2004/02/skos/core#",),
        URIRef("http://www.w3.org/TR/skos-reference/skos-owl1-dl.rdf"),
        "skos-makaao-reasoning-module.owl",
        "skos",
    ),
)

SKOS_NOTATION = URIRef(str(SKOS) + "notation")
SKOS_PREF_LABEL = URIRef(str(SKOS) + "prefLabel")

LOGICAL_PREDICATES = {
    RDF.type,
    RDFS.subClassOf,
    RDFS.subPropertyOf,
    RDFS.domain,
    RDFS.range,
    OWL.equivalentClass,
    OWL.equivalentProperty,
    OWL.disjointWith,
    OWL.propertyDisjointWith,
    OWL.inverseOf,
    OWL.intersectionOf,
    OWL.unionOf,
    OWL.complementOf,
    OWL.oneOf,
    OWL.onProperty,
    OWL.someValuesFrom,
    OWL.allValuesFrom,
    OWL.hasValue,
    OWL.hasSelf,
    OWL.cardinality,
    OWL.minCardinality,
    OWL.maxCardinality,
    OWL.qualifiedCardinality,
    OWL.minQualifiedCardinality,
    OWL.maxQualifiedCardinality,
    OWL.onClass,
    OWL.onDataRange,
    OWL.propertyChainAxiom,
    OWL.members,
    OWL.distinctMembers,
    OWL.disjointUnionOf,
    OWL.hasKey,
    RDF.first,
    RDF.rest,
}

OBJECT_PROPERTY_CHARACTERISTICS = {
    OWL.SymmetricProperty,
    OWL.TransitiveProperty,
    OWL.FunctionalProperty,
    OWL.InverseFunctionalProperty,
    OWL.AsymmetricProperty,
    OWL.ReflexiveProperty,
    OWL.IrreflexiveProperty,
}

# These characteristics are legal only for object properties. FunctionalProperty
# is deliberately excluded because it is shared by object and datatype
# properties. Structural recognition is required for sources such as the pinned
# Biolink OWL export, where an object-hierarchy property may also carry an
# erroneous owl:DatatypeProperty declaration but is typed owl:SymmetricProperty.
OBJECT_ONLY_PROPERTY_CHARACTERISTICS = {
    OWL.SymmetricProperty,
    OWL.TransitiveProperty,
    OWL.InverseFunctionalProperty,
    OWL.AsymmetricProperty,
    OWL.ReflexiveProperty,
    OWL.IrreflexiveProperty,
}

# Named built-ins that may legally occur as datatype-property ranges without a
# local rdfs:Datatype declaration. XSD datatypes are recognized by namespace.
BUILTIN_DATATYPES = {
    RDFS.Literal,
    RDF.XMLLiteral,
    RDF.HTML,
    RDF.langString,
    RDF.PlainLiteral,
    OWL.real,
    OWL.rational,
}


def is_builtin_datatype_iri(node: Any) -> bool:
    """Return whether an IRI is a built-in RDF/OWL/XSD datatype."""
    return (
        isinstance(node, URIRef)
        and (node in BUILTIN_DATATYPES or str(node).startswith(str(XSD)))
    )

METADATA_NAMESPACES = (
    str(DCTERMS),
    str(DCAT),
    str(VCARD),
    str(FOAF),
    str(VOID),
    str(SCHEMA),
    str(ODRL),
)


@dataclass(frozen=True)
class ModuleResult:
    name: str
    path: str
    source_path: str
    source_sha256: str
    output_sha256: str
    used_terms: int
    used_predicates: int
    selected_classes: int
    selected_object_properties: int
    selected_datatype_properties: int
    selected_annotation_properties: int
    selected_datatypes: int
    triples: int
    direct_term_iris: tuple[str, ...]
    direct_predicate_iris: tuple[str, ...]
    selected_class_iris: tuple[str, ...]
    selected_object_property_iris: tuple[str, ...]
    selected_datatype_property_iris: tuple[str, ...]
    selected_annotation_property_iris: tuple[str, ...]
    selected_datatype_iris: tuple[str, ...]
    infrastructure_annotation_property_iris: tuple[str, ...]
    resolved_annotation_logical_overlap_iris: tuple[str, ...]
    retained_axiom_counts: dict[str, int]


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portable_source_path(path: str | Path) -> str:
    """Return a stable project-relative source label for manifests/modules."""
    resolved = Path(path).resolve()
    parts = resolved.parts
    if "external_ontologies" in parts:
        index = parts.index("external_ontologies")
        return str(Path("data", *parts[index:]))
    return resolved.name


def parse_graph(path: str | Path) -> Graph:
    path = str(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    errors: list[str] = []
    for rdf_format in (None, "xml", "turtle", "nt", "json-ld"):
        graph = Graph()
        try:
            if rdf_format is None:
                graph.parse(path)
            else:
                graph.parse(path, format=rdf_format)
            return graph
        except Exception as exc:  # pragma: no cover - error aggregation
            graph.close()
            errors.append(f"{rdf_format or 'auto'}: {exc}")
    raise ValueError(f"Could not parse RDF source {path}:\n  " + "\n  ".join(errors))


def iri_in_prefixes(node: Any, prefixes: Sequence[str]) -> bool:
    return isinstance(node, URIRef) and any(str(node).startswith(prefix) for prefix in prefixes)


def collect_signature(graph: Graph, prefixes: Sequence[str]) -> tuple[set[URIRef], set[URIRef]]:
    terms: set[URIRef] = set()
    predicates: set[URIRef] = set()
    for subject, predicate, obj in graph:
        for node in (subject, predicate, obj):
            if iri_in_prefixes(node, prefixes):
                terms.add(node)
        if iri_in_prefixes(predicate, prefixes):
            predicates.add(predicate)
    return terms, predicates


def source_entity_sets(source: Graph) -> dict[str, set[URIRef]]:
    """Collect declared and structurally implied named entity categories.

    RDF/OWL exports occasionally omit declarations or assign conflicting
    property categories. Category evidence is propagated across named
    subproperty and equivalent-property structures so actual KG usage can
    resolve a strict OWL 2 DL projection without silently discarding a source
    contradiction. Object-only OWL constructs are tracked separately because
    they can never be represented as datatype or annotation properties.

    Named range axioms also provide category evidence: a class range indicates
    an object property, while a datatype range indicates a datatype property.
    When an explicitly single-category property uses an undeclared range IRI,
    the range is inferred as a class or custom datatype accordingly.
    """
    explicit_object_properties = {
        node
        for node in source.subjects(RDF.type, OWL.ObjectProperty)
        if isinstance(node, URIRef)
    }
    explicit_datatype_properties = {
        node
        for node in source.subjects(RDF.type, OWL.DatatypeProperty)
        if isinstance(node, URIRef)
    }
    explicit_annotation_properties = {
        node
        for node in source.subjects(RDF.type, OWL.AnnotationProperty)
        if isinstance(node, URIRef)
    }
    named_individuals = {
        node
        for node in source.subjects(RDF.type, OWL.NamedIndividual)
        if isinstance(node, URIRef)
    }

    datatypes = {
        node
        for node in source.subjects(RDF.type, RDFS.Datatype)
        if isinstance(node, URIRef)
    }
    classes = {
        node
        for class_type in (OWL.Class, RDFS.Class)
        for node in source.subjects(RDF.type, class_type)
        if isinstance(node, URIRef)
    }
    for predicate in (RDFS.subClassOf, OWL.equivalentClass, OWL.disjointWith):
        for left, right in source.subject_objects(predicate):
            if isinstance(left, URIRef):
                classes.add(left)
            if isinstance(right, URIRef):
                classes.add(right)

    object_only_properties: set[URIRef] = set()
    for characteristic in OBJECT_ONLY_PROPERTY_CHARACTERISTICS:
        object_only_properties.update(
            node
            for node in source.subjects(RDF.type, characteristic)
            if isinstance(node, URIRef)
        )
    for left, right in source.subject_objects(OWL.inverseOf):
        if isinstance(left, URIRef):
            object_only_properties.add(left)
        if isinstance(right, URIRef):
            object_only_properties.add(right)
    object_only_properties.update(
        node
        for node in source.subjects(OWL.propertyChainAxiom, None)
        if isinstance(node, URIRef)
    )

    # Domains are class expressions for both object and datatype properties.
    # A datatype-valued domain is retained as contradictory evidence and is
    # rejected later if the property enters the selected module.
    for target in source.objects(None, RDFS.domain):
        if (
            isinstance(target, URIRef)
            and not is_builtin_datatype_iri(target)
            and target not in datatypes
        ):
            classes.add(target)

    range_object_properties: set[URIRef] = set()
    range_datatype_properties: set[URIRef] = set()
    range_pairs = [
        (prop, target)
        for prop, target in source.subject_objects(RDFS.range)
        if isinstance(prop, URIRef) and isinstance(target, URIRef)
    ]

    # First use declarations already present on the range target.
    for prop, target in range_pairs:
        if is_builtin_datatype_iri(target) or target in datatypes:
            range_datatype_properties.add(prop)
        elif target in classes:
            range_object_properties.add(prop)

    # For a property with exactly one explicit logical category, an otherwise
    # undeclared named range can be classified without guessing.
    for prop, target in range_pairs:
        if is_builtin_datatype_iri(target) or target in datatypes or target in classes:
            continue
        object_evidence = prop in explicit_object_properties or prop in object_only_properties
        datatype_evidence = prop in explicit_datatype_properties
        if object_evidence and not datatype_evidence:
            classes.add(target)
            range_object_properties.add(prop)
        elif datatype_evidence and not object_evidence:
            datatypes.add(target)
            range_datatype_properties.add(prop)

    logical_properties: set[URIRef] = (
        explicit_object_properties
        | explicit_datatype_properties
        | object_only_properties
        | range_object_properties
        | range_datatype_properties
    )
    logical_properties.update(
        node
        for node in source.subjects(RDF.type, OWL.FunctionalProperty)
        if isinstance(node, URIRef)
    )
    for predicate in (OWL.equivalentProperty, OWL.propertyDisjointWith):
        for left, right in source.subject_objects(predicate):
            if isinstance(left, URIRef):
                logical_properties.add(left)
            if isinstance(right, URIRef):
                logical_properties.add(right)

    object_properties = (
        set(explicit_object_properties)
        | object_only_properties
        | range_object_properties
    )
    datatype_properties = (
        set(explicit_datatype_properties)
        | range_datatype_properties
    )
    annotation_properties = set(explicit_annotation_properties)

    # In OWL 2 DL both endpoints of a named subproperty/equivalence axiom must
    # belong to the same property category. Propagate already known category
    # evidence in both directions. Conflicting evidence is intentionally kept
    # in multiple sets; direct KG usage resolves it later or causes a focused
    # failure.
    changed = True
    while changed:
        before = (
            len(object_properties),
            len(datatype_properties),
            len(annotation_properties),
            len(object_only_properties),
            len(logical_properties),
        )
        for child, parent in source.subject_objects(RDFS.subPropertyOf):
            if not isinstance(child, URIRef) or not isinstance(parent, URIRef):
                continue
            pair = {child, parent}
            if pair & object_properties:
                object_properties.update(pair)
            if pair & datatype_properties:
                datatype_properties.update(pair)
            if pair & annotation_properties:
                # Annotation-property hierarchy evidence must not contaminate a
                # logical object/datatype hierarchy. Some published vocabularies
                # (notably PROV-O) contain annotation/logical overlap for the same
                # IRI or neighboring terms. Retain explicit overlap as evidence,
                # but do not propagate annotation status across an edge that
                # already has logical category evidence.
                if not pair & (
                    object_properties
                    | datatype_properties
                    | object_only_properties
                    | logical_properties
                ):
                    annotation_properties.update(pair)
            if pair & object_only_properties:
                object_only_properties.update(pair)
                object_properties.update(pair)
            if pair & logical_properties:
                logical_properties.update(pair)
        for left, right in source.subject_objects(OWL.equivalentProperty):
            if not isinstance(left, URIRef) or not isinstance(right, URIRef):
                continue
            pair = {left, right}
            if pair & object_properties:
                object_properties.update(pair)
            if pair & datatype_properties:
                datatype_properties.update(pair)
            if pair & annotation_properties:
                # Annotation-property hierarchy evidence must not contaminate a
                # logical object/datatype hierarchy. Some published vocabularies
                # (notably PROV-O) contain annotation/logical overlap for the same
                # IRI or neighboring terms. Retain explicit overlap as evidence,
                # but do not propagate annotation status across an edge that
                # already has logical category evidence.
                if not pair & (
                    object_properties
                    | datatype_properties
                    | object_only_properties
                    | logical_properties
                ):
                    annotation_properties.update(pair)
            if pair & object_only_properties:
                object_only_properties.update(pair)
                object_properties.update(pair)
            logical_properties.update(pair)
        for left, right in source.subject_objects(OWL.inverseOf):
            if isinstance(left, URIRef) and isinstance(right, URIRef):
                pair = {left, right}
                object_only_properties.update(pair)
                object_properties.update(pair)
                logical_properties.update(pair)
        after = (
            len(object_properties),
            len(datatype_properties),
            len(annotation_properties),
            len(object_only_properties),
            len(logical_properties),
        )
        changed = after != before

    return {
        "classes": classes,
        "object_properties": object_properties,
        "object_only_properties": object_only_properties,
        "datatype_properties": datatype_properties,
        "annotation_properties": annotation_properties,
        "logical_properties": logical_properties,
        "datatypes": datatypes,
        "named_individuals": named_individuals,
    }

def predicate_usage_kind(graph: Graph, predicate: URIRef) -> str | None:
    objects = list(graph.objects(None, predicate))
    if not objects:
        return None
    literal_count = sum(isinstance(obj, Literal) for obj in objects)
    resource_count = len(objects) - literal_count
    if literal_count and resource_count:
        return "mixed"
    if literal_count:
        return "literal"
    return "resource"


def assert_source_property_category_compatible(
    *,
    spec: VocabularySpec,
    predicate: URIRef,
    desired_category: str,
    source_object_properties: set[URIRef],
    source_object_only_properties: set[URIRef],
    source_datatype_properties: set[URIRef],
    source_annotation_properties: set[URIRef],
    source_logical_properties: set[URIRef],
) -> None:
    """Prevent hierarchy closure from silently changing source semantics."""
    if (
        spec.name == "skos"
        and predicate in {SKOS_NOTATION, SKOS_PREF_LABEL}
        and desired_category == "annotation"
    ):
        return

    is_object = predicate in source_object_properties
    is_datatype = predicate in source_datatype_properties
    is_annotation = predicate in source_annotation_properties
    is_logical = predicate in source_logical_properties

    # OWL 2 DL does not permit one IRI to be both an annotation property and a
    # logical property. Some upstream RDF serializations nevertheless contain
    # that overlap. For the strict projection, corresponding logical evidence
    # takes precedence over the annotation declaration; the conflicting
    # annotation declaration is omitted and recorded in the module audit.
    if desired_category == "object":
        if is_datatype and not is_object:
            raise RuntimeError(
                f"{spec.name}: property closure would reclassify {predicate} as an "
                "object property contrary to its pinned source declaration"
            )
        if is_annotation and not (is_object or predicate in source_object_only_properties):
            raise RuntimeError(
                f"{spec.name}: annotation-only property {predicate} cannot be "
                "reclassified as an object property"
            )
    elif desired_category == "datatype":
        if predicate in source_object_only_properties:
            raise RuntimeError(
                f"{spec.name}: property {predicate} has OWL structure legal only for "
                "object properties and cannot be resolved as a datatype property"
            )
        if is_object and not is_datatype:
            raise RuntimeError(
                f"{spec.name}: property closure would reclassify {predicate} as a "
                "datatype property contrary to its pinned source declaration"
            )
        if is_annotation and not is_datatype:
            raise RuntimeError(
                f"{spec.name}: annotation-only property {predicate} cannot be "
                "reclassified as a datatype property"
            )
    elif desired_category == "annotation":
        if is_logical or is_object or is_datatype:
            raise RuntimeError(
                f"{spec.name}: property closure would reclassify {predicate} as an "
                "annotation property contrary to its pinned source declaration"
            )
    else:  # pragma: no cover - internal programming error
        raise ValueError(f"Unknown property category: {desired_category}")


def add_module_header(
    module: Graph,
    spec: VocabularySpec,
    source_path: str,
    source_hash: str,
    kg_version: str,
) -> None:
    module.add((spec.module_iri, RDF.type, OWL.Ontology))
    module.add(
        (
            spec.module_iri,
            OWL.versionInfo,
            Literal(f"MAKAAO {kg_version} strict reasoning module for {spec.name}"),
        )
    )
    module.add((DCTERMS.source, RDF.type, OWL.AnnotationProperty))
    module.add((DCTERMS.identifier, RDF.type, OWL.AnnotationProperty))
    module.add(
        (
            spec.module_iri,
            DCTERMS.source,
            Literal(portable_source_path(source_path)),
        )
    )
    module.add((spec.module_iri, DCTERMS.identifier, Literal(f"sha256:{source_hash}")))


def expand_class_closure(
    source: Graph,
    selected: set[URIRef],
    classes: set[URIRef],
    prefixes: Sequence[str],
) -> None:
    changed = True
    while changed:
        previous = len(selected)
        for child, parent in source.subject_objects(RDFS.subClassOf):
            if (
                child in selected
                and isinstance(parent, URIRef)
                and parent in classes
                and parent != OWL.Thing
                and iri_in_prefixes(parent, prefixes)
            ):
                selected.add(parent)
        for left, right in source.subject_objects(OWL.equivalentClass):
            if not isinstance(left, URIRef) or not isinstance(right, URIRef):
                continue
            if left in selected and right in classes and iri_in_prefixes(right, prefixes):
                selected.add(right)
            if right in selected and left in classes and iri_in_prefixes(left, prefixes):
                selected.add(left)
        changed = len(selected) != previous


def validate_skos_source_identity(source: Graph, source_path: str | Path) -> None:
    """Require the configured SKOS file to identify the W3C SKOS ontology."""
    ontology_iris = {
        subject
        for subject in source.subjects(RDF.type, OWL.Ontology)
        if isinstance(subject, URIRef)
    }
    if SKOS_ONTOLOGY_IRI not in ontology_iris:
        declared = ", ".join(sorted(map(str, ontology_iris))) or "none"
        raise RuntimeError(
            f"Configured SKOS source {source_path} does not declare the W3C SKOS "
            f"ontology IRI {SKOS_ONTOLOGY_IRI}; declared ontology IRIs: {declared}"
        )


def validate_biolink_biomarker_axioms(source: Graph, module: Graph) -> None:
    """Verify critical Biolink biomarker axioms survive strict extraction.

    Only named domain/range targets are relevant to this strict module policy.
    Anonymous targets are intentionally unsupported and therefore cause a
    focused failure rather than being silently omitted.
    """
    properties = (BIOLINK.biomarker_for, BIOLINK.has_biomarker)
    errors: list[str] = []
    for prop in properties:
        if not any(source.triples((prop, None, None))):
            errors.append(f"source does not define {prop}")
            continue
        for predicate in (RDFS.domain, RDFS.range):
            source_all = set(source.objects(prop, predicate))
            source_named = {value for value in source_all if isinstance(value, URIRef)}
            module_values = set(module.objects(prop, predicate))
            if not source_all:
                errors.append(f"source has no {predicate.n3()} axiom for {prop}")
            elif len(source_named) != len(source_all):
                errors.append(
                    f"source uses anonymous {predicate.n3()} target(s) for {prop}, "
                    "which cannot be retained by the strict named-axiom policy"
                )
            elif source_named != module_values:
                errors.append(
                    f"{prop} {predicate.n3()}: source={sorted(map(str, source_named))}, "
                    f"module={sorted(map(str, module_values))}"
                )

    inverse_pair = {
        (BIOLINK.biomarker_for, BIOLINK.has_biomarker),
        (BIOLINK.has_biomarker, BIOLINK.biomarker_for),
    }
    source_inverse = {
        (subject, obj)
        for subject in properties
        for obj in source.objects(subject, OWL.inverseOf)
        if (subject, obj) in inverse_pair
    }
    module_inverse = {
        (subject, obj)
        for subject in properties
        for obj in module.objects(subject, OWL.inverseOf)
        if (subject, obj) in inverse_pair
    }
    if not source_inverse:
        errors.append(
            "source declares no inverse axiom between biomarker_for and has_biomarker"
        )
    elif source_inverse != module_inverse:
        errors.append(
            "inverse axioms differ: "
            f"source={sorted((str(a), str(b)) for a, b in source_inverse)}, "
            f"module={sorted((str(a), str(b)) for a, b in module_inverse)}"
        )

    if errors:
        raise RuntimeError(
            "Biolink module validation failed; critical biomarker axioms were "
            "missing or not preserved:\n  " + "\n  ".join(errors)
        )


def validate_bao_target_properties(source: Graph, module: Graph) -> None:
    """Require the pinned BAO source to define both target predicates logically."""
    errors: list[str] = []
    for prop in (BAO_HAS_TARGET, BAO_IS_TARGET_FOR):
        if not any(source.triples((prop, None, None))) and not any(
            source.triples((None, prop, None))
        ):
            errors.append(f"source does not define {prop}")
        if (prop, RDF.type, OWL.ObjectProperty) not in module:
            errors.append(f"strict module does not declare {prop} as owl:ObjectProperty")
    if errors:
        raise RuntimeError(
            "BAO module validation failed; required target predicates were missing "
            "or not classified as object properties:\n  " + "\n  ".join(errors)
        )


def build_reasoning_module(
    spec: VocabularySpec,
    source_path: str,
    signature_graph: Graph,
    destination: Path,
    kg_version: str,
) -> ModuleResult:
    source = parse_graph(source_path)
    if spec.name == "skos":
        validate_skos_source_identity(source, source_path)
    entities = source_entity_sets(source)
    classes = entities["classes"]
    source_object_properties = entities["object_properties"]
    source_object_only_properties = entities["object_only_properties"]
    source_datatype_properties = entities["datatype_properties"]
    source_annotation_properties = entities["annotation_properties"]
    source_logical_properties = entities["logical_properties"]
    source_datatypes = entities["datatypes"]
    source_named_individuals = entities["named_individuals"]

    used_terms, used_predicates = collect_signature(signature_graph, spec.prefixes)
    selected_classes = {term for term in used_terms if term in classes}
    selected_object_properties: set[URIRef] = set()
    selected_datatype_properties: set[URIRef] = set()
    selected_annotation_properties: set[URIRef] = set()
    resolved_annotation_logical_overlaps: set[URIRef] = set()

    undeclared_used_terms = {
        term
        for term in used_terms
        if term not in classes
        and term not in source_object_properties
        and term not in source_datatype_properties
        and term not in source_annotation_properties
        and not any(source.triples((term, None, None)))
        and not any(source.triples((None, term, None)))
        and not any(source.triples((None, None, term)))
    }
    if undeclared_used_terms:
        preview = "\n  ".join(map(str, sorted(undeclared_used_terms, key=str)[:30]))
        source_hint = ""
        if spec.name == "bao":
            source_hint = (
                f"\nConfigured BAO source: {source_path}"
                "\nRequired BAO source: "
                "data/external_ontologies/bao_complete_merged.owl containing "
                "BAO_0000211 and BAO_0000598. The modular "
                "bao_complete.owl entry file is not sufficient because imported "
                "BAO modules are not followed automatically."
            )
        raise RuntimeError(
            f"{spec.name}: {len(undeclared_used_terms)} used terms are absent from the pinned source:\n  "
            f"{preview}{source_hint}"
        )

    for predicate in used_predicates:
        if spec.name == "skos" and predicate in {SKOS_NOTATION, SKOS_PREF_LABEL}:
            selected_annotation_properties.add(predicate)
            continue

        is_object = predicate in source_object_properties
        is_datatype = predicate in source_datatype_properties
        is_annotation = predicate in source_annotation_properties
        is_logical = predicate in source_logical_properties
        usage = predicate_usage_kind(signature_graph, predicate)

        # Annotation/logical overlap is resolved into one strict logical
        # category when KG usage and source structure are decisive. The
        # annotation declaration is not copied into the OWL 2 DL module.
        if is_object and is_datatype:
            if usage == "resource":
                selected_object_properties.add(predicate)
            elif usage == "literal":
                selected_datatype_properties.add(predicate)
            else:
                raise RuntimeError(
                    f"{spec.name}: property {predicate} is both object and datatype in the source "
                    f"and has {usage or 'no'} decisive KG usage"
                )
        elif is_object:
            if usage in {"literal", "mixed"}:
                raise RuntimeError(
                    f"{spec.name}: object property {predicate} is used with "
                    f"{usage} values in the KG; object properties require resource values only"
                )
            selected_object_properties.add(predicate)
        elif is_datatype:
            if usage in {"resource", "mixed"}:
                raise RuntimeError(
                    f"{spec.name}: datatype property {predicate} is used with "
                    f"{usage} values in the KG; datatype properties require literal values only"
                )
            selected_datatype_properties.add(predicate)
        elif is_annotation and not is_logical:
            selected_annotation_properties.add(predicate)
        else:
            if usage == "resource":
                selected_object_properties.add(predicate)
            elif usage == "literal":
                if is_logical:
                    selected_datatype_properties.add(predicate)
                else:
                    selected_annotation_properties.add(predicate)
            else:
                raise RuntimeError(f"{spec.name}: cannot classify used predicate {predicate}")

        if predicate in selected_object_properties:
            assert_source_property_category_compatible(
                spec=spec,
                predicate=predicate,
                desired_category="object",
                source_object_properties=source_object_properties,
                source_object_only_properties=source_object_only_properties,
                source_datatype_properties=source_datatype_properties,
                source_annotation_properties=source_annotation_properties,
                source_logical_properties=source_logical_properties,
            )
        elif predicate in selected_datatype_properties:
            assert_source_property_category_compatible(
                spec=spec,
                predicate=predicate,
                desired_category="datatype",
                source_object_properties=source_object_properties,
                source_object_only_properties=source_object_only_properties,
                source_datatype_properties=source_datatype_properties,
                source_annotation_properties=source_annotation_properties,
                source_logical_properties=source_logical_properties,
            )
        elif predicate in selected_annotation_properties:
            assert_source_property_category_compatible(
                spec=spec,
                predicate=predicate,
                desired_category="annotation",
                source_object_properties=source_object_properties,
                source_object_only_properties=source_object_only_properties,
                source_datatype_properties=source_datatype_properties,
                source_annotation_properties=source_annotation_properties,
                source_logical_properties=source_logical_properties,
            )

        if (
            predicate in source_annotation_properties
            and predicate in source_logical_properties
            and predicate not in selected_annotation_properties
        ):
            resolved_annotation_logical_overlaps.add(predicate)

    directly_supported = (
        selected_classes
        | selected_object_properties
        | selected_datatype_properties
        | selected_annotation_properties
    )
    unsupported_used_terms = used_terms - directly_supported
    if unsupported_used_terms:
        preview = "\n  ".join(map(str, sorted(unsupported_used_terms, key=str)[:30]))
        raise RuntimeError(
            f"{spec.name}: {len(unsupported_used_terms)} directly used terms are not "
            "supported classes or properties in the strict module policy:\n  "
            + preview
        )

    expand_class_closure(source, selected_classes, classes, spec.prefixes)

    # Object-property closure. A property reached from an object property stays
    # in the object category when the source also supplies independent
    # object-only structural evidence; otherwise conflicting source categories
    # cause a focused failure.
    changed = True
    while changed:
        previous = len(selected_object_properties)
        for child, parent in source.subject_objects(RDFS.subPropertyOf):
            if (
                child in selected_object_properties
                and isinstance(parent, URIRef)
                and iri_in_prefixes(parent, spec.prefixes)
            ):
                assert_source_property_category_compatible(
                    spec=spec,
                    predicate=parent,
                    desired_category="object",
                    source_object_properties=source_object_properties,
                    source_object_only_properties=source_object_only_properties,
                    source_datatype_properties=source_datatype_properties,
                    source_annotation_properties=source_annotation_properties,
                    source_logical_properties=source_logical_properties,
                )
                selected_object_properties.add(parent)
        for relation in (OWL.inverseOf, OWL.equivalentProperty):
            for left, right in source.subject_objects(relation):
                if not isinstance(left, URIRef) or not isinstance(right, URIRef):
                    continue
                if not (iri_in_prefixes(left, spec.prefixes) and iri_in_prefixes(right, spec.prefixes)):
                    continue
                if left in selected_object_properties or right in selected_object_properties:
                    for endpoint in (left, right):
                        assert_source_property_category_compatible(
                            spec=spec,
                            predicate=endpoint,
                            desired_category="object",
                            source_object_properties=source_object_properties,
                            source_object_only_properties=source_object_only_properties,
                            source_datatype_properties=source_datatype_properties,
                            source_annotation_properties=source_annotation_properties,
                            source_logical_properties=source_logical_properties,
                        )
                    selected_object_properties.update((left, right))
        changed = len(selected_object_properties) != previous

    # Datatype-property closure. Preserve named superproperties and named
    # equivalent datatype properties, while rejecting object/data overlap.
    changed = True
    while changed:
        previous = len(selected_datatype_properties)
        for child, parent in source.subject_objects(RDFS.subPropertyOf):
            if (
                child in selected_datatype_properties
                and isinstance(parent, URIRef)
                and iri_in_prefixes(parent, spec.prefixes)
            ):
                if parent in selected_object_properties:
                    raise RuntimeError(
                        f"{spec.name}: datatype superproperty {parent} is already selected as object property"
                    )
                assert_source_property_category_compatible(
                    spec=spec,
                    predicate=parent,
                    desired_category="datatype",
                    source_object_properties=source_object_properties,
                    source_object_only_properties=source_object_only_properties,
                    source_datatype_properties=source_datatype_properties,
                    source_annotation_properties=source_annotation_properties,
                    source_logical_properties=source_logical_properties,
                )
                selected_datatype_properties.add(parent)
        for left, right in source.subject_objects(OWL.equivalentProperty):
            if not isinstance(left, URIRef) or not isinstance(right, URIRef):
                continue
            if not (
                iri_in_prefixes(left, spec.prefixes)
                and iri_in_prefixes(right, spec.prefixes)
            ):
                continue
            if left in selected_datatype_properties or right in selected_datatype_properties:
                if left in selected_object_properties or right in selected_object_properties:
                    raise RuntimeError(
                        f"{spec.name}: equivalent property pair mixes object and datatype categories: "
                        f"{left}, {right}"
                    )
                for endpoint in (left, right):
                    assert_source_property_category_compatible(
                        spec=spec,
                        predicate=endpoint,
                        desired_category="datatype",
                        source_object_properties=source_object_properties,
                        source_object_only_properties=source_object_only_properties,
                        source_datatype_properties=source_datatype_properties,
                        source_annotation_properties=source_annotation_properties,
                        source_logical_properties=source_logical_properties,
                    )
                selected_datatype_properties.update((left, right))
        changed = len(selected_datatype_properties) != previous

    selected_datatype_properties.difference_update(selected_object_properties)
    selected_annotation_properties.difference_update(selected_object_properties)
    selected_annotation_properties.difference_update(selected_datatype_properties)

    # Annotation-property closure follows named rdfs:subPropertyOf ancestors.
    changed = True
    while changed:
        previous = len(selected_annotation_properties)
        for child, parent in source.subject_objects(RDFS.subPropertyOf):
            if (
                child in selected_annotation_properties
                and isinstance(parent, URIRef)
                and iri_in_prefixes(parent, spec.prefixes)
            ):
                if parent in selected_object_properties or parent in selected_datatype_properties:
                    raise RuntimeError(
                        f"{spec.name}: annotation superproperty {parent} conflicts with a logical property category"
                    )
                assert_source_property_category_compatible(
                    spec=spec,
                    predicate=parent,
                    desired_category="annotation",
                    source_object_properties=source_object_properties,
                    source_object_only_properties=source_object_only_properties,
                    source_datatype_properties=source_datatype_properties,
                    source_annotation_properties=source_annotation_properties,
                    source_logical_properties=source_logical_properties,
                )
                selected_annotation_properties.add(parent)
        changed = len(selected_annotation_properties) != previous

    # Recheck actual KG usage after hierarchy/equivalence closure. A property
    # can enter a different category indirectly through a selected child or
    # equivalent property, so direct-signature checks alone are insufficient.
    for prop in selected_object_properties:
        usage = predicate_usage_kind(signature_graph, prop)
        if usage in {"literal", "mixed"}:
            raise RuntimeError(
                f"{spec.name}: property closure selected {prop} as an object property, "
                f"but its KG usage is {usage}; object properties require resource values only"
            )
    for prop in selected_datatype_properties:
        usage = predicate_usage_kind(signature_graph, prop)
        if usage in {"resource", "mixed"}:
            raise RuntimeError(
                f"{spec.name}: property closure selected {prop} as a datatype property, "
                f"but its KG usage is {usage}; datatype properties require literal values only"
            )

    # Reject named domain/range category contradictions instead of silently
    # omitting them from the strict module. Anonymous expressions remain excluded
    # by policy and are handled later by the OWL 2 DL validation of the result.
    for prop in selected_object_properties:
        for predicate in (RDFS.domain, RDFS.range):
            for target in source.objects(prop, predicate):
                if not isinstance(target, URIRef):
                    continue
                if (
                    target in BUILTIN_DATATYPES
                    or str(target).startswith(str(XSD))
                    or target in source_datatypes
                ):
                    raise RuntimeError(
                        f"{spec.name}: object property {prop} has datatype-valued "
                        f"{predicate} target {target} in the pinned source"
                    )
    for prop in selected_datatype_properties:
        for target in source.objects(prop, RDFS.domain):
            if not isinstance(target, URIRef):
                continue
            if (
                target in BUILTIN_DATATYPES
                or str(target).startswith(str(XSD))
                or target in source_datatypes
            ):
                raise RuntimeError(
                    f"{spec.name}: datatype property {prop} has datatype-valued domain {target}"
                )
        for target in source.objects(prop, RDFS.range):
            if not isinstance(target, URIRef):
                continue
            if (
                target not in BUILTIN_DATATYPES
                and not str(target).startswith(str(XSD))
                and target not in source_datatypes
            ):
                category = "class" if target in classes else "undeclared non-datatype"
                raise RuntimeError(
                    f"{spec.name}: datatype property {prop} has {category} range {target}; "
                    "a named datatype range is required by the strict module policy"
                )

    # Add named domain/range classes and their ancestors.
    for prop in selected_object_properties:
        for predicate in (RDFS.domain, RDFS.range):
            for target in source.objects(prop, predicate):
                if (
                    isinstance(target, URIRef)
                    and target in classes
                    and iri_in_prefixes(target, spec.prefixes)
                ):
                    selected_classes.add(target)
    for prop in selected_datatype_properties:
        for target in source.objects(prop, RDFS.domain):
            if (
                isinstance(target, URIRef)
                and target in classes
                and iri_in_prefixes(target, spec.prefixes)
            ):
                selected_classes.add(target)
    expand_class_closure(source, selected_classes, classes, spec.prefixes)

    selected_datatypes: set[URIRef] = set()
    for prop in selected_datatype_properties:
        for target in source.objects(prop, RDFS.range):
            if (
                isinstance(target, URIRef)
                and target in source_datatypes
                and target not in BUILTIN_DATATYPES
                and not str(target).startswith(str(XSD))
            ):
                selected_datatypes.add(target)

    datatype_entity_overlap = selected_datatypes & (
        selected_classes
        | selected_object_properties
        | selected_datatype_properties
        | selected_annotation_properties
    )
    if datatype_entity_overlap:
        preview = "\n  ".join(
            map(str, sorted(datatype_entity_overlap, key=str)[:30])
        )
        raise RuntimeError(
            f"{spec.name}: selected IRIs are used as both datatypes and classes/properties, "
            f"which is outside this OWL 2 DL module policy:\n  {preview}"
        )

    class_property_overlap = selected_classes & (
        selected_object_properties
        | selected_datatype_properties
        | selected_annotation_properties
    )
    if class_property_overlap:
        preview = "\n  ".join(map(str, sorted(class_property_overlap, key=str)[:30]))
        raise RuntimeError(
            f"{spec.name}: selected IRIs are used as both classes and properties, "
            f"which is outside this OWL 2 DL module policy:\n  {preview}"
        )

    module = Graph()
    module.namespace_manager = source.namespace_manager
    module.bind("owl", OWL)
    module.bind("rdfs", RDFS)
    module.bind("skos", SKOS)
    module.bind("dcterms", DCTERMS)
    source_hash = file_sha256(source_path)
    add_module_header(module, spec, source_path, source_hash, kg_version)

    for cls in selected_classes:
        module.add((cls, RDF.type, OWL.Class))
    for prop in selected_object_properties:
        module.add((prop, RDF.type, OWL.ObjectProperty))
    for prop in selected_datatype_properties:
        module.add((prop, RDF.type, OWL.DatatypeProperty))
    for prop in selected_annotation_properties:
        module.add((prop, RDF.type, OWL.AnnotationProperty))
    for datatype in selected_datatypes:
        module.add((datatype, RDF.type, RDFS.Datatype))

    for child, parent in source.subject_objects(RDFS.subClassOf):
        if child in selected_classes and parent in selected_classes:
            module.add((child, RDFS.subClassOf, parent))
    for predicate in (OWL.equivalentClass, OWL.disjointWith):
        for left, right in source.subject_objects(predicate):
            if left in selected_classes and right in selected_classes:
                module.add((left, predicate, right))

    for predicate in (RDFS.subPropertyOf, OWL.inverseOf, OWL.equivalentProperty):
        for left, right in source.subject_objects(predicate):
            if left in selected_object_properties and right in selected_object_properties:
                module.add((left, predicate, right))
    for child, parent in source.subject_objects(RDFS.subPropertyOf):
        if child in selected_datatype_properties and parent in selected_datatype_properties:
            module.add((child, RDFS.subPropertyOf, parent))
        if child in selected_annotation_properties and parent in selected_annotation_properties:
            module.add((child, RDFS.subPropertyOf, parent))
    for left, right in source.subject_objects(OWL.equivalentProperty):
        if left in selected_datatype_properties and right in selected_datatype_properties:
            module.add((left, OWL.equivalentProperty, right))

    for prop in selected_object_properties:
        for characteristic in OBJECT_PROPERTY_CHARACTERISTICS:
            if (prop, RDF.type, characteristic) in source:
                module.add((prop, RDF.type, characteristic))
        for predicate in (RDFS.domain, RDFS.range):
            for target in source.objects(prop, predicate):
                if isinstance(target, URIRef) and target in selected_classes:
                    module.add((prop, predicate, target))

    for prop in selected_datatype_properties:
        if (prop, RDF.type, OWL.FunctionalProperty) in source:
            module.add((prop, RDF.type, OWL.FunctionalProperty))
        for target in source.objects(prop, RDFS.domain):
            if isinstance(target, URIRef) and target in selected_classes:
                module.add((prop, RDFS.domain, target))
        for target in source.objects(prop, RDFS.range):
            if (
                isinstance(target, URIRef)
                and (
                    target in BUILTIN_DATATYPES
                    or str(target).startswith(str(XSD))
                    or target in selected_datatypes
                )
            ):
                module.add((prop, RDFS.range, target))

    resolved_annotation_logical_overlaps.update(
        (selected_object_properties | selected_datatype_properties)
        & source_annotation_properties
        & source_logical_properties
    )

    selected_entities = (
        selected_classes
        | selected_object_properties
        | selected_datatype_properties
        | selected_annotation_properties
    )
    annotation_predicates: set[URIRef] = set(selected_annotation_properties)
    annotation_literal_datatypes: set[URIRef] = set()
    for entity in selected_entities:
        for predicate, value in source.predicate_objects(entity):
            if predicate in LOGICAL_PREDICATES or not isinstance(value, Literal):
                continue
            if predicate in selected_object_properties or predicate in selected_datatype_properties:
                continue
            if (
                predicate in source_object_properties
                or predicate in source_datatype_properties
                or predicate in source_logical_properties
            ):
                continue
            module.add((entity, predicate, value))
            if isinstance(predicate, URIRef):
                annotation_predicates.add(predicate)
            if (
                isinstance(value.datatype, URIRef)
                and not is_builtin_datatype_iri(value.datatype)
            ):
                annotation_literal_datatypes.add(value.datatype)

    for predicate in annotation_predicates:
        if predicate not in selected_object_properties and predicate not in selected_datatype_properties:
            module.add((predicate, RDF.type, OWL.AnnotationProperty))

    class_annotation_overlap = selected_classes & annotation_predicates
    if class_annotation_overlap:
        preview = "\n  ".join(
            map(str, sorted(class_annotation_overlap, key=str)[:30])
        )
        raise RuntimeError(
            f"{spec.name}: source annotation copying would use selected classes as "
            f"annotation properties, which is outside this OWL 2 DL module policy:\n  {preview}"
        )

    datatype_annotation_overlap = (selected_datatypes | annotation_literal_datatypes) & annotation_predicates
    if datatype_annotation_overlap:
        preview = "\n  ".join(
            map(str, sorted(datatype_annotation_overlap, key=str)[:30])
        )
        raise RuntimeError(
            f"{spec.name}: source annotation copying would use the same IRI as a "
            f"custom datatype and annotation property:\n  {preview}"
        )

    datatype_entity_conflicts = annotation_literal_datatypes & (
        classes
        | source_object_properties
        | source_datatype_properties
        | source_annotation_properties
        | source_logical_properties
        | source_named_individuals
        | selected_classes
        | selected_object_properties
        | selected_datatype_properties
        | selected_annotation_properties
    )
    if datatype_entity_conflicts:
        preview = "\n  ".join(
            map(str, sorted(datatype_entity_conflicts, key=str)[:30])
        )
        raise RuntimeError(
            f"{spec.name}: literal annotations use datatype IRIs that are also "
            f"classes or properties in the pinned source/module:\n  {preview}"
        )

    selected_datatypes.update(annotation_literal_datatypes)
    for datatype in annotation_literal_datatypes:
        module.add((datatype, RDF.type, RDFS.Datatype))

    if spec.name == "biolink" and {
        BIOLINK.biomarker_for,
        BIOLINK.has_biomarker,
    } & used_terms:
        validate_biolink_biomarker_axioms(source, module)

    if spec.name == "bao" and {BAO_HAS_TARGET, BAO_IS_TARGET_FOR} & used_terms:
        validate_bao_target_properties(source, module)

    # Strong structural postconditions for the strict module.
    if any(module.subjects(RDF.type, OWL.Restriction)):
        raise RuntimeError(f"{spec.name}: strict module unexpectedly contains owl:Restriction")
    for predicate in (
        RDFS.subClassOf,
        OWL.equivalentClass,
        OWL.disjointWith,
        RDFS.subPropertyOf,
        OWL.inverseOf,
        OWL.equivalentProperty,
        RDFS.domain,
        RDFS.range,
    ):
        for subject, obj in module.subject_objects(predicate):
            if isinstance(subject, BNode) or isinstance(obj, BNode):
                raise RuntimeError(f"{spec.name}: blank-node logical axiom remains: {(subject, predicate, obj)}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    module.serialize(destination=str(destination), format="xml")
    verified = Graph()
    verified.parse(str(destination), format="xml")
    if len(verified) != len(module):
        raise RuntimeError(f"{spec.name}: serialization round-trip changed triple count")

    output_terms = {
        node
        for triple in verified
        for node in triple
        if iri_in_prefixes(node, spec.prefixes)
    }
    missing_output = used_terms - output_terms
    if missing_output:
        preview = "\n  ".join(map(str, sorted(missing_output, key=str)[:30]))
        raise RuntimeError(f"{spec.name}: used terms missing from strict module:\n  {preview}")

    retained_axiom_counts = {
        "class_declarations": len(set(verified.subjects(RDF.type, OWL.Class))),
        "object_property_declarations": len(set(verified.subjects(RDF.type, OWL.ObjectProperty))),
        "datatype_property_declarations": len(set(verified.subjects(RDF.type, OWL.DatatypeProperty))),
        "annotation_property_declarations": len(set(verified.subjects(RDF.type, OWL.AnnotationProperty))),
        "datatype_declarations": len(set(verified.subjects(RDF.type, RDFS.Datatype))),
        "subclass_axioms": len(list(verified.triples((None, RDFS.subClassOf, None)))),
        "equivalent_class_axioms": len(list(verified.triples((None, OWL.equivalentClass, None)))),
        "disjoint_class_axioms": len(list(verified.triples((None, OWL.disjointWith, None)))),
        "subproperty_axioms": len(list(verified.triples((None, RDFS.subPropertyOf, None)))),
        "inverse_property_axioms": len(list(verified.triples((None, OWL.inverseOf, None)))),
        "equivalent_property_axioms": len(list(verified.triples((None, OWL.equivalentProperty, None)))),
        "domain_axioms": len(list(verified.triples((None, RDFS.domain, None)))),
        "range_axioms": len(list(verified.triples((None, RDFS.range, None)))),
        "restrictions": len(set(verified.subjects(RDF.type, OWL.Restriction))),
    }

    result = ModuleResult(
        name=spec.name,
        path=str(Path("modules") / destination.name),
        source_path=portable_source_path(source_path),
        source_sha256=source_hash,
        output_sha256=file_sha256(destination),
        used_terms=len(used_terms),
        used_predicates=len(used_predicates),
        selected_classes=len(selected_classes),
        selected_object_properties=len(selected_object_properties),
        selected_datatype_properties=len(selected_datatype_properties),
        selected_annotation_properties=len(annotation_predicates),
        selected_datatypes=len(selected_datatypes),
        triples=len(verified),
        direct_term_iris=tuple(sorted(map(str, used_terms))),
        direct_predicate_iris=tuple(sorted(map(str, used_predicates))),
        selected_class_iris=tuple(sorted(map(str, selected_classes))),
        selected_object_property_iris=tuple(sorted(map(str, selected_object_properties))),
        selected_datatype_property_iris=tuple(sorted(map(str, selected_datatype_properties))),
        selected_annotation_property_iris=tuple(sorted(map(str, annotation_predicates))),
        selected_datatype_iris=tuple(sorted(map(str, selected_datatypes))),
        infrastructure_annotation_property_iris=tuple(
            sorted(map(str, {DCTERMS.source, DCTERMS.identifier}))
        ),
        resolved_annotation_logical_overlap_iris=tuple(
            sorted(map(str, resolved_annotation_logical_overlaps))
        ),
        retained_axiom_counts=retained_axiom_counts,
    )
    # Close large RDFLib stores before processing the next pinned ontology.
    # The serialized module and immutable audit result are already complete.
    verified.close()
    module.close()
    source.close()
    return result


def write_catalog(path: Path, results: Iterable[ModuleResult]) -> None:
    """Write and verify the reasoning-local XML catalog."""
    result_list = list(results)
    result_by_name = {result.name: result for result in result_list}
    expected_names = {spec.name for spec in VOCABULARIES}
    if len(result_by_name) != len(result_list) or set(result_by_name) != expected_names:
        raise RuntimeError(
            "Reasoning module results do not match the configured vocabulary set: "
            f"expected={sorted(expected_names)}, actual={sorted(result_by_name)}"
        )

    root = ET.Element(
        "catalog",
        {
            "xmlns": "urn:oasis:names:tc:entity:xmlns:xml:catalog",
            "prefer": "public",
        },
    )
    expected_mappings: dict[str, str] = {}
    for spec in VOCABULARIES:
        result = result_by_name[spec.name]
        relative_target = Path(result.path).as_posix()
        expected_mappings[str(spec.module_iri)] = relative_target
        ET.SubElement(
            root,
            "uri",
            {"name": str(spec.module_iri), "uri": relative_target},
        )
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)

    parsed = ET.parse(path).getroot()
    found: dict[str, str] = {}
    for element in parsed:
        if not element.tag.endswith("uri"):
            continue
        name = element.get("name")
        target = element.get("uri")
        if name and target:
            found[name] = target
    if found != expected_mappings:
        raise RuntimeError(
            f"Reasoning catalog mappings are inconsistent: expected={expected_mappings}, "
            f"found={found}"
        )
    for target in found.values():
        resolved = path.parent / target
        if not resolved.is_file() or resolved.stat().st_size == 0:
            raise RuntimeError(f"Reasoning catalog target is missing or empty: {resolved}")


def validate_expected_imports(tbox_graph: Graph) -> dict[str, list[str]]:
    """Require the standalone TBox imports to match the generated module set."""
    expected = {spec.module_iri for spec in VOCABULARIES}
    import_triples = list(tbox_graph.triples((None, OWL.imports, None)))
    ontology_headers = set(tbox_graph.subjects(RDF.type, OWL.Ontology))
    invalid_subjects = sorted(
        map(
            str,
            {
                subject
                for subject, _, _ in import_triples
                if subject not in ontology_headers
            },
        )
    )
    actual_nodes = {obj for _, _, obj in import_triples}
    non_iri = sorted(map(str, (node for node in actual_nodes if not isinstance(node, URIRef))))
    actual = {node for node in actual_nodes if isinstance(node, URIRef)}
    missing = sorted(map(str, expected - actual))
    unexpected = sorted(map(str, actual - expected))
    if invalid_subjects or non_iri or missing or unexpected:
        details = []
        if invalid_subjects:
            details.append(
                "imports asserted by non-ontology subjects:\n  "
                + "\n  ".join(invalid_subjects)
            )
        if non_iri:
            details.append("non-IRI imports:\n  " + "\n  ".join(non_iri))
        if missing:
            details.append("missing expected imports:\n  " + "\n  ".join(missing))
        if unexpected:
            details.append("unsupported imports that would otherwise be dropped:\n  " + "\n  ".join(unexpected))
        raise RuntimeError(
            "Standalone MAKAAO TBox import set does not match the strict module set:\n"
            + "\n".join(details)
        )
    return {
        "expected": sorted(map(str, expected)),
        "actual": sorted(map(str, actual)),
    }


def ontology_subjects(graph: Graph) -> set[URIRef | BNode]:
    return set(graph.subjects(RDF.type, OWL.Ontology))


def copy_without_headers_and_imports(source: Graph, destination: Graph) -> None:
    headers = ontology_subjects(source)
    for subject, predicate, obj in source:
        if predicate == OWL.imports or subject in headers:
            continue
        destination.add((subject, predicate, obj))


def add_reasoning_metadata_declarations(graph: Graph) -> None:
    # Any URI used as an rdf:type object is a class in the RDF graph. Add the
    # missing declarations needed for the OWL 2 DL parser, except for OWL/RDF
    # metatypes and property characteristics.
    non_class_types = {
        OWL.Ontology,
        OWL.Class,
        RDFS.Class,
        OWL.ObjectProperty,
        OWL.DatatypeProperty,
        OWL.AnnotationProperty,
        OWL.NamedIndividual,
        RDF.Property,
        RDFS.Datatype,
        RDF.Statement,
        *OBJECT_PROPERTY_CHARACTERISTICS,
    }
    declared_classes = set(graph.subjects(RDF.type, OWL.Class)) | set(
        graph.subjects(RDF.type, RDFS.Class)
    )
    for cls in set(graph.objects(None, RDF.type)):
        if not isinstance(cls, URIRef):
            continue
        if str(cls).startswith((str(RDF), str(RDFS), str(OWL), str(XSD))):
            continue
        if cls not in non_class_types and cls not in declared_classes:
            graph.add((cls, RDF.type, OWL.Class))

    # Dataset/contact/rights metadata is non-logical for the reasoning release.
    # Declare undeclared predicates from those vocabularies as annotations.
    declared_properties = {
        prop
        for prop_type in (OWL.ObjectProperty, OWL.DatatypeProperty, OWL.AnnotationProperty)
        for prop in graph.subjects(RDF.type, prop_type)
        if isinstance(prop, URIRef)
    }

    # OWL 2 DL requires non-built-in literal datatype IRIs to be declared and
    # disjoint from class/property/named-individual entity categories.  Module
    # extraction normally supplies these declarations, but this graph-level
    # pass also covers custom datatypes used only in ABox or metadata literals.
    custom_literal_datatypes = {
        obj.datatype
        for obj in graph.objects()
        if isinstance(obj, Literal)
        and isinstance(obj.datatype, URIRef)
        and not is_builtin_datatype_iri(obj.datatype)
    }
    datatype_conflicting_entities = {
        entity
        for entity_type in (
            OWL.Class,
            RDFS.Class,
            OWL.ObjectProperty,
            OWL.DatatypeProperty,
            OWL.AnnotationProperty,
            OWL.NamedIndividual,
        )
        for entity in graph.subjects(RDF.type, entity_type)
        if isinstance(entity, URIRef)
    }
    builtin_predicates = {
        RDF.type,
        RDFS.label,
        RDFS.comment,
        RDFS.seeAlso,
        RDFS.isDefinedBy,
        OWL.versionInfo,
        OWL.versionIRI,
        OWL.imports,
        SKOS.prefLabel,
        SKOS.altLabel,
        SKOS.hiddenLabel,
    }
    metadata_annotation_candidates = {
        predicate
        for predicate in graph.predicates()
        if isinstance(predicate, URIRef)
        and predicate not in declared_properties
        and predicate not in builtin_predicates
        and str(predicate).startswith(METADATA_NAMESPACES)
    }
    datatype_conflicts = custom_literal_datatypes & (
        datatype_conflicting_entities | metadata_annotation_candidates
    )
    if datatype_conflicts:
        preview = "\n  ".join(map(str, sorted(datatype_conflicts, key=str)[:30]))
        raise RuntimeError(
            "Custom literal datatype IRIs overlap class/property/individual or "
            f"metadata-annotation entity categories in the DL graph:\n  {preview}"
        )
    for datatype in custom_literal_datatypes:
        graph.add((datatype, RDF.type, RDFS.Datatype))

    for predicate in metadata_annotation_candidates:
        graph.add((predicate, RDF.type, OWL.AnnotationProperty))

    # Integration-specific SKOS profile used by MAKAAO.
    graph.remove((SKOS_NOTATION, RDF.type, OWL.DatatypeProperty))
    graph.remove((SKOS_NOTATION, RDFS.domain, None))
    graph.remove((SKOS_NOTATION, RDFS.range, None))
    graph.add((SKOS_NOTATION, RDF.type, OWL.AnnotationProperty))


def remove_generic_reification(graph: Graph) -> Graph:
    """Remove RDF/MAKAAO reification record instances, retaining schema classes."""
    reification_predicates = {RDF.subject, RDF.predicate, RDF.object}
    reification_record_classes = {RDF.Statement, MAKAAO.Relation, MAKAAO.Document}
    nodes: set[Any] = set()
    for cls in reification_record_classes:
        nodes.update(graph.subjects(RDF.type, cls))
    for predicate in reification_predicates:
        nodes.update(graph.subjects(predicate, None))

    output = Graph()
    output.namespace_manager = graph.namespace_manager
    for subject, predicate, obj in graph:
        if subject in nodes or obj in nodes:
            continue
        if predicate in reification_predicates:
            continue
        output.add((subject, predicate, obj))
    return output

def graph_postconditions(graph: Graph, *, require_no_imports: bool) -> dict[str, int]:
    """Check structural conditions promised for retained reasoning artifacts."""
    reification_predicates = {RDF.subject, RDF.predicate, RDF.object}
    statement_nodes = set(graph.subjects(RDF.type, RDF.Statement))
    relation_nodes = set(graph.subjects(RDF.type, MAKAAO.Relation))
    document_nodes = set(graph.subjects(RDF.type, MAKAAO.Document))
    reification_triples = sum(
        1
        for predicate in reification_predicates
        for _ in graph.triples((None, predicate, None))
    )
    imports = len(list(graph.triples((None, OWL.imports, None))))

    if statement_nodes or relation_nodes or document_nodes or reification_triples:
        raise RuntimeError(
            "Reification/provenance records remain in a reasoning artifact: "
            f"rdf:Statement instances={len(statement_nodes)}, "
            f"makaao:Relation instances={len(relation_nodes)}, "
            f"makaao:Document instances={len(document_nodes)}, "
            f"rdf:subject/predicate/object triples={reification_triples}"
        )
    if require_no_imports and imports:
        raise RuntimeError(f"Reasoning artifact still contains {imports} owl:imports triples")

    return {
        "rdf_statement_instances": len(statement_nodes),
        "makaao_relation_instances": len(relation_nodes),
        "makaao_document_instances": len(document_nodes),
        "rdf_reification_predicate_triples": reification_triples,
        "owl_imports": imports,
    }

def serialize_verified(graph: Graph, path: Path) -> Graph:
    path.parent.mkdir(parents=True, exist_ok=True)
    graph.serialize(destination=str(path), format="xml")
    verified = Graph()
    try:
        verified.parse(str(path), format="xml")
        if len(verified) != len(graph):
            raise RuntimeError(f"Serialization round-trip changed triple count for {path}")
        return verified
    except BaseException:
        verified.close()
        raise


def run_command_tail(
    command: Sequence[str],
    *,
    timeout_seconds: int,
    log_path: Path,
    tail_lines: int = 200,
    env: dict[str, str] | None = None,
) -> tuple[int, str]:
    """Run a command with a real wall-clock timeout and retain only its tail."""
    tail: deque[str] = deque(maxlen=tail_lines)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    process = subprocess.Popen(
        list(command),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
        env=env,
        start_new_session=(os.name == "posix"),
    )

    def _read_stdout() -> None:
        assert process.stdout is not None
        try:
            for line in process.stdout:
                tail.append(line.rstrip("\n"))
        except (OSError, ValueError):
            # The main thread may close the stream after a timeout if a wrapper
            # or descendant process keeps the pipe open.
            pass

    def _kill_process_group() -> None:
        if process.poll() is not None:
            return
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:  # pragma: no cover - Windows fallback
            process.kill()

    reader = threading.Thread(target=_read_stdout, daemon=True)
    reader.start()
    timed_out = False
    try:
        return_code = process.wait(timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        timed_out = True
        _kill_process_group()
        return_code = process.wait()
    except BaseException:
        # Ctrl-C or SystemExit must not leave HermiT/Java running after the
        # Python build has stopped. The finally block still records the bounded
        # output tail before the original exception is re-raised.
        _kill_process_group()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:  # pragma: no cover - SIGKILL fallback
            process.kill()
            process.wait()
        raise
    finally:
        reader.join(timeout=10)
        if reader.is_alive() and process.stdout is not None:
            process.stdout.close()
            reader.join(timeout=2)
        if process.stdout is not None and not process.stdout.closed:
            process.stdout.close()
        with log_path.open("w", encoding="utf-8") as log:
            log.write("\n".join(tail) + ("\n" if tail else ""))

    tail_text = "\n".join(tail)
    if timed_out:
        raise RuntimeError(
            f"Command timed out after {timeout_seconds}s: {' '.join(command)}"
            + (f"\nLast output lines:\n{tail_text}" if tail_text else "")
        )
    return return_code, tail_text


def validate_profile(
    robot_command: Sequence[str],
    ontology_path: Path,
    report_path: Path,
    logs_dir: Path,
    timeout_seconds: int,
    subprocess_env: dict[str, str] | None,
) -> None:
    command = [
        *robot_command,
        "--strict",
        "validate-profile",
        "--profile",
        "DL",
        "--input",
        str(ontology_path),
        "--output",
        str(report_path),
    ]
    return_code, tail = run_command_tail(
        command,
        timeout_seconds=timeout_seconds,
        log_path=logs_dir / f"{ontology_path.stem}-profile.log",
        env=subprocess_env,
    )
    if return_code != 0:
        raise RuntimeError(
            f"OWL 2 DL validation failed for {ontology_path}.\n"
            f"Command: {' '.join(command)}\nLast output lines:\n{tail}"
        )
    if not report_path.is_file():
        raise RuntimeError(f"ROBOT did not create profile report {report_path}")


def run_reasoner(
    robot_command: Sequence[str],
    ontology_path: Path,
    output_path: Path,
    logs_dir: Path,
    timeout_seconds: int,
    reasoner_name: str,
    subprocess_env: dict[str, str] | None,
) -> None:
    command = [
        *robot_command,
        "reason",
        "--reasoner",
        reasoner_name,
        "--input",
        str(ontology_path),
        "--output",
        str(output_path),
    ]
    return_code, tail = run_command_tail(
        command,
        timeout_seconds=timeout_seconds,
        log_path=logs_dir / f"{ontology_path.stem}-{reasoner_name.lower()}.log",
        env=subprocess_env,
    )
    if return_code != 0:
        raise RuntimeError(
            f"{reasoner_name} failed for {ontology_path}.\n"
            f"Command: {' '.join(command)}\nLast output lines:\n{tail}"
        )
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Reasoner did not create output {output_path}")



def preserve_asserted_axioms(
    asserted_path: Path,
    reasoned_path: Path,
) -> dict[str, int]:
    """Make the retained reasoned graph an assertion-preserving union.

    ROBOT's ``reason`` command may normalize a classified ontology by omitting
    asserted schema axioms that remain logically entailed.  The MAKAAO release
    contract is stricter: the retained reasoned file must contain every RDF
    triple from the asserted DL input plus the inferred triples produced by the
    reasoner.  The current strict release is named-node-only; blank-node axioms
    are rejected here because their parser-local identifiers cannot support a
    reliable triple-by-triple preservation check across serializations.
    """
    asserted = Graph()
    reasoned = Graph()
    verified: Graph | None = None
    try:
        asserted.parse(str(asserted_path), format="xml")
        reasoned.parse(str(reasoned_path), format="xml")

        asserted_blank_nodes = {
            node
            for triple in asserted
            for node in triple
            if isinstance(node, BNode)
        }
        reasoned_blank_nodes = {
            node
            for triple in reasoned
            for node in triple
            if isinstance(node, BNode)
        }
        if asserted_blank_nodes or reasoned_blank_nodes:
            raise RuntimeError(
                "Assertion-preserving reasoned publication requires named-node-only "
                "graphs; blank nodes were found in the asserted or reasoned ontology"
            )

        reasoner_output_triples = len(reasoned)
        missing = [triple for triple in asserted if triple not in reasoned]
        for triple in missing:
            reasoned.add(triple)

        verified = serialize_verified(reasoned, reasoned_path)
        remaining = [triple for triple in asserted if triple not in verified]
        if remaining:
            preview = "\n  ".join(map(str, remaining[:20]))
            raise RuntimeError(
                f"Reasoned ontology is missing {len(remaining)} asserted triples "
                f"after assertion restoration:\n  {preview}"
            )

        return {
            "asserted_input_triples": len(asserted),
            "reasoner_output_triples_before_restore": reasoner_output_triples,
            "asserted_triples_restored": len(missing),
            "asserted_triples_missing_after_restore": len(remaining),
            "reasoned_triples_after_restore": len(verified),
        }
    finally:
        asserted.close()
        reasoned.close()
        if verified is not None:
            verified.close()

def check_reasoned_output(path: Path) -> dict[str, Any]:
    graph = Graph()
    graph.parse(str(path), format="xml")
    structural = graph_postconditions(graph, require_no_imports=True)
    unsatisfiable: set[URIRef] = set()
    for left, right in graph.subject_objects(OWL.equivalentClass):
        if right == OWL.Nothing and isinstance(left, URIRef):
            unsatisfiable.add(left)
        if left == OWL.Nothing and isinstance(right, URIRef):
            unsatisfiable.add(right)
    for cls in graph.subjects(RDFS.subClassOf, OWL.Nothing):
        if isinstance(cls, URIRef):
            unsatisfiable.add(cls)
    unsatisfiable.discard(OWL.Nothing)
    nothing_instances = {
        individual
        for individual in graph.subjects(RDF.type, OWL.Nothing)
        if isinstance(individual, (URIRef, BNode))
    }
    result = {
        "triples": len(graph),
        "unsatisfiable_named_classes": sorted(map(str, unsatisfiable)),
        "individuals_typed_owl_nothing": sorted(map(str, nothing_instances)),
        "structural_postconditions": structural,
    }
    graph.close()
    return result


def extraction_policy() -> dict[str, Any]:
    """Machine-readable definition of the strict module extraction policy."""
    return {
        "signature_source": "complete generated MAKAAO KG before reification removal",
        "direct_signature": {
            "terms": "IRIs in subject, predicate, or object position matching the vocabulary namespace",
            "predicates": "matching IRIs used in predicate position",
        },
        "retained": [
            "declarations for directly used classes and properties",
            "named superclass ancestors within the same vocabulary namespace",
            "named equivalent classes reached from selected classes",
            "named superproperties within the same vocabulary namespace",
            "named inverse/equivalent object properties and named equivalent datatype properties",
            "named domain and range classes plus their named superclass ancestors",
            "named custom datatypes used as retained datatype-property ranges or retained literal annotations",
            "named subclass, equivalence, disjointness, property hierarchy, inverse, domain and range axioms when both endpoints are retained",
            "OWL object-property characteristics and functional datatype-property characteristics present on retained properties",
            "literal annotations on retained entities when the annotation predicate is not an object or datatype property",
        ],
        "excluded": [
            "anonymous class expressions and owl:Restriction structures",
            "RDF lists and n-ary OWL axiom structures",
            "property chains and keys",
            "URI-valued mapping annotations not required as logical axioms",
            "entities and ontology branches outside the generated KG signature and its named closure",
            "unresolved object/datatype property overlap; category is resolved only when KG usage, named range evidence, or OWL object-only structural evidence is decisive",
            "annotation/logical dual declarations are projected to the decisive logical category; the conflicting annotation declaration is omitted and audited",
            "any category selected through hierarchy/equivalence closure whose actual KG value kind conflicts with that category",
            "named object/datatype property domain or range axioms whose target category contradicts the selected property category",
        ],
        "source_integrity_checks": {
            "SKOS": "the configured file must declare the W3C SKOS ontology IRI",
            "Biolink biomarker axioms": "when biomarker_for/has_biomarker are selected, named domains, ranges, and their inverse relation must be present and preserved",
            "imports": "all nine imports must be IRI-valued and asserted by owl:Ontology subjects",
        },
        "integration_overrides": {
            "skos:notation": "owl:AnnotationProperty",
            "skos:prefLabel": "owl:AnnotationProperty",
            "generic RDF and MAKAAO reification records": "excluded from the DL graph, retained in the original KG; schema class declarations may remain",
            "dataset metadata vocabularies": "minimal declarations added directly to the DL graph; no full metadata ontology modules",
            "object-only OWL structures": "symmetric/transitive/inverse-functional/asymmetric/reflexive/irreflexive types, inverseOf endpoints, and property-chain subjects provide object-property evidence; property chains themselves are not retained",
        },
        "formal_module_type": "deterministic at the selected RDF-graph level, signature-driven syntactic module; not a formal locality module or canonical byte serialization",
    }


def write_extraction_documentation(output_root: Path) -> None:
    policy = extraction_policy()
    lines = [
        "MAKAAO strict reasoning module extraction policy",
        "=" * 49,
        "",
        "Source",
        "------",
        "Each module is extracted directly from the pinned file in data/external_ontologies.",
        "No previously extracted import module is used as source material.",
        "",
        "Direct signature",
        "----------------",
        "For each vocabulary namespace, every matching IRI occurring as subject, predicate,",
        "or object in the complete generated MAKAAO KG is a direct term. Matching IRIs used",
        "in predicate position form the direct predicate signature.",
        "",
        "Retained",
        "--------",
    ]
    lines.extend(f"- {item}" for item in policy["retained"])
    lines.extend(["", "Excluded", "--------"])
    lines.extend(f"- {item}" for item in policy["excluded"])
    lines.extend(["", "Source-integrity checks", "-----------------------"])
    lines.extend(f"- {key}: {value}" for key, value in policy["source_integrity_checks"].items())
    lines.extend(["", "Integration overrides", "---------------------"])
    lines.extend(f"- {key}: {value}" for key, value in policy["integration_overrides"].items())
    lines.extend([
        "",
        "Module type",
        "-----------",
        policy["formal_module_type"],
        "",
        "Exact selections",
        "----------------",
        "reports/<vocabulary>-selection.json lists every direct IRI, every selected IRI,",
        "the pinned source path and SHA-256 checksum, retained axiom counts, and the output checksum.",
    ])
    (output_root / "EXTRACTION_POLICY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_reasoning_release(
    *,
    kg_graph: Graph,
    tbox_graph: Graph,
    non_reified_graph: Graph,
    output_dir: str | Path,
    source_files: dict[str, str],
    skos_source_file: str,
    kg_version: str,
    robot_command: Sequence[str] | None,
    reasoner_name: str = "HermiT",
    timeout_seconds: int = 1800,
    require_robot: bool = True,
    canonical_kg_path: str | Path | None = None,
    canonical_tbox_path: str | Path | None = None,
    subprocess_env: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Build and validate a complete strict reasoning release."""
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be greater than zero")
    output_root = Path(output_dir)
    modules_dir = output_root / "modules"
    reports_dir = output_root / "reports"
    modules_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    write_extraction_documentation(output_root)
    import_validation = validate_expected_imports(tbox_graph)

    configured_sources = dict(source_files)
    configured_sources["skos"] = skos_source_file

    results: list[ModuleResult] = []
    for spec in VOCABULARIES:
        source_path = configured_sources.get(spec.source_key)
        if not source_path:
            raise RuntimeError(f"No source configured for {spec.name}")
        result = build_reasoning_module(
            spec,
            source_path,
            kg_graph,
            modules_dir / spec.filename,
            kg_version,
        )
        results.append(result)
        selection_report = {
            "module": result.__dict__,
            "policy": extraction_policy(),
        }
        (reports_dir / f"{spec.name}-selection.json").write_text(
            json.dumps(selection_report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(
            f"Reasoning module {spec.name}: used={result.used_terms} "
            f"classes={result.selected_classes} object_properties={result.selected_object_properties} "
            f"datatype_properties={result.selected_datatype_properties} "
            f"datatypes={result.selected_datatypes} triples={result.triples}"
        )
        # RDFLib stores can participate in reference cycles. Force collection
        # after each large pinned ontology has been processed so the next source
        # does not inherit hundreds of megabytes of unreachable graph state.
        gc.collect()

    write_catalog(output_root / "catalog-v001.xml", results)

    if robot_command is None and require_robot:
        raise RuntimeError(
            "ROBOT is required for the strict reasoning release but was not found. "
            "Install `robot` or place robot.jar in project/tools, project/kg, or the project root."
        )

    # ROBOT command tails and the classified TBox are build-time diagnostics.
    # They live in a temporary work directory and are removed after validation.
    status = "BUILT_NOT_VALIDATED"
    tbox_reasoning: dict[str, Any] | None = None
    kg_reasoning: dict[str, Any] | None = None
    assertion_preservation: dict[str, int] | None = None

    merged_tbox = Graph()
    merged_tbox.namespace_manager = tbox_graph.namespace_manager
    copy_without_headers_and_imports(tbox_graph, merged_tbox)
    for result in results:
        module = Graph()
        module.parse(str(output_root / result.path), format="xml")
        copy_without_headers_and_imports(module, merged_tbox)
        module.close()
    tbox_iri = URIRef(f"http://makaao.inria.fr/kg/{kg_version}/curated-tbox")
    merged_tbox.add((tbox_iri, RDF.type, OWL.Ontology))
    merged_tbox.add((tbox_iri, OWL.versionInfo, Literal(kg_version)))
    merged_tbox_path = output_root / f"makaao_kg_{kg_version}_curated-tbox-merged.owl"
    merged_tbox = serialize_verified(merged_tbox, merged_tbox_path)
    merged_tbox_structural = graph_postconditions(
        merged_tbox, require_no_imports=True
    )

    dl_graph = Graph()
    dl_graph.namespace_manager = non_reified_graph.namespace_manager
    copy_without_headers_and_imports(merged_tbox, dl_graph)
    copy_without_headers_and_imports(non_reified_graph, dl_graph)
    dl_graph = remove_generic_reification(dl_graph)
    add_reasoning_metadata_declarations(dl_graph)
    dl_iri = URIRef(f"http://makaao.inria.fr/kg/{kg_version}/curated-reasoning")
    dl_graph.add((dl_iri, RDF.type, OWL.Ontology))
    dl_graph.add((dl_iri, OWL.versionInfo, Literal(kg_version)))
    dl_path = output_root / f"makaao_kg_{kg_version}_curated-full-kg-dl.owl"
    dl_graph = serialize_verified(dl_graph, dl_path)
    dl_graph_structural = graph_postconditions(dl_graph, require_no_imports=True)

    reasoned_kg_path = output_root / f"makaao_kg_{kg_version}_curated-full-kg-dl-reasoned.owl"

    with tempfile.TemporaryDirectory(prefix=".reasoning-work-", dir=output_root) as work:
        work_dir = Path(work)

        if robot_command is not None:
            for result in results:
                validate_profile(
                    robot_command,
                    output_root / result.path,
                    reports_dir / f"{result.name}-profile.txt",
                    work_dir,
                    timeout_seconds,
                    subprocess_env,
                )

            validate_profile(
                robot_command,
                merged_tbox_path,
                reports_dir / "curated-tbox-profile.txt",
                work_dir,
                timeout_seconds,
                subprocess_env,
            )
            classified_tbox_path = work_dir / f"makaao_kg_{kg_version}_curated-tbox-classified.owl"
            run_reasoner(
                robot_command,
                merged_tbox_path,
                classified_tbox_path,
                work_dir,
                timeout_seconds,
                reasoner_name,
                subprocess_env,
            )
            tbox_reasoning = check_reasoned_output(classified_tbox_path)
            if tbox_reasoning["unsatisfiable_named_classes"]:
                raise RuntimeError(
                    "Curated TBox contains unsatisfiable named classes:\n  "
                    + "\n  ".join(tbox_reasoning["unsatisfiable_named_classes"][:30])
                )

            validate_profile(
                robot_command,
                dl_path,
                reports_dir / "curated-full-kg-dl-profile.txt",
                work_dir,
                timeout_seconds,
                subprocess_env,
            )
            run_reasoner(
                robot_command,
                dl_path,
                reasoned_kg_path,
                work_dir,
                timeout_seconds,
                reasoner_name,
                subprocess_env,
            )
            assertion_preservation = preserve_asserted_axioms(
                dl_path,
                reasoned_kg_path,
            )
            kg_reasoning = check_reasoned_output(reasoned_kg_path)
            kg_reasoning["assertion_preservation"] = assertion_preservation
            if kg_reasoning["unsatisfiable_named_classes"]:
                raise RuntimeError(
                    "Reasoned KG contains unsatisfiable named classes:\n  "
                    + "\n  ".join(kg_reasoning["unsatisfiable_named_classes"][:30])
                )
            if kg_reasoning["individuals_typed_owl_nothing"]:
                raise RuntimeError(
                    "Reasoned KG contains individuals typed owl:Nothing:\n  "
                    + "\n  ".join(kg_reasoning["individuals_typed_owl_nothing"][:30])
                )
            validate_profile(
                robot_command,
                reasoned_kg_path,
                reports_dir / "curated-full-kg-dl-reasoned-profile.txt",
                work_dir,
                timeout_seconds,
                subprocess_env,
            )
            status = "PASSED_OWL2_DL_AND_HERMIT"

    manifest = {
        "script_version": SCRIPT_VERSION,
        "script_iteration": SCRIPT_ITERATION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "kg_version": kg_version,
        "status": status,
        "design": {
            "module_set": "single strict module set built directly from pinned external ontology files",
            "anonymous_logical_axioms_retained": False,
            "generic_rdf_reification_retained_in_dl_graph": False,
            "asserted_dl_triples_preserved_in_reasoned_graph": (
                assertion_preservation["asserted_triples_missing_after_restore"] == 0
                if assertion_preservation is not None
                else None
            ),
            "extraction_policy": extraction_policy(),
            "validated_import_set": import_validation,
        },
        "modules": [result.__dict__ for result in results],
        "canonical_inputs": {
            "kg": (
                {
                    "path": f"../{Path(canonical_kg_path).name}",
                    "sha256": file_sha256(canonical_kg_path),
                }
                if canonical_kg_path is not None
                else None
            ),
            "tbox": (
                {
                    "path": f"../{Path(canonical_tbox_path).name}",
                    "sha256": file_sha256(canonical_tbox_path),
                }
                if canonical_tbox_path is not None
                else None
            ),
        },
        "outputs": {
            "merged_tbox": str(merged_tbox_path.name),
            "dl_kg": str(dl_path.name),
            "reasoned_dl_kg": str(reasoned_kg_path.name) if reasoned_kg_path.exists() else None,
        },
        "structural_postconditions": {
            "merged_tbox": merged_tbox_structural,
            "dl_kg": dl_graph_structural,
            "reasoned_dl_kg": (
                kg_reasoning["structural_postconditions"] if kg_reasoning else None
            ),
        },
        "triple_counts": {
            "original_kg": len(kg_graph),
            "non_reified_projection_in_memory": len(non_reified_graph),
            "merged_tbox": len(merged_tbox),
            "dl_kg": len(dl_graph),
            "reasoned_dl_kg": kg_reasoning["triples"] if kg_reasoning else None,
        },
        "tbox_reasoning": tbox_reasoning,
        "kg_reasoning": kg_reasoning,
    }
    manifest_path = output_root / "reasoning-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    readme = output_root / "README.txt"
    readme.write_text(
        "\n".join(
            [
                f"MAKAAO {kg_version} strict reasoning release",
                "=" * 48,
                "",
                f"Status: {status}",
                "",
                "The modules in modules/ are the only external ontology modules produced",
                "by the build. They were extracted directly from the pinned files in",
                "data/external_ontologies; no intermediate modules were used.",
                "",
                "See EXTRACTION_POLICY.txt and reports/*-selection.json for the exact",
                "selection rules and complete retained IRI lists.",
                "",
                f"Canonical reasoning input: {dl_path.name}",
                f"Reasoned output: {reasoned_kg_path.name if reasoned_kg_path.exists() else 'not generated'}",
                f"Merged TBox: {merged_tbox_path.name}",
                "",
                "Only permanent release artifacts are retained. The non-reified projection,",
                "classified TBox, raw curated union, and ROBOT command logs are temporary.",
                "The original provenance/reification graph remains in the normal KG outputs.",
                "The asserted and reasoned DL files are both checked to contain no generic RDF reification.",
                "The retained reasoned DL file contains every asserted DL triple plus reasoner inferences.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    checksum_path = output_root / "SHA256SUMS"
    files = sorted(
        path
        for path in output_root.rglob("*")
        if path.is_file() and path != checksum_path
    )
    with checksum_path.open("w", encoding="utf-8") as stream:
        for path in files:
            stream.write(f"{file_sha256(path)}  {path.relative_to(output_root)}\n")

    return manifest


def main() -> None:
    raise SystemExit(
        "This is an internal helper. Run 03_build_kg_from_tables.py; it stages the "
        "raw KG and the single strict reasoning module/release set transactionally."
    )


if __name__ == "__main__":
    main()