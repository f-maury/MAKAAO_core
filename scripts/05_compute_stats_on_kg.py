#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from rdflib import Graph, RDF, RDFS, OWL, URIRef, Literal
from collections import defaultdict
from pathlib import Path
from typing import Union, Optional, List, Tuple, Dict, Set
import csv
import os
import time
import re
import requests
from gprofiler import GProfiler
import pandas as pd  # only for convenience to write the enrichment file

UMLS_API_KEY = "67bd1b8b-87f0-40de-bf20-255e6f1721a3"


# =====================================================================
# 1. Generic RDF summary
# =====================================================================

def summarize_rdf(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:
    """
    Read RDF graph from `in_path` and write two CSV files into `out_dir`:

      - <prefix>_classes.csv
          class_uri,
          direct_instances,
          instances_with_subclasses,
          direct_subclasses,
          subclasses_with_descendants

      - <prefix>_predicates.csv
          predicate_uri,
          direct_instances,               # used as predicate
          instances_with_subproperties,   # this predicate + all its subproperties
          direct_subproperties,
          subproperties_with_descendants,
          as_object_triples               # used as object (URIRef)

    If `prefix` is None, the stem of `in_path` is used as base name.
    """
    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    # ---------- rdf:type usage: classes and instance types ----------
    classes = {}                    # class URI (str) -> "owl:Class" or "rdfs:Class"
    type_counts = defaultdict(int)  # type URI (str) -> number of subjects

    for s, o in g.subject_objects(RDF.type):
        if o == OWL.Class or o == RDFS.Class:
            classes[str(s)] = "owl:Class" if o == OWL.Class else "rdfs:Class"
        else:
            type_counts[str(o)] += 1

    # ---------- subclass graph (for transitive instance/subclass counts) ----------
    subclass_children = defaultdict(set)  # parent_uri -> set(child_uri)

    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    def instances_with_descendants(root_uri: str) -> int:
        """Sum of direct instances of root_uri and all its subclasses."""
        total = 0
        stack = [root_uri]
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            total += type_counts.get(c, 0)
            for child in subclass_children.get(c, ()):
                if child not in visited:
                    stack.append(child)
        return total

    def subclasses_descendants_count(root_uri: str) -> int:
        """Number of distinct subclasses (descendants) of root_uri."""
        count = 0
        stack = list(subclass_children.get(root_uri, ()))
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            count += 1
            for child in subclass_children.get(c, ()):
                if child not in visited:
                    stack.append(child)
        return count

    # ---------- predicates: counts as predicate and as object ----------
    pred_total = defaultdict(int)    # direct instances: used as predicate
    uri_as_object = defaultdict(int) # used as object (URIRef)

    for s, p, o in g:
        pred_total[str(p)] += 1
        if isinstance(o, URIRef):
            uri_as_object[str(o)] += 1

    # property hierarchy: rdfs:subPropertyOf
    prop_children = defaultdict(set)  # parent_property_uri -> set(child_property_uri)
    for s, p, o in g.triples((None, RDFS.subPropertyOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        prop_children[parent_uri].add(child_uri)

    def predicate_instances_with_descendants(root_uri: str) -> int:
        """Sum of direct instances of root_uri and all its subproperties."""
        total = 0
        stack = [root_uri]
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            total += pred_total.get(c, 0)
            for child in prop_children.get(c, ()):
                if child not in visited:
                    stack.append(child)
        return total

    def subproperties_descendants_count(root_uri: str) -> int:
        """Number of distinct subproperties (descendants) of root_uri."""
        count = 0
        stack = list(prop_children.get(root_uri, ()))
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            count += 1
            for child in prop_children.get(c, ()):
                if child not in visited:
                    stack.append(child)
        return count

    # set of all predicate URIs we care about:
    all_predicates = set(pred_total.keys()) \
        | set(prop_children.keys()) \
        | {c for children in prop_children.values() for c in children}

    # ---------- write <out_dir>/<base>_classes.csv ----------
    class_rows = []
    for c_uri in classes.keys():
        direct_inst = type_counts.get(c_uri, 0)
        with_sub_inst = instances_with_descendants(c_uri)
        direct_sub = len(subclass_children.get(c_uri, set()))
        sub_with_desc = subclasses_descendants_count(c_uri)
        class_rows.append({
            "class_uri": c_uri,
            "direct_instances": direct_inst,
            "instances_with_subclasses": with_sub_inst,
            "direct_subclasses": direct_sub,
            "subclasses_with_descendants": sub_with_desc,
        })
    class_rows.sort(
        key=lambda r: (-r["instances_with_subclasses"],
                       -r["direct_instances"],
                       r["class_uri"])
    )

    classes_path = out_dir / f"{base}_classes.csv"
    with classes_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "class_uri",
                "direct_instances",
                "instances_with_subclasses",
                "direct_subclasses",
                "subclasses_with_descendants",
            ],
        )
        writer.writeheader()
        writer.writerows(class_rows)

    # ---------- write <out_dir>/<base>_predicates.csv ----------
    pred_rows = []
    for p_uri in all_predicates:
        direct_inst = pred_total.get(p_uri, 0)
        with_sub_inst = predicate_instances_with_descendants(p_uri)
        direct_sub = len(prop_children.get(p_uri, set()))
        sub_with_desc = subproperties_descendants_count(p_uri)
        as_obj = uri_as_object.get(p_uri, 0)
        pred_rows.append({
            "predicate_uri": p_uri,
            "direct_instances": direct_inst,
            "instances_with_subproperties": with_sub_inst,
            "direct_subproperties": direct_sub,
            "subproperties_with_descendants": sub_with_desc,
            "as_object_triples": as_obj,
        })
    pred_rows.sort(
        key=lambda r: (-r["instances_with_subproperties"],
                       -r["direct_instances"],
                       r["predicate_uri"])
    )

    preds_path = out_dir / f"{base}_predicates.csv"
    with preds_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "predicate_uri",
                "direct_instances",
                "instances_with_subproperties",
                "direct_subproperties",
                "subproperties_with_descendants",
                "as_object_triples",
            ],
        )
        writer.writeheader()
        writer.writerows(pred_rows)


# =====================================================================
# 2. AAB → Target stats
# =====================================================================

def summarize_aab_targets(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:

    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    AAB_PREFIX = MAK + "aab_"
    TARGET_CLASS = MAK + "Target"
    AUTOANTIBODY_CLASS = MAK + "Autoantibody"  # <- add this

    # ---------- 1) Collect classes and find AAB classes ----------
    class_kinds = {}
    aab_classes = set()
    instance_types = defaultdict(set)

    for s, o in g.subject_objects(RDF.type):
        s_uri = str(s)
        o_uri = str(o)

        instance_types[s_uri].add(o_uri)

        if o in (OWL.Class, RDFS.Class):
            class_kinds[s_uri] = o
            if s_uri.startswith(AAB_PREFIX):
                aab_classes.add(s_uri)

    # Explicitly include Autoantibody as an AAB class/root if present
    if AUTOANTIBODY_CLASS in class_kinds:
        aab_classes.add(AUTOANTIBODY_CLASS)


    # ---------- 2) Build subclass graph (all classes) ----------
    subclass_children = defaultdict(set)  # parent_uri -> set(child_uri)

    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- 3) AAB instances and Target instances ----------
    aab_instances_by_class = defaultdict(set)  # AAB class URI -> set(instances)
    target_instances = set()                   # set of instance URIs

    for inst_uri, types in instance_types.items():
        # mark Target instances
        if TARGET_CLASS in types:
            target_instances.add(inst_uri)

        # mark AAB instances (can be multiple AAB types, just in case)
        for cls in types:
            if cls in aab_classes:
                aab_instances_by_class[cls].add(inst_uri)

    # ---------- 4) Map AAB class -> Targets linked to its instances ----------
    targets_by_aab = {cls: set() for cls in aab_classes}

    for s, p, o in g:
        s_uri = str(s)
        o_uri = str(o)

        # case 1: s is AAB instance, o is Target instance
        if o_uri in target_instances:
            for cls in instance_types.get(s_uri, ()):
                if cls in aab_classes:
                    targets_by_aab[cls].add(o_uri)

        # case 2: o is AAB instance, s is Target instance
        if s_uri in target_instances:
            for cls in instance_types.get(o_uri, ()):
                if cls in aab_classes:
                    targets_by_aab[cls].add(s_uri)

    # ---------- 5) Targets including subclasses ----------
    def targets_with_descendants(aab_class_uri: str) -> set:
        """
        Return the set of distinct Target instances linked to instances of
        aab_class_uri or any of its AAB subclasses (via rdfs:subClassOf*).
        """
        result = set()
        stack = [aab_class_uri]
        visited = set()

        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)

            # add direct targets of this class
            result.update(targets_by_aab.get(c, set()))

            # go down to subclasses, but only within AAB classes
            for child in subclass_children.get(c, ()):
                if child in aab_classes and child not in visited:
                    stack.append(child)

        return result

    # ---------- 6) Build rows and write CSV ----------
    rows = []
    for cls in aab_classes:
        direct_targets = len(targets_by_aab.get(cls, set()))
        all_targets = len(targets_with_descendants(cls))
        rows.append({
            "aab_class_uri": cls,
            "direct_target_instances": direct_targets,
            "target_instances_with_subclasses": all_targets,
        })

    rows.sort(
        key=lambda r: (-r["target_instances_with_subclasses"],
                       -r["direct_target_instances"],
                       r["aab_class_uri"])
    )

    out_path = out_dir / f"{base}_aab_targets.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "aab_class_uri",
                "direct_target_instances",
                "target_instances_with_subclasses",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


# =====================================================================
# 2b. AAB leaf stats (terminal AAB subclasses)
# =====================================================================

def summarize_aab_leaf_stats(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:
    """
    For each AAB class, count how many AAB leaf subclasses it has.

    A leaf AAB class = AAB class that has no AAB subclasses
                       (i.e., no child in subclass_children that is an AAB class).

    Output CSV: <prefix>_aab_leaf_stats.csv

    Columns:
      aab_class_uri
      aab_label                # main English name (rdfs:label@en if possible)
      n_leaf_aab_subclasses
      leaf_aab_classes         # pipe-separated URIs of leaf AABs in its subtree
    """
    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    AAB_PREFIX = MAK + "aab_"
    AUTOANTIBODY_CLASS = MAK + "Autoantibody"

    # ---------- helper: get main English label for a URI ----------
    def get_en_label(uri_str: str) -> str:
        """
        Prefer rdfs:label with @en language tag.
        If none, take the first rdfs:label (any language).
        If no label, return empty string.
        """
        u = URIRef(uri_str)
        best_any = None
        for lbl in g.objects(u, RDFS.label):
            if isinstance(lbl, Literal):
                if lbl.language == "en":
                    return str(lbl)
                if best_any is None:
                    best_any = str(lbl)
        return best_any or ""

    # ---------- 1) Collect classes and AAB classes ----------
    classes: Set[str] = set()
    aab_classes: Set[str] = set()

    for s, o in g.subject_objects(RDF.type):
        if o in (OWL.Class, RDFS.Class):
            s_uri = str(s)
            classes.add(s_uri)
            if s_uri.startswith(AAB_PREFIX):
                aab_classes.add(s_uri)

    # Explicitly include Autoantibody as an AAB root class if present
    if AUTOANTIBODY_CLASS in classes:
        aab_classes.add(AUTOANTIBODY_CLASS)

    if not aab_classes:
        print("[AAB LEAVES] No AAB classes found; nothing to do.")
        return

    # ---------- 2) Subclass graph over all classes ----------
    subclass_children: Dict[str, Set[str]] = defaultdict(set)
    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- 3) Helper: compute leaf AABs under an AAB class ----------
    def aab_leaves_for_class(root_uri: str) -> Set[str]:
        # 1) Traverse ALL subclasses, not only AAB-labeled ones
        subtree: Set[str] = set()
        stack = [root_uri]
        visited: Set[str] = set()

        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)

            # Keep only AAB classes in the subtree set
            if c in aab_classes:
                subtree.add(c)

            # But always traverse through any child class, AAB or not
            for child in subclass_children.get(c, ()):
                if child not in visited:
                    stack.append(child)

        # 2) Leaves = nodes in subtree without any AAB child in subtree
        leaves: Set[str] = set()
        for node in subtree:
            has_aab_child = any(
                (child in subtree)
                for child in subclass_children.get(node, ())
            )
            if not has_aab_child:
                leaves.add(node)

        return leaves

    # ---------- 4) Build rows and write CSV ----------
    rows = []
    for cls in aab_classes:
        leaves = aab_leaves_for_class(cls)
        rows.append({
            "aab_class_uri": cls,
            "aab_label": get_en_label(cls),   # NEW COLUMN
            "n_leaf_aab_subclasses": len(leaves),
            "leaf_aab_classes": "|".join(sorted(leaves)),
        })

    rows.sort(
        key=lambda r: (-r["n_leaf_aab_subclasses"], r["aab_class_uri"])
    )

    out_path = out_dir / f"{base}_aab_leaf_stats.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "aab_class_uri",
                "aab_label",              # 2nd position
                "n_leaf_aab_subclasses",
                "leaf_aab_classes",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[AAB LEAVES] Wrote {out_path}")





# =====================================================================
# 3. AAB ↔ Disease class stats
# =====================================================================

def summarize_aab_diseases(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:
    
    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    AAB_PREFIX = MAK + "aab_"
    DISEASE_ROOT = MAK + "AutoimmuneDisease"

    # ---------- 1) Collect instance types and AAB classes ----------
    instance_types = defaultdict(set)  # inst_uri -> set of class URIs
    aab_classes = set()
    classes = set()

    for s, o in g.subject_objects(RDF.type):
        s_uri = str(s)
        o_uri = str(o)
        instance_types[s_uri].add(o_uri)

        # record classes
        if o in (OWL.Class, RDFS.Class):
            classes.add(s_uri)
            if s_uri.startswith(AAB_PREFIX):
                aab_classes.add(s_uri)

    # ---------- 2) Subclass graph over all classes ----------
    subclass_children = defaultdict(set)  # parent_uri -> set(child_uri)
    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- 3) Disease classes: AutoimmuneDisease and its subclasses ----------
    disease_classes = set()
    if DISEASE_ROOT in classes:
        stack = [DISEASE_ROOT]
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            disease_classes.add(c)
            for child in subclass_children.get(c, ()):
                stack.append(child)

    # ---------- 4) Disease instances and AAB instances ----------
    disease_instances = set()
    aab_instances_by_class = defaultdict(set)

    for inst_uri, types in instance_types.items():
        # disease if any type is in disease_classes
        if disease_classes and any(t in disease_classes for t in types):
            disease_instances.add(inst_uri)

        # AAB instances: inst with types in aab_classes
        for cls in types:
            if cls in aab_classes:
                aab_instances_by_class[cls].add(inst_uri)

    # ---------- 5) Map AAB class -> diseases linked to its instances ----------
    diseases_by_aab = {cls: set() for cls in aab_classes}

    for s, p, o in g:
        s_uri = str(s)
        o_uri = str(o)

        # case 1: s is AAB instance, o is disease instance
        if o_uri in disease_instances:
            for cls in instance_types.get(s_uri, ()):
                if cls in aab_classes:
                    diseases_by_aab[cls].add(o_uri)

        # case 2: o is AAB instance, s is disease instance
        if s_uri in disease_instances:
            for cls in instance_types.get(o_uri, ()):
                if cls in aab_classes:
                    diseases_by_aab[cls].add(s_uri)

    # ---------- 6) Disease counts including AAB subclasses ----------
    def diseases_with_descendants(aab_class_uri: str) -> set:
        """
        Return distinct disease instances linked to instances of
        aab_class_uri or any of its AAB subclasses (via rdfs:subClassOf*).
        """
        result = set()
        stack = [aab_class_uri]
        visited = set()

        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)

            result.update(diseases_by_aab.get(c, set()))

            # restrict to AAB subclasses
            for child in subclass_children.get(c, ()):
                if child in aab_classes and child not in visited:
                    stack.append(child)

        return result

    # ---------- 7) Build rows and write CSV ----------
    rows = []
    for cls in aab_classes:
        direct_diseases = len(diseases_by_aab.get(cls, set()))
        all_diseases = len(diseases_with_descendants(cls))
        rows.append({
            "aab_class_uri": cls,
            "direct_disease_instances": direct_diseases,
            "disease_instances_with_subclasses": all_diseases,
        })

    rows.sort(
        key=lambda r: (-r["disease_instances_with_subclasses"],
                       -r["direct_disease_instances"],
                       r["aab_class_uri"])
    )

    out_path = out_dir / f"{base}_aab_diseases.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "aab_class_uri",
                "direct_disease_instances",
                "disease_instances_with_subclasses",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def summarize_disease_aab(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:
    
    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    AAB_PREFIX = MAK + "aab_"
    DISEASE_ROOT = MAK + "AutoimmuneDisease"

    # ---------- 1) Collect instance types, classes, AAB classes ----------
    instance_types = defaultdict(set)  # inst_uri -> set of class URIs
    classes = set()
    aab_classes = set()

    for s, o in g.subject_objects(RDF.type):
        s_uri = str(s)
        o_uri = str(o)
        instance_types[s_uri].add(o_uri)

        if o in (OWL.Class, RDFS.Class):
            classes.add(s_uri)
            if s_uri.startswith(AAB_PREFIX):
                aab_classes.add(s_uri)

    # ---------- 2) Subclass graph for all classes ----------
    subclass_children = defaultdict(set)  # parent_uri -> set(child_uri)
    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- 3) Disease classes: AutoimmuneDisease and its subclasses ----------
    disease_classes = set()
    if DISEASE_ROOT in classes:
        stack = [DISEASE_ROOT]
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            disease_classes.add(c)
            for child in subclass_children.get(c, ()):
                stack.append(child)

    # ---------- 4) Disease instances and AAB instances ----------
    disease_instances_by_class = defaultdict(set)  # disease_class -> set(instances)
    aab_instances = set()

    for inst_uri, types in instance_types.items():
        # AAB instance?
        if any(t in aab_classes for t in types):
            aab_instances.add(inst_uri)

        # disease instance: collect per disease class
        for cls in types:
            if cls in disease_classes:
                disease_instances_by_class[cls].add(inst_uri)

    # ---------- 5) Map disease class -> AABs linked to its instances ----------
    aab_by_disease = {cls: set() for cls in disease_classes}

    for s, p, o in g:
        s_uri = str(s)
        o_uri = str(o)

        # case 1: s is disease instance, o is AAB instance
        if o_uri in aab_instances:
            # all disease classes of s_uri
            for cls in instance_types.get(s_uri, ()):
                if cls in disease_classes:
                    aab_by_disease[cls].add(o_uri)

        # case 2: o is disease instance, s is AAB instance
        if s_uri in aab_instances:
            for cls in instance_types.get(o_uri, ()):
                if cls in disease_classes:
                    aab_by_disease[cls].add(s_uri)

    # ---------- 6) AAB counts including disease subclasses ----------
    def aab_with_descendants(disease_class_uri: str) -> set:
        """
        Return distinct AAB instances linked to instances of disease_class_uri
        or any of its disease subclasses (via rdfs:subClassOf*).
        """
        result = set()
        stack = [disease_class_uri]
        visited = set()

        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)

            result.update(aab_by_disease.get(c, set()))

            # restrict traversal to disease subclasses
            for child in subclass_children.get(c, ()):
                if child in disease_classes and child not in visited:
                    stack.append(child)

        return result

    # ---------- 7) Build rows and write CSV ----------
    rows = []
    for cls in disease_classes:
        direct_aab = len(aab_by_disease.get(cls, set()))
        all_aab = len(aab_with_descendants(cls))
        rows.append({
            "disease_class_uri": cls,
            "direct_aab_instances": direct_aab,
            "aab_instances_with_subclasses": all_aab,
        })

    rows.sort(
        key=lambda r: (-r["aab_instances_with_subclasses"],
                       -r["direct_aab_instances"],
                       r["disease_class_uri"])
    )

    out_path = out_dir / f"{base}_disease_aab.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "disease_class_uri",
                "direct_aab_instances",
                "aab_instances_with_subclasses",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def summarize_disease_instances_aab(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:

    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    AAB_PREFIX = MAK + "aab_"
    DISEASE_ROOT = MAK + "AutoimmuneDisease"

    # ---------- 1) Collect instance types and classes ----------
    instance_types = defaultdict(set)  # inst_uri -> set of class URIs
    classes = set()
    aab_classes = set()

    for s, o in g.subject_objects(RDF.type):
        s_uri = str(s)
        o_uri = str(o)
        instance_types[s_uri].add(o_uri)

        if o in (OWL.Class, RDFS.Class):
            classes.add(s_uri)
            if s_uri.startswith(AAB_PREFIX):
                aab_classes.add(s_uri)

    # ---------- 2) Subclass graph over all classes ----------
    subclass_children = defaultdict(set)  # parent_uri -> set(child_uri)
    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- 3) Disease classes: AutoimmuneDisease and its subclasses ----------
    disease_classes = set()
    if DISEASE_ROOT in classes:
        stack = [DISEASE_ROOT]
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            disease_classes.add(c)
            for child in subclass_children.get(c, ()):
                stack.append(child)

    # ---------- 4) Determine disease instances and AAB instances ----------
    disease_instances = set()
    aab_instances = set()

    for inst_uri, types in instance_types.items():
        # AAB instance?
        if any(t in aab_classes for t in types):
            aab_instances.add(inst_uri)

        # disease instance?
        if any(t in disease_classes for t in types):
            disease_instances.add(inst_uri)

    # ---------- 5) For each disease instance, collect linked AAB instances ----------
    aab_by_disease_instance = {d: set() for d in disease_instances}

    for s, p, o in g:
        s_uri = str(s)
        o_uri = str(o)

        # case 1: s is disease instance, o is AAB instance
        if s_uri in disease_instances and o_uri in aab_instances:
            aab_by_disease_instance[s_uri].add(o_uri)

        # case 2: o is disease instance, s is AAB instance
        if o_uri in disease_instances and s_uri in aab_instances:
            aab_by_disease_instance[o_uri].add(s_uri)

    # ---------- 6) Build rows and write CSV ----------
    rows = []
    for d in disease_instances:
        direct_aab = len(aab_by_disease_instance.get(d, set()))
        rows.append({
            "disease_instance_uri": d,
            "direct_aab_instances": direct_aab,
        })

    rows.sort(
        key=lambda r: (-r["direct_aab_instances"], r["disease_instance_uri"])
    )

    out_path = out_dir / f"{base}_disease_instances_aab.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["disease_instance_uri", "direct_aab_instances"],
        )
        writer.writeheader()
        writer.writerows(rows)


# =====================================================================
# 4. AAB → LOINC parts
# =====================================================================

def summarize_aab_loinc_parts(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:
    
    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    AAB_PREFIX = MAK + "aab_"

    def is_loinc_part(uri: str) -> bool:
        return uri.startswith("https://loinc.org/LP") or uri.startswith("http://loinc.org/LP")

    # ---------- 1) Collect instance types, classes, and AAB classes ----------
    instance_types = defaultdict(set)  # inst_uri -> set of class URIs
    classes = set()
    aab_classes = set()

    for s, o in g.subject_objects(RDF.type):
        s_uri = str(s)
        o_uri = str(o)
        instance_types[s_uri].add(o_uri)

        if o in (OWL.Class, RDFS.Class):
            classes.add(s_uri)
            if s_uri.startswith(AAB_PREFIX):
                aab_classes.add(s_uri)

    # ---------- 2) Subclass graph for all classes ----------
    subclass_children = defaultdict(set)  # parent_uri -> set(child_uri)
    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- 3) AAB instances per AAB class ----------
    aab_instances_by_class = defaultdict(set)  # AAB class -> set(instances)

    for inst_uri, types in instance_types.items():
        for cls in types:
            if cls in aab_classes:
                aab_instances_by_class[cls].add(inst_uri)

    # ---------- 4) Map AAB class -> LOINC parts linked to it ----------
    loinc_by_aab_class = {cls: set() for cls in aab_classes}

    for s, p, o in g:
        s_uri = str(s)
        o_uri = str(o)

        # case 1: subject AAB / AAB instance, object LOINC part
        if is_loinc_part(o_uri):
            # s is an AAB class directly
            if s_uri in aab_classes:
                loinc_by_aab_class[s_uri].add(o_uri)
            # s is an AAB instance
            for cls in instance_types.get(s_uri, ()):
                if cls in aab_classes:
                    loinc_by_aab_class[cls].add(o_uri)

        # case 2: subject LOINC part, object AAB / AAB instance
        if is_loinc_part(s_uri):
            # o is an AAB class directly
            if o_uri in aab_classes:
                loinc_by_aab_class[o_uri].add(s_uri)
            # o is an AAB instance
            for cls in instance_types.get(o_uri, ()):
                if cls in aab_classes:
                    loinc_by_aab_class[cls].add(s_uri)

    # ---------- 5) LOINC counts including AAB subclasses ----------
    def loinc_with_descendants(aab_class_uri: str) -> set:
        """
        Return distinct LOINC parts linked to instances of aab_class_uri
        or any of its AAB subclasses (via rdfs:subClassOf*).
        """
        result = set()
        stack = [aab_class_uri]
        visited = set()

        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)

            result.update(loinc_by_aab_class.get(c, set()))

            # restrict traversal to AAB subclasses
            for child in subclass_children.get(c, ()):
                if child in aab_classes and child not in visited:
                    stack.append(child)

        return result

    # ---------- 6) Build rows and write CSV ----------
    rows = []
    for cls in aab_classes:
        direct_loinc = len(loinc_by_aab_class.get(cls, set()))
        all_loinc = len(loinc_with_descendants(cls))
        rows.append({
            "aab_class_uri": cls,
            "direct_loinc_parts": direct_loinc,
            "loinc_parts_with_subclasses": all_loinc,
        })

    rows.sort(
        key=lambda r: (-r["loinc_parts_with_subclasses"],
                       -r["direct_loinc_parts"],
                       r["aab_class_uri"])
    )

    out_path = out_dir / f"{base}_aab_loinc_parts.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames[
                "aab_class_uri",
                "direct_loinc_parts",
                "loinc_parts_with_subclasses",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


# =====================================================================
# 5. Orpha → phenotypes
# =====================================================================

def summarize_orpha_phenotypes(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:
    
    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    ORPHA_PREFIX = MAK + "orpha_"
    HP_INSTANCE_PREFIX = MAK + "HP_"
    HPO_CLASS_PREFIX = "http://purl.obolibrary.org/obo/HP_"

    def is_orpha_instance(uri: str) -> bool:
        return uri.startswith(ORPHA_PREFIX)

    def is_phenotype(uri: str) -> bool:
        # match your HP_*_instance URIs, and optionally raw HPO class URIs
        return uri.startswith(HP_INSTANCE_PREFIX) or uri.startswith(HPO_CLASS_PREFIX)

    # ---------- 1) Collect all Orpha disease instances (subject or object) ----------
    orpha_instances = set()

    for s in g.subjects():
        s_uri = str(s)
        if is_orpha_instance(s_uri):
            orpha_instances.add(s_uri)

    for o in g.objects():
        o_uri = str(o)
        if isinstance(o, URIRef) and is_orpha_instance(o_uri):
            orpha_instances.add(o_uri)

    # ---------- 2) For each Orpha instance, collect linked phenotype instances ----------
    phenos_by_orpha = {o: set() for o in orpha_instances}

    for s, p, o in g:
        s_uri = str(s)
        o_uri = str(o)

        # case 1: Orpha -> phenotype
        if s_uri in orpha_instances and isinstance(o, URIRef) and is_phenotype(o_uri):
            phenos_by_orpha[s_uri].add(o_uri)

        # case 2: phenotype -> Orpha
        if isinstance(s, URIRef) and is_phenotype(s_uri) and o_uri in orpha_instances:
            phenos_by_orpha[o_uri].add(s_uri)

    # ---------- 3) Build rows and write CSV ----------
    rows = []
    for d in orpha_instances:
        count = len(phenos_by_orpha.get(d, set()))
        rows.append({
            "disease_instance_uri": d,
            "direct_phenotype_instances": count,
        })

    rows.sort(
        key=lambda r: (-r["direct_phenotype_instances"], r["disease_instance_uri"])
    )

    out_path = out_dir / f"{base}_orpha_phenotypes.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["disease_instance_uri", "direct_phenotype_instances"],
        )
        writer.writeheader()
        writer.writerows(rows)


# =====================================================================
# 6. Target origin summary
# =====================================================================

def summarize_target_origins(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
) -> None:
    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    TARGET_CLASS = MAK + "Target"

    # ---------- 1) collect all classes and subclass graph ----------
    classes = set()
    subclass_children = defaultdict(set)  # parent_uri -> set(child_uri)

    for s, o in g.subject_objects(RDF.type):
        if o in (OWL.Class, RDFS.Class):
            classes.add(str(s))

    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- 2) compute Target and all its subclasses ----------
    target_classes = set()
    if TARGET_CLASS in classes:
        stack = [TARGET_CLASS]
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            target_classes.add(c)
            for child in subclass_children.get(c, ()):
                stack.append(child)

    # ---------- 3) collect instances of Target or any subclass ----------
    target_instances = set()

    for s, o in g.subject_objects(RDF.type):
        s_uri = str(s)
        o_uri = str(o)
        if o_uri in target_classes:
            target_instances.add(s_uri)

    # ---------- 4) classify origin ----------
    def classify_origin(uri: str) -> str:
        if uri.startswith(MAK + "UP_"):
            return "UniProt"
        if uri.startswith(MAK + "CHEBI_"):
            return "ChEBI"
        if uri.startswith("https://uts.nlm.nih.gov/uts/umls/concept/"):
            return "UMLS"
        return "Other"

    counts = defaultdict(int)
    for inst in target_instances:
        origin = classify_origin(inst)
        counts[origin] += 1

    # ---------- 5) write summary CSV ----------
    summary_path = out_dir / f"{base}_target_origins_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["origin_type", "n_targets"],
        )
        writer.writeheader()
        for origin in sorted(counts.keys()):
            writer.writerow({"origin_type": origin, "n_targets": counts[origin]})


# =====================================================================
# 7. UniProt annotations for Target instances
# =====================================================================

def fetch_uniprot_annotations_for_targets(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    prefix: Optional[str] = None,
    email: Optional[str] = None,
    sleep_sec: float = 0.2,
) -> None:
    """
    1. Parse RDF graph, find Target instances and associated UniProt accessions.
    2. For each accession, fetch UniProt annotations (GO, keywords, pathways, tissue).
    3. Write a per-target annotation CSV: <prefix>_uniprot_annotations.csv
    4. Run g:Profiler GO enrichment on all unique UniProt accessions and
       write a GO p-value file: <prefix>_go_enrichment_gprofiler.csv
    """
    g = Graph()
    g.parse(in_path)

    MAK = "http://makaao.inria.fr/kg/"
    TARGET_CLASS = MAK + "Target"

    # ---------- collect classes and subclass graph ----------
    classes = set()
    subclass_children = defaultdict(set)

    for s, o in g.subject_objects(RDF.type):
        if o in (OWL.Class, RDFS.Class):
            classes.add(str(s))

    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- Target and all subclasses ----------
    target_classes = set()
    if TARGET_CLASS in classes:
        stack = [TARGET_CLASS]
        visited = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            target_classes.add(c)
            for child in subclass_children.get(c, ()):
                stack.append(child)

    # ---------- collect Target instances ----------
    target_instances: List[str] = []
    for s, o in g.subject_objects(RDF.type):
        s_uri = str(s)
        o_uri = str(o)
        if o_uri in target_classes:
            target_instances.append(s_uri)

    # ---------- extract UniProt accessions ----------
    def extract_uniprot_accession(uri: str) -> Optional[str]:
        marker = "UP_"
        if marker not in uri:
            return None
        after = uri.split(marker, 1)[1]  # e.g. "P01308_instance"
        if "_instance" in after:
            acc = after.split("_instance", 1)[0]
        else:
            acc = after
        acc = acc.strip()
        return acc or None

    uniprot_targets: List[Tuple[str, str]] = []
    for t_uri in target_instances:
        acc = extract_uniprot_accession(t_uri)
        if acc is not None:
            uniprot_targets.append((t_uri, acc))

    # ---------- UniProt REST client ----------
    base_url = "https://rest.uniprot.org/uniprotkb/"
    headers: Dict[str, str] = {}
    if email:
        headers["User-Agent"] = f"makg-uniprot-client/1.0 ({email})"

    def fetch_uniprot_annotations(
        accession: str,
    ) -> Tuple[List[str], List[str], List[str], List[str], List[str],
               List[str], List[str], List[str], List[str]]:
        """
        Return:
          go_ids, go_terms, aspect_codes, aspect_names, evidence_types,
          keyword_names, keyword_categories,
          pathways, tissue_specificity
        """
        url = base_url + accession
        params = {"format": "json"}
        resp = requests.get(url, params=params, headers=headers, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        # ---- GO ----
        go_ids: List[str] = []
        go_terms: List[str] = []
        aspect_codes: List[str] = []
        aspect_names: List[str] = []
        evidence_types: List[str] = []

        for ann in data.get("uniProtKBCrossReferences", []) or []:
            if ann.get("database") != "GO":
                continue

            go_id = ann.get("id", "")
            if not go_id:
                continue

            term = ""
            aspect_code = ""
            aspect_name = ""
            evidence = ""

            for prop in ann.get("properties", []) or []:
                key = prop.get("key", "")
                val = prop.get("value", "") or ""
                if key == "GoTerm":
                    term = val
                elif key == "term":
                    parts = val.split(":", 1)
                    aspect_code = parts[0].strip() if parts else ""
                    if len(parts) > 1:
                        aspect_name = parts[1].strip()
                elif key == "GoEvidenceType":
                    evidence = val

            go_ids.append(go_id)
            go_terms.append(term)
            aspect_codes.append(aspect_code)
            aspect_names.append(aspect_name)
            evidence_types.append(evidence)

        # ---- keywords ----
        keyword_names: List[str] = []
        keyword_categories: List[str] = []
        for kw in data.get("keywords", []) or []:
            name = kw.get("name", "") or ""
            cat = kw.get("category", "") or ""
            if name:
                keyword_names.append(name)
                keyword_categories.append(cat)

        # ---- pathways ----
        pathways: List[str] = []
        pathway_dbs = {"Reactome", "KEGG", "BioCyc"}
        for xref in data.get("uniProtKBCrossReferences", []) or []:
            db = xref.get("database", "")
            if db not in pathway_dbs:
                continue
            pid = xref.get("id", "")
            if pid:
                pathways.append(f"{db}:{pid}")

        # ---- tissue specificity ----
        tissues: List[str] = []
        for c in data.get("comments", []) or []:
            if c.get("commentType") == "TISSUE SPECIFICITY":
                for txt in c.get("texts", []) or []:
                    val = txt.get("value", "") or ""
                    if val:
                        tissues.append(val)

        return (
            go_ids,
            go_terms,
            aspect_codes,
            aspect_names,
            evidence_types,
            keyword_names,
            keyword_categories,
            pathways,
            tissues,
        )

    # ---------- fetch annotations for each UniProt target ----------
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = prefix if prefix is not None else Path(in_path).stem

    rows = []
    all_uniprot_accessions: set[str] = set()

    print(f"[UniProt] Found {len(uniprot_targets)} UniProt-based targets")
    for i, (target_uri, acc) in enumerate(sorted(uniprot_targets, key=lambda x: x[1]), start=1):
        all_uniprot_accessions.add(acc)

        # per-category containers
        mf_ids: List[str] = []
        mf_terms: List[str] = []
        cc_ids: List[str] = []
        cc_terms: List[str] = []
        bp_ids: List[str] = []
        bp_terms: List[str] = []

        try:
            (
                go_ids,
                go_terms,
                aspect_codes,
                aspect_names,
                evidence_types,
                keyword_names,
                keyword_categories,
                pathways,
                tissues,
            ) = fetch_uniprot_annotations(acc)

            # split GO annotations by aspect:
            # F = molecular function, C = cellular component, P = biological process
            for gid, term, acode in zip(go_ids, go_terms, aspect_codes):
                if acode == "F":
                    mf_ids.append(gid)
                    mf_terms.append(term)
                elif acode == "C":
                    cc_ids.append(gid)
                    cc_terms.append(term)
                elif acode == "P":
                    bp_ids.append(gid)
                    bp_terms.append(term)

            api_error = ""
        except Exception as e:
            go_ids = []
            go_terms = []
            aspect_codes = []
            aspect_names = []
            evidence_types = []
            keyword_names = []
            keyword_categories = []
            pathways = []
            tissues = []
            api_error = str(e)

        rows.append({
            "target_instance_uri": target_uri,
            "uniprot_accession": acc,
            "n_go_terms": len(go_ids),
            "n_mf_go_terms": len(mf_ids),
            "n_cc_go_terms": len(cc_ids),
            "n_bp_go_terms": len(bp_ids),
            "go_ids": "|".join(go_ids),
            "go_terms": "|".join(go_terms),
            "go_aspect_codes": "|".join(aspect_codes),
            "go_aspect_names": "|".join(aspect_names),
            "go_evidence_types": "|".join(evidence_types),
            "mf_go_ids": "|".join(mf_ids),
            "mf_go_terms": "|".join(mf_terms),
            "cc_go_ids": "|".join(cc_ids),
            "cc_go_terms": "|".join(cc_terms),
            "bp_go_ids": "|".join(bp_ids),
            "bp_go_terms": "|".join(bp_terms),
            "keywords": "|".join(keyword_names),
            "keyword_categories": "|".join(keyword_categories),
            "pathways": "|".join(pathways),
            "tissue_specificity": "|".join(tissues),
            "api_error": api_error,
        })

        if i % 20 == 0:
            print(f"[UniProt] Processed {i}/{len(uniprot_targets)} targets ...")

        time.sleep(sleep_sec)

    # ---------- write per-target annotation CSV ("count file") ----------
    out_path = out_dir / f"{base}_uniprot_annotations.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "target_instance_uri",
                "uniprot_accession",
                "n_go_terms",
                "n_mf_go_terms",
                "n_cc_go_terms",
                "n_bp_go_terms",
                "go_ids",
                "go_terms",
                "go_aspect_codes",
                "go_aspect_names",
                "go_evidence_types",
                "mf_go_ids",
                "mf_go_terms",
                "cc_go_ids",
                "cc_go_terms",
                "bp_go_ids",
                "bp_go_terms",
                "keywords",
                "keyword_categories",
                "pathways",
                "tissue_specificity",
                "api_error",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[UniProt] Wrote {out_path}")

    # ---------- g:Profiler GO enrichment (p-value file) ----------
    all_uniprot_list = sorted(all_uniprot_accessions)
    if not all_uniprot_list:
        print("[g:Profiler] No UniProt accessions found; skipping enrichment.")
        return

    print(f"[g:Profiler] Running GO enrichment on {len(all_uniprot_list)} unique UniProt accessions ...")

    gp = GProfiler(return_dataframe=True)
    enrichment_df = gp.profile(
        organism="hsapiens",
        query=all_uniprot_list,
        sources=["GO:BP", "GO:MF", "GO:CC"],  # only GO terms
        user_threshold=1.0,                   # keep all, you can threshold later on p_value
    )

    enr_path = out_dir / f"{base}_go_enrichment_gprofiler.csv"
    enrichment_df.to_csv(enr_path, index=False)

    print(f"[g:Profiler] Wrote {enr_path}")


# =====================================================================
# 8. MONDO annotations for CUIs (helper)
# =====================================================================

def build_mondo_annotations_for_cuis(
    mondo_owl_path: Union[str, Path],
    cuis_of_interest: Set[str],
) -> Dict[str, Dict[str, Set[str]]]:
    """
    Load MONDO OWL file and, for the given set of CUIs, return a mapping:

      CUI -> {
          "mondo_ids":        { "MONDO:0001234", ... },
          "mondo_labels":     { "Some disease", ... },
          "locations":        { "liver", "thyroid gland", ... },          # RO:0004026
          "material_basis_in":{ "autoimmunity", "genetic variant", ... }, # RO:0004028
          "features":         { "myositis", "encephalitis", ... },        # RO:0004025
          "dysfunction_of":   { "T-cell", "complement system", ... },     # RO:0004007
          "modifiers":        { "autoimmune", "inherited", "viral", ... } # RO:0002200
      }

    Labels are taken from rdfs:label when possible, otherwise a local ID is used.
    """

    mondo_owl_path = Path(mondo_owl_path)
    mondo_g = Graph()
    print(f"[MONDO] Parsing MONDO from {mondo_owl_path} ...")
    mondo_g.parse(mondo_owl_path)
    print(f"[MONDO] Parsed MONDO with {len(mondo_g)} triples")

    # --- Predicates / namespaces we care about ---
    HAS_DBXREF = URIRef("http://www.geneontology.org/formats/oboInOwl#hasDbXref")

    RO_HAS_LOCATION            = URIRef("http://purl.obolibrary.org/obo/RO_0004026")
    RO_HAS_MATERIAL_BASIS_IN   = URIRef("http://purl.obolibrary.org/obo/RO_0004028")
    RO_HAS_FEATURE             = URIRef("http://purl.obolibrary.org/obo/RO_0004025")
    RO_DISEASE_HAS_DYSFUNCTION = URIRef("http://purl.obolibrary.org/obo/RO_0004007")
    RO_HAS_MODIFIER            = URIRef("http://purl.obolibrary.org/obo/RO_0002200")

    MONDO_PREFIX = "http://purl.obolibrary.org/obo/MONDO_"

    # --- helper: get a readable label for a URI ---
    def get_label(u: URIRef) -> str:
        for lbl in mondo_g.objects(u, RDFS.label):
            return str(lbl)
        s = str(u)
        if "#" in s:
            return s.rsplit("#", 1)[1]
        return s.rsplit("/", 1)[-1]

    # --- helper: collect fillers from an OWL class expression / restriction ---
    def collect_restriction_fillers(expr, buckets: Dict[str, Set[str]]) -> None:
        """
        expr: BNode or class expression reachable from a MONDO class.

        Look for owl:Restriction with owl:onProperty = RO_xxx and
        owl:someValuesFrom / owl:allValuesFrom = filler URI.
        Handles both simple restrictions and those inside owl:intersectionOf lists.
        """
        from rdflib.term import BNode

        if not isinstance(expr, BNode):
            return

        # Direct restriction: onProperty + some/allValuesFrom
        on_prop = mondo_g.value(expr, OWL.onProperty)
        filler = (
            mondo_g.value(expr, OWL.someValuesFrom)
            or mondo_g.value(expr, OWL.allValuesFrom)
        )

        if isinstance(on_prop, URIRef) and isinstance(filler, URIRef):
            filler_label = get_label(filler)

            if on_prop == RO_HAS_LOCATION:
                buckets["locations"].add(filler_label)
            elif on_prop == RO_HAS_MATERIAL_BASIS_IN:
                buckets["material_basis_in"].add(filler_label)
            elif on_prop == RO_HAS_FEATURE:
                buckets["features"].add(filler_label)
            elif on_prop == RO_DISEASE_HAS_DYSFUNCTION:
                buckets["dysfunction_of"].add(filler_label)
            elif on_prop == RO_HAS_MODIFIER:
                buckets["modifiers"].add(filler_label)

        # intersectionOf list case (e.g., equivalentClass [ owl:intersectionOf (...) ])
        for lst in mondo_g.objects(expr, OWL.intersectionOf):
            node = lst
            while node and node != RDF.nil:
                head = mondo_g.value(node, RDF.first)
                if head is None:
                    break
                collect_restriction_fillers(head, buckets)
                node = mondo_g.value(node, RDF.rest)

    # --- step 1: map CUI -> MONDO classes via hasDbXref ---
    cui_to_mondo_uris: Dict[str, Set[URIRef]] = defaultdict(set)

    print("[MONDO] Scanning hasDbXref for UMLS CUIs ...")
    for s, p, o in mondo_g.triples((None, HAS_DBXREF, None)):
        if not isinstance(o, Literal):
            continue
        val = str(o)

        # Typical patterns: "UMLS:C0001234", "UMLS_CUI:C0001234", or just "...C0001234"
        m = re.search(r"(C\d{7})", val)
        if not m:
            continue
        cui = m.group(1)
        if cui not in cuis_of_interest:
            continue

        cui_to_mondo_uris[cui].add(s)

    print(
        f"[MONDO] Found MONDO mappings for {len(cui_to_mondo_uris)} CUIs "
        f"(out of {len(cuis_of_interest)} requested)"
    )

    # --- helper: collect RO-based annotations for a given MONDO class ---
    mondo_cache: Dict[URIRef, Dict[str, Set[str]]] = {}

    def annotations_for_mondo_class(cls: URIRef) -> Dict[str, Set[str]]:
        if cls in mondo_cache:
            return mondo_cache[cls]

        ann: Dict[str, Set[str]] = {
            "mondo_ids": set(),
            "mondo_labels": set(),
            "locations": set(),
            "material_basis_in": set(),
            "features": set(),
            "dysfunction_of": set(),
            "modifiers": set(),
        }

        # MONDO ID as CURIE + main label
        s = str(cls)
        if s.startswith(MONDO_PREFIX):
            local = s.rsplit("/", 1)[-1]      # e.g. "MONDO_0001234"
            num = local.split("MONDO_")[-1]   # "0001234"
            ann["mondo_ids"].add(f"MONDO:{num}")
        ann["mondo_labels"].add(get_label(cls))

        # Direct RO triples (rare, but handle them)
        for p, o in mondo_g.predicate_objects(cls):
            if not isinstance(o, URIRef):
                continue
            filler_label = get_label(o)
            if p == RO_HAS_LOCATION:
                ann["locations"].add(filler_label)
            elif p == RO_HAS_MATERIAL_BASIS_IN:
                ann["material_basis_in"].add(filler_label)
            elif p == RO_HAS_FEATURE:
                ann["features"].add(filler_label)
            elif p == RO_DISEASE_HAS_DYSFUNCTION:
                ann["dysfunction_of"].add(filler_label)
            elif p == RO_HAS_MODIFIER:
                ann["modifiers"].add(filler_label)

        # RO restrictions under rdfs:subClassOf
        for sup in mondo_g.objects(cls, RDFS.subClassOf):
            collect_restriction_fillers(sup, ann)

        # RO restrictions under owl:equivalentClass
        for eq in mondo_g.objects(cls, OWL.equivalentClass):
            collect_restriction_fillers(eq, ann)

        mondo_cache[cls] = ann
        return ann

    # --- step 2: aggregate annotations per CUI ---
    result: Dict[str, Dict[str, Set[str]]] = {}

    for cui, mondo_uris in cui_to_mondo_uris.items():
        agg = {
            "mondo_ids": set(),
            "mondo_labels": set(),
            "locations": set(),
            "material_basis_in": set(),
            "features": set(),
            "dysfunction_of": set(),
            "modifiers": set(),
        }

        for cls in mondo_uris:
            ann = annotations_for_mondo_class(cls)
            for key in agg.keys():
                agg[key].update(ann.get(key, set()))

        result[cui] = agg

    return result


# =====================================================================
# 9. UMLS annotations for disease entities (+ MONDO)
# =====================================================================

def summarize_disease_umls_annotations(
    in_path: Union[str, Path],
    out_dir: Union[str, Path],
    umls_apikey: str,
    prefix: Optional[str] = None,
    orpha_source: str = "ORPHANET",
    umls_base_url: str = "https://uts-ws.nlm.nih.gov/rest",
    mondo_owl: Optional[Union[str, Path]] = None,
    umls_orphanet_csv: Optional[Union[str, Path]] = "../data/enrichment_tables/umls_orphanet.csv",
) -> None:  
    """
    For each disease entity related to AutoimmuneDisease, collect UMLS CUIs
    and aggregate them per CUI, then optionally enrich with MONDO annotations.

    Disease entities =

      - instances having rdf:type AutoimmuneDisease or any subclass
      - Orphanet disease classes (URIs under ORDO) that are subclasses
        of AutoimmuneDisease.

    Orphanet → CUI mapping is taken primarily from the local
    `umls_orphanet.csv` (extracted from MRCONSO.RRF, with SAB=ORPHANET),
    using the CODE column as Orphanet numeric ID.

    Output: <prefix>_disease_umls_cui_mondo_annotations.csv

    Columns:
      cui
      n_disease_instances
      pref_label                 (from UMLS)
      semantic_types             (from UMLS)
      mondo_ids                  (MONDO:nnnnnn)
      mondo_labels               (primary labels of mapped MONDO classes)
      mondo_locations            (labels from RO:0004026)
      mondo_material_basis_in    (labels from RO:0004028)
      mondo_features             (labels from RO:0004025)
      mondo_dysfunction_of       (labels from RO:0004007)
      mondo_modifiers            (labels from RO:0002200)
    """
    if not umls_apikey:
        raise ValueError("summarize_disease_umls_annotations: `umls_apikey` is empty")

    in_path = Path(in_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prefix if prefix is not None else in_path.stem

    g = Graph()
    print(f"[UMLS] Parsing KG from {in_path} ...")
    g.parse(in_path)
    print(f"[UMLS] Parsed KG with {len(g)} triples")

    MAK = "http://makaao.inria.fr/kg/"
    DISEASE_ROOT = MAK + "AutoimmuneDisease"
    MAK_ORPHA_PREFIX = MAK + "orpha_"
    ORDO_PREFIX = "http://www.orpha.net/ORDO/Orphanet_"
    UMLS_PREFIX = "https://uts.nlm.nih.gov/uts/umls/concept/"

    UMLS_AUTH_ENDPOINT = "https://utslogin.nlm.nih.gov/cas/v1/api-key"
    umls_base_url = umls_base_url.rstrip("/")

    # --------- helpers for Orpha / UMLS parsing ---------

    def norm_orpha_id(s: str) -> str:
        """
        Normalize an Orphanet identifier to the numeric code expected by UMLS / MRCONSO.

        Examples:
        "84"                -> "84"
        "000084"            -> "84"
        "Orphanet_000084"   -> "84"
        "orpha_84_instance" -> "84"
        """
        digits = "".join(ch for ch in s if ch.isdigit())
        stripped = digits.lstrip("0")
        return stripped or digits

    def extract_cui_from_uri(uri: str) -> Optional[str]:
        if not uri.startswith(UMLS_PREFIX):
            return None
        tail = uri[len(UMLS_PREFIX):].strip()
        if tail.startswith("C") and tail[1:].isdigit():
            return tail
        return None

    # --------- minimal UMLS REST client ---------

    class UMLSClient:
        def __init__(self, apikey: str, base_url: str, service_url: str):
            self.apikey = apikey
            self.base_url = base_url.rstrip("/")
            self.service_url = service_url
            self.tgt_url: Optional[str] = None

        def _ensure_tgt(self) -> None:
            if self.tgt_url is not None:
                return
            data = {"apikey": self.apikey}
            headers = {"Content-Type": "application/x-www-form-urlencoded"}
            resp = requests.post(UMLS_AUTH_ENDPOINT, data=data, headers=headers, timeout=30)
            resp.raise_for_status()
            m = re.search(r'action="(.+?)"', resp.text)
            if not m:
                raise RuntimeError("Could not find TGT URL in UMLS auth response")
            self.tgt_url = m.group(1)

        def get_ticket(self) -> str:
            self._ensure_tgt()
            assert self.tgt_url is not None
            data = {"service": self.service_url}
            headers = {"Content-Type": "application/x-www-form-urlencoded"}
            resp = requests.post(self.tgt_url, data=data, headers=headers, timeout=30)
            resp.raise_for_status()
            return resp.text.strip()

        def get_source_concept(self, source: str, code: str) -> Dict:
            ticket = self.get_ticket()
            url = f"{self.base_url}/content/current/source/{source}/{code}"
            params = {"ticket": ticket}
            resp = requests.get(url, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            return data.get("result", {}) or {}

        def get_cui_annotations(self, cui: str) -> Tuple[str, List[str]]:
            ticket = self.get_ticket()
            url = f"{self.base_url}/content/current/CUI/{cui}"
            params = {"ticket": ticket}
            resp = requests.get(url, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            result = data.get("result", {}) or {}
            label = (result.get("name") or "").strip()
            semtype_names: List[str] = []
            for st in result.get("semanticTypes", []) or []:
                nm = (st.get("name") or "").strip()
                if nm:
                    semtype_names.append(nm)
            return label, semtype_names

    umls_client = UMLSClient(
        apikey=umls_apikey,
        base_url=umls_base_url,
        service_url=umls_base_url,
    )

    # ---------- Orphanet → CUI mapping from local CSV ----------
    orpha_code_to_cui: Dict[str, str] = {}
    if umls_orphanet_csv is not None:
        try:
            csv_path = Path(umls_orphanet_csv)
            print(f"[UMLS] Loading Orphanet↔CUI mapping from {csv_path} ...")
            with csv_path.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sab = (row.get("SAB") or "").strip()
                    if sab != orpha_source:
                        continue
                    code_raw = (row.get("CODE") or "").strip()
                    cui = (row.get("CUI") or "").strip()
                    if not code_raw or not cui:
                        continue
                    norm_code = norm_orpha_id(code_raw)
                    # If multiple CUIs exist for same code, keep the first one seen.
                    if norm_code and norm_code not in orpha_code_to_cui:
                        orpha_code_to_cui[norm_code] = cui
            print(f"[UMLS] Loaded {len(orpha_code_to_cui)} Orphanet codes with CUI mapping from CSV")
        except FileNotFoundError:
            print(f"[UMLS] WARNING: umls_orphanet_csv not found at {umls_orphanet_csv}; Orphanet mapping will rely only on UMLS API")
        except Exception as e:
            print(f"[UMLS] WARNING: error while reading umls_orphanet_csv {umls_orphanet_csv}: {e}")
    else:
        print("[UMLS] No umls_orphanet_csv provided; Orphanet mapping will rely only on UMLS API")

    # caches
    orpha_to_cui_cache: Dict[str, Optional[str]] = {}
    cui_ann_cache: Dict[str, Tuple[str, List[str]]] = {}

    def map_orpha_to_cui(orpha_id: str) -> Optional[str]:
        """
        Map a normalized Orphanet numeric code to a CUI.

        Priority:
          1) Local MRCONSO-derived CSV (umls_orphanet.csv)
          2) UMLS source API lookup (ORPHANET / orpha_source)
        """
        if not orpha_id:
            return None
        if orpha_id in orpha_to_cui_cache:
            return orpha_to_cui_cache[orpha_id]

        # 1) local CSV
        cui = orpha_code_to_cui.get(orpha_id)
        if cui:
            orpha_to_cui_cache[orpha_id] = cui
            return cui

        # 2) fallback: UMLS source lookup
        try:
            result = umls_client.get_source_concept(orpha_source, orpha_id)
        except Exception:
            orpha_to_cui_cache[orpha_id] = None
            return None

        mapped_cui: Optional[str] = None
        if isinstance(result, dict):
            for key in ("concept", "defaultPreferredConcept", "cui"):
                val = result.get(key)
                if isinstance(val, str) and val.startswith("C") and val[1:].isdigit():
                    mapped_cui = val
                    break

            if mapped_cui is None:
                for coll_key in ("relatedConcepts", "concepts", "atoms"):
                    coll = result.get(coll_key)
                    if not isinstance(coll, list):
                        continue
                    for item in coll:
                        if not isinstance(item, dict):
                            continue
                        for k2 in ("ui", "cui"):
                            v2 = item.get(k2)
                            if isinstance(v2, str) and v2.startswith("C") and v2[1:].isdigit():
                                mapped_cui = v2
                                break
                        if mapped_cui:
                            break
                    if mapped_cui:
                        break

        orpha_to_cui_cache[orpha_id] = mapped_cui
        return mapped_cui

    def get_cui_annotations_cached(cui: str) -> Tuple[str, List[str]]:
        if cui in cui_ann_cache:
            return cui_ann_cache[cui]
        try:
            label, semtypes = umls_client.get_cui_annotations(cui)
        except Exception:
            label, semtypes = "", []
        cui_ann_cache[cui] = (label, semtypes)
        return label, semtypes
    
    def extract_orpha_id_from_uri(uri: str) -> Optional[str]:
        tail = None
        if uri.startswith(MAK_ORPHA_PREFIX):
            tail = uri[len(MAK_ORPHA_PREFIX):]
            if tail.endswith("_instance"):
                tail = tail[: -len("_instance")]
        elif uri.startswith(ORDO_PREFIX):
            tail = uri[len(ORDO_PREFIX):]
        if not tail:
            return None
        oid = norm_orpha_id(tail)
        return oid or None

    # ---------- collect instance types and class hierarchy ----------
    instance_types: Dict[str, Set[str]] = defaultdict(set)
    classes: Set[str] = set()
    subclass_children: Dict[str, Set[str]] = defaultdict(set)

    for s, o in g.subject_objects(RDF.type):
        s_uri = str(s)
        o_uri = str(o)
        instance_types[s_uri].add(o_uri)
        if o in (OWL.Class, RDFS.Class):
            classes.add(s_uri)

    for s, p, o in g.triples((None, RDFS.subClassOf, None)):
        parent_uri = str(o)
        child_uri = str(s)
        subclass_children[parent_uri].add(child_uri)

    # ---------- disease classes ----------
    disease_classes: Set[str] = set()
    if DISEASE_ROOT in classes:
        stack: List[str] = [DISEASE_ROOT]
        visited: Set[str] = set()
        while stack:
            c = stack.pop()
            if c in visited:
                continue
            visited.add(c)
            disease_classes.add(c)
            for child in subclass_children.get(c, ()):
                stack.append(child)

    print(f"[UMLS] Found {len(disease_classes)} disease classes under AutoimmuneDisease")

    # ---------- disease instances ----------
    disease_instances: Set[str] = set()
    for inst_uri, types in instance_types.items():
        if any(t in disease_classes for t in types):
            disease_instances.add(inst_uri)

    print(f"[UMLS] Found {len(disease_instances)} disease instances")

    # ---------- Orphanet disease classes under AutoimmuneDisease ----------
    orpha_disease_classes: Set[str] = {
        c for c in disease_classes
        if c.startswith(ORDO_PREFIX)
    }
    print(f"[UMLS] Found {len(orpha_disease_classes)} Orphanet disease classes under AutoimmuneDisease")

    # Unified set of disease entities (instances + Orphanet disease classes)
    disease_entities: Set[str] = set(disease_instances) | orpha_disease_classes
    print(f"[UMLS] Total disease entities (instances + Orphanet classes): {len(disease_entities)}")

    # ---------- map disease entity -> CUIs ----------
    umls_by_disease: Dict[str, Set[str]] = {d: set() for d in disease_entities}

    # 1) from the disease entity URI itself (Orpha or UMLS)
    n_orpha_mapped = 0
    for d in disease_entities:
        oid = extract_orpha_id_from_uri(d)
        if oid:
            cui = map_orpha_to_cui(oid)
            if cui:
                umls_by_disease[d].add(cui)
                n_orpha_mapped += 1

        cui_direct = extract_cui_from_uri(d)
        if cui_direct:
            umls_by_disease[d].add(cui_direct)

    print(f"[UMLS] Mapped {n_orpha_mapped} Orpha disease entities to at least one CUI")

    # 2) from neighbouring Orpha / UMLS nodes
    for s, p, o in g:
        s_uri = str(s)
        o_uri = str(o)

        # disease as subject
        if s_uri in disease_entities and isinstance(o, URIRef):
            cui = extract_cui_from_uri(o_uri)
            if cui:
                umls_by_disease[s_uri].add(cui)
            oid = extract_orpha_id_from_uri(o_uri)
            if oid:
                cui2 = map_orpha_to_cui(oid)
                if cui2:
                    umls_by_disease[s_uri].add(cui2)

        # disease as object
        if o_uri in disease_entities and isinstance(s, URIRef):
            cui = extract_cui_from_uri(s_uri)
            if cui:
                umls_by_disease[o_uri].add(cui)
            oid = extract_orpha_id_from_uri(s_uri)
            if oid:
                cui2 = map_orpha_to_cui(oid)
                if cui2:
                    umls_by_disease[o_uri].add(cui2)

    # ---------- invert: CUI → set of disease entities ----------
    cui_to_diseases: Dict[str, Set[str]] = defaultdict(set)
    for d, cuis in umls_by_disease.items():
        for cui in cuis:
            cui_to_diseases[cui].add(d)

    print(f"[UMLS] Total distinct CUIs: {len(cui_to_diseases)}")

    all_cuis = sorted(cui_to_diseases.keys())

    # ---------- UMLS annotations per CUI ----------
    cui_pref_labels: Dict[str, str] = {}
    cui_semtypes: Dict[str, List[str]] = {}

    print("[UMLS] Fetching label + semantic types for CUIs ...")
    for i, cui in enumerate(all_cuis, start=1):
        lbl, semtypes = get_cui_annotations_cached(cui)
        cui_pref_labels[cui] = lbl
        cui_semtypes[cui] = semtypes
        if i % 20 == 0:
            print(f"[UMLS] Annotated {i}/{len(all_cuis)} CUIs")

    # ---------- MONDO annotations per CUI (optional) ----------
    mondo_ann: Dict[str, Dict[str, Set[str]]] = {}
    if mondo_owl is not None:
        mondo_ann = build_mondo_annotations_for_cuis(mondo_owl, set(all_cuis))
    else:
        print("[MONDO] No MONDO OWL path provided; MONDO columns will be empty.")

    # ---------- write CUI-centric CSV with MONDO columns ----------
    out_path = out_dir / f"{base}_disease_umls_cui_mondo_annotations.csv"
    print(f"[OUT] Writing CUI+MONDO summary to {out_path}")

    with out_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "cui",
            "n_disease_instances",
            "pref_label",
            "semantic_types",
            "mondo_ids",
            "mondo_labels",
            "mondo_locations",
            "mondo_material_basis_in",
            "mondo_features",
            "mondo_dysfunction_of",
            "mondo_modifiers",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for cui in all_cuis:
            n_d = len(cui_to_diseases[cui])
            pref = cui_pref_labels.get(cui, "")
            semtypes = "|".join(sorted(set(cui_semtypes.get(cui, []))))

            ma = mondo_ann.get(cui, {})
            mondo_ids       = "|".join(sorted(ma.get("mondo_ids", set())))
            mondo_labels    = "|".join(sorted(ma.get("mondo_labels", set())))
            mondo_locations = "|".join(sorted(ma.get("locations", set())))
            mondo_mat_basis = "|".join(sorted(ma.get("material_basis_in", set())))
            mondo_features  = "|".join(sorted(ma.get("features", set())))
            mondo_dysfunc   = "|".join(sorted(ma.get("dysfunction_of", set())))
            mondo_modifiers = "|".join(sorted(ma.get("modifiers", set())))

            writer.writerow({
                "cui": cui,
                "n_disease_instances": n_d,
                "pref_label": pref,
                "semantic_types": semtypes,
                "mondo_ids": mondo_ids,
                "mondo_labels": mondo_labels,
                "mondo_locations": mondo_locations,
                "mondo_material_basis_in": mondo_mat_basis,
                "mondo_features": mondo_features,
                "mondo_dysfunction_of": mondo_dysfunc,
                "mondo_modifiers": mondo_modifiers,
            })

    print("[OUT] Done.")


# =====================================================================
# 10. Main
# =====================================================================

if __name__ == "__main__":
    # Uncomment what you want to run
    # summarize_rdf("../kg/makg-core_v1.rdf", "../plots_and_stats/")
    # summarize_aab_targets("../kg/makg-core_v1.rdf", "../plots_and_stats/")
    summarize_aab_leaf_stats("../kg/makg-core_v1.rdf", "../plots_and_stats/")
    # summarize_aab_diseases("../kg/makg-core_v1.rdf", "../plots_and_stats/")
    # summarize_disease_aab("../kg/makg-core_v1.rdf", "../plots_and_stats/")
    # summarize_disease_instances_aab("../kg/makg-core_v1.rdf", "../plots_and_stats/")
    # summarize_aab_loinc_parts("../kg/makg-core_v1.rdf", "../plots_and_stats/")
    # summarize_orpha_phenotypes("../kg/makg-core_v1.rdf", "../plots_and_stats/")
    # summarize_target_origins("../kg/makg-core_v1.rdf", "../plots_and_stats/")

    #fetch_uniprot_annotations_for_targets(
    #    "../kg/makg-core_v1.rdf",
    #    "../plots_and_stats/",
    #    email="fabien.maury@inserm.fr",
    #)

    # summarize_disease_umls_annotations(
    #     "../kg/makg-core_v1.rdf",
    #     "../plots_and_stats/",
    #     UMLS_API_KEY,
    #     mondo_owl="../data/mondo.owl",
    #     # override if needed:
    #     # umls_orphanet_csv="../data/enrichment_tables/umls_orphanet.csv",
    # )
