#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import json
import math
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Union
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib_venn import venn2


# ---------------------------------------------------------------------------
# Uniprot annotation stats
# ---------------------------------------------------------------------------

def plot_uniprot_annotation_stats(
    csv_path: Union[str, Path],
    out_prefix: Union[str, Path],
    top_n_go_terms: int = 20,
    top_n_keywords: int = 20,
    top_n_pathway_dbs: int = 5,
    enrichment_csv_path: Union[str, Path, None] = None,
    top_n_go_pvals: int = 20,
) -> None:
    csv_path = Path(csv_path)
    out_prefix = Path(out_prefix)
    out_dir = out_prefix.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    base = out_prefix.stem

    # global GO + per-aspect GO counters
    go_term_counter = Counter()
    mf_go_term_counter = Counter()  # F = molecular function
    cc_go_term_counter = Counter()  # C = cellular component
    bp_go_term_counter = Counter()  # P = biological process

    keyword_counter = Counter()
    pathway_db_counter = Counter()

    # ---- parse CSV and accumulate counts ----
    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # GO terms: "C:cytosol|P:apoptotic process|F:DNA binding|..."
            raw_go_terms = (row.get("go_terms") or "").strip()
            if raw_go_terms:
                for item in raw_go_terms.split("|"):
                    item = item.strip()
                    if not item:
                        continue

                    aspect_code = ""
                    label = ""

                    if ":" in item:
                        aspect_code, label = item.split(":", 1)
                        aspect_code = aspect_code.strip()
                        label = label.strip()
                    else:
                        label = item

                    if not label:
                        continue

                    # global GO counts
                    go_term_counter[label] += 1

                    # per-aspect counts
                    if aspect_code == "F":
                        mf_go_term_counter[label] += 1
                    elif aspect_code == "C":
                        cc_go_term_counter[label] += 1
                    elif aspect_code == "P":
                        bp_go_term_counter[label] += 1

            # Keywords (names)
            raw_keywords = (row.get("keywords") or "").strip()
            if raw_keywords:
                for kw in raw_keywords.split("|"):
                    kw = kw.strip()
                    if kw:
                        keyword_counter[kw] += 1

            # Pathways: e.g. "KEGG:hsa:402665|Reactome:R-HSA-211728"
            raw_paths = (row.get("pathways") or "").strip()
            if raw_paths:
                for pth in raw_paths.split("|"):
                    pth = pth.strip()
                    if not pth:
                        continue
                    db = pth.split(":", 1)[0].strip()
                    if db:
                        pathway_db_counter[db] += 1

    def _plot_bar(counter: Counter, title: str, out_file: Path, top_n: int) -> None:
        if not counter:
            return
        most_common = counter.most_common(top_n)
        labels = [k for k, _ in most_common]
        counts = [v for _, v in most_common]

        plt.figure(figsize=(10, 6))
        plt.bar(range(len(labels)), counts)
        plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
        plt.ylabel("Count")
        plt.title(title)
        plt.tight_layout()
        plt.savefig(out_file, dpi=150)
        plt.close()

    # aggregated (all aspects together, as before)
    _plot_bar(
        go_term_counter,
        title="Most frequent GO terms (all aspects)",
        out_file=out_dir / f"{base}_top_go_terms.png",
        top_n=top_n_go_terms,
    )

    # Molecular function (F)
    _plot_bar(
        mf_go_term_counter,
        title="Most frequent GO terms – Molecular function (F)",
        out_file=out_dir / f"{base}_top_mf_go_terms.png",
        top_n=top_n_go_terms,
    )

    # Cellular component (C)
    _plot_bar(
        cc_go_term_counter,
        title="Most frequent GO terms – Cellular component (C)",
        out_file=out_dir / f"{base}_top_cc_go_terms.png",
        top_n=top_n_go_terms,
    )

    # Biological process (P)
    _plot_bar(
        bp_go_term_counter,
        title="Most frequent GO terms – Biological process (P)",
        out_file=out_dir / f"{base}_top_bp_go_terms.png",
        top_n=top_n_go_terms,
    )

    _plot_bar(
        keyword_counter,
        title="Most frequent UniProt keywords",
        out_file=out_dir / f"{base}_top_keywords.png",
        top_n=top_n_keywords,
    )

    _plot_bar(
        pathway_db_counter,
        title="Pathway databases represented (Targets)",
        out_file=out_dir / f"{base}_pathway_dbs.png",
        top_n=top_n_pathway_dbs,
    )

    # ------------------------------------------------------------------
    # Optional: GO enrichment plots using g:Profiler p-values
    # ------------------------------------------------------------------
    if enrichment_csv_path is not None:
        enr_path = Path(enrichment_csv_path)
        if enr_path.is_file():
            # collect terms and p-values
            all_terms = []   # (name, p_value, source)
            bp_terms = []
            cc_terms = []
            mf_terms = []

            with enr_path.open("r", newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    name = (row.get("name") or "").strip()
                    src = (row.get("source") or "").strip()
                    p_str = (row.get("p_value") or "").strip()
                    if not name or not p_str:
                        continue
                    try:
                        pval = float(p_str)
                    except ValueError:
                        continue
                    if pval <= 0.0 or math.isnan(pval):
                        continue

                    all_terms.append((name, pval, src))
                    if src == "GO:BP":
                        bp_terms.append((name, pval, src))
                    elif src == "GO:CC":
                        cc_terms.append((name, pval, src))
                    elif src == "GO:MF":
                        mf_terms.append((name, pval, src))

            def _plot_go_pvals(terms, title: str, outfile: Path, top_n: int) -> None:
                if not terms:
                    return
                # sort by p-value ascending (most significant first)
                terms_sorted = sorted(terms, key=lambda x: x[1])[:top_n]
                labels = [t[0] for t in terms_sorted]
                pvals = [t[1] for t in terms_sorted]
                # transform to -log10(p)
                scores = [-math.log10(p) for p in pvals]

                plt.figure(figsize=(10, 6))
                plt.bar(range(len(labels)), scores)
                plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
                plt.ylabel(r"-log$_{10}$(adjusted p-value)")
                plt.title(title)
                plt.tight_layout()
                plt.savefig(outfile, dpi=150)
                plt.close()

            # all GO terms
            _plot_go_pvals(
                all_terms,
                title="Top GO terms by enrichment (−log10 p-value)",
                outfile=out_dir / f"{base}_go_enrichment_all.png",
                top_n=top_n_go_pvals,
            )
            # per aspect
            _plot_go_pvals(
                bp_terms,
                title="Top GO:BP terms by enrichment (−log10 p-value)",
                outfile=out_dir / f"{base}_go_enrichment_bp.png",
                top_n=top_n_go_pvals,
            )
            _plot_go_pvals(
                cc_terms,
                title="Top GO:CC terms by enrichment (−log10 p-value)",
                outfile=out_dir / f"{base}_go_enrichment_cc.png",
                top_n=top_n_go_pvals,
            )
            _plot_go_pvals(
                mf_terms,
                title="Top GO:MF terms by enrichment (−log10 p-value)",
                outfile=out_dir / f"{base}_go_enrichment_mf.png",
                top_n=top_n_go_pvals,
            )




# ---------------------------------------------------------------------------
# HPO descendants of HP:0030057 from hp.json
# ---------------------------------------------------------------------------

HP_JSON = Path("../data/hp.json")                          # input HPO graph (OBO Graph JSON)
ROOT_ID = "HP:0030057"                             # Autoimmune antibody positivity
DESC_CSV = Path("../plots_and_stats/hp0030057_descendants.csv")       # all descendants (ID + label)
MAKAAO_CSV = Path("../data/makaao_core.csv")       # MAKAAO core table
NOT_IN_MAKAAO_CSV = Path("../plots_and_stats/hp0030057_not_in_makaao.csv")
VENN_PNG = Path("../plots_and_stats/hpo_coverage_venn.png")


def normalize_id(s: str) -> str:
    if s.startswith("HP:"):
        return s
    prefix = "http://purl.obolibrary.org/obo/HP_"
    if s.startswith(prefix):
        return "HP:" + s[len(prefix):]
    return s


def load_hpo_graph_and_labels(path: Path):
    if not path.is_file():
        print(f"ERROR: {path} not found. Download hp.json first.", file=sys.stderr)
        sys.exit(1)

    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    graphs = data.get("graphs", [])
    if not graphs:
        print("ERROR: No 'graphs' found in hp.json", file=sys.stderr)
        sys.exit(1)

    g = graphs[0]

    labels = {}
    for n in g.get("nodes", []):
        raw_id = n.get("id")
        if not raw_id:
            continue
        nid = normalize_id(raw_id)
        lbl = n.get("lbl")
        if lbl:
            labels[nid] = lbl

    edges = g.get("edges", [])
    children = defaultdict(set)

    for e in edges:
        if e.get("pred") != "is_a":
            continue
        raw_child = e.get("sub") or e.get("subj")
        raw_parent = e.get("obj")
        if not raw_child or not raw_parent:
            continue
        child = normalize_id(raw_child)
        parent = normalize_id(raw_parent)
        children[parent].add(child)

    return children, labels


def descendants(root_id: str, children_graph):
    root_id = normalize_id(root_id)
    visited = set()
    queue = deque([root_id])

    while queue:
        cur = queue.popleft()
        if cur in visited:
            continue
        visited.add(cur)
        for ch in children_graph.get(cur, ()):
            if ch not in visited:
                queue.append(ch)

    return visited


def write_descendants_csv(hp_json: Path, root_id: str, out_csv: Path) -> None:
    children_graph, labels = load_hpo_graph_and_labels(hp_json)
    all_terms = descendants(root_id, children_graph)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["hpo_id", "label"])
        for term in sorted(all_terms):
            writer.writerow([term, labels.get(term, "")])

    print(f"Wrote {len(all_terms)} descendant terms to {out_csv}")


# ---------------------------------------------------------------------------
# Helpers to read descendants and MAKAAO
# ---------------------------------------------------------------------------

def load_descendants(path: Path):
    ids = set()
    labels = {}
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if "hpo_id" not in reader.fieldnames:
            raise ValueError("descendants CSV must contain 'hpo_id'")
        for row in reader:
            hpo_id = (row.get("hpo_id") or "").strip()
            if not hpo_id:
                continue
            ids.add(hpo_id)
            labels[hpo_id] = (row.get("label") or "").strip()
    return ids, labels


def load_makaao_rows_and_hpos(path: Path):
    """
    Returns:
      total_rows: number of non-empty rows
      mak_ids: set of distinct HPO IDs seen (split on '|')
      row_hpo_lists: list[set[str]] of HPO IDs per row
    """
    total_rows = 0
    mak_ids = set()
    row_hpo_lists = []

    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

        id_like_cols = {"aab_id", "parent_id"}
        non_id_cols = [k for k in fieldnames if k not in id_like_cols]

        for row in reader:
            non_id_values = [(row.get(k) or "").strip() for k in non_id_cols]
            if all(v == "" for v in non_id_values):
                continue

            total_rows += 1

            raw = (row.get("hpo_id") or "").strip()
            row_ids = set()
            if raw:
                for token in raw.split("|"):
                    tok = token.strip()
                    if tok:
                        mak_ids.add(tok)
                        row_ids.add(tok)

            row_hpo_lists.append(row_ids)

    return total_rows, mak_ids, row_hpo_lists


def write_not_in_makaao(
    descendants_csv: Path,
    makaao_csv: Path,
    out_csv: Path,
):
    desc_ids, labels = load_descendants(descendants_csv)
    _, mak_ids, _ = load_makaao_rows_and_hpos(makaao_csv)

    missing = sorted(desc_ids - mak_ids)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["hpo_id", "label"])
        for h in missing:
            writer.writerow([h, labels.get(h, "")])

    print(f"{len(missing)} descendants of {ROOT_ID} not in makaao_core.csv.")
    return desc_ids, mak_ids, set(missing)


# ---------------------------------------------------------------------------
# Venn diagram using matplotlib-venn
#   Left set  = distinct HPO IDs in MAKAAO
#   Right set = descendants of HP:0030057
#   Intersection = descendant HPO terms mentioned in MAKAAO
#   Extra text under the plot: row counts and rows with ≥1 descendant term
# ---------------------------------------------------------------------------

def plot_hpo_venn(
    descendants_csv: Union[str, Path],
    makaao_csv: Union[str, Path],
    out_path: Union[str, Path],
) -> None:
    descendants_csv = Path(descendants_csv)
    makaao_csv = Path(makaao_csv)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # HPO term sets
    desc_ids, _labels = load_descendants(descendants_csv)
    total_rows, mak_ids, row_hpo_lists = load_makaao_rows_and_hpos(makaao_csv)

    desc_used = desc_ids & mak_ids
    rows_with_desc = sum(
        1 for row_ids in row_hpo_lists if any(h in desc_ids for h in row_ids)
    )

    nA = len(mak_ids)        # HPO IDs in MAKAAO
    nB = len(desc_ids)       # descendants of HP:0030057
    nAB = len(desc_used)     # descendants used in MAKAAO

    # sizes for venn2: (only A, only B, both)
    subsets = (nA - nAB, nB - nAB, nAB)

    fig, ax = plt.subplots(figsize=(6, 6))
    v = venn2(
        subsets=subsets,
        set_labels=("HPO IDs in MAKAAO", "HP:0030057 descendants"),
        ax=ax,
    )

    # Make numbers readable
    for region_id in ("10", "01", "11"):
        lbl = v.get_label_by_id(region_id)
        if lbl is not None:
            lbl.set_fontsize(11)

    # Title
    ax.set_title(
        "Overlap between HPO terms in MAKAAO and HP:0030057 descendants",
        fontsize=11,
        pad=20,
    )

    # Extra text below the diagram
    ax.text(
        0.0,
        -0.8,
        f"Rows in MAKAAO: {total_rows}   "
        f"Distinct HPO IDs in MAKAAO: {nA}   "
        f"Descendant HPO terms (incl. HP:0030057): {nB}",
        fontsize=9,
        ha="center",
        va="center",
    )
    ax.text(
        0.0,
        -0.95,
        f"Rows with ≥1 descendant term: {rows_with_desc}",
        fontsize=9,
        ha="center",
        va="center",
    )

    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved Venn diagram to {out_path}")

#########################################################################################
# bar plot for aab taxonomy leaves

def plot_aab_leaf_barplots(
    csv_path: Union[str, Path],
    out_prefix: Union[str, Path],
    top_n: int = 30,
) -> None:
    """
    Read makg-core_v1_aab_leaf_stats.csv and produce:
      1) Barplot of top N AAB classes by number of leaf subclasses.
      2) Barplot of the distribution of leaf-subclass counts across classes.

    The class(es) with the largest number of leaves are excluded from BOTH plots
    (typically the root Autoantibody class).
    """
    csv_path = Path(csv_path)
    out_prefix = Path(out_prefix)
    out_dir = out_prefix.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    base = out_prefix.stem

    # ---- read CSV ----
    records = []  # (label, n_leaf)
    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = (row.get("aab_label") or "").strip()
            n_leaf_str = (row.get("n_leaf_aab_subclasses") or "").strip()
            if not n_leaf_str:
                continue
            try:
                n_leaf = int(n_leaf_str)
            except ValueError:
                continue

            # If label is empty, fall back to local name from URI
            if not label:
                uri = (row.get("aab_class_uri") or "").strip()
                if uri:
                    label = uri.rsplit("/", 1)[-1]
                else:
                    label = "<?>"

            records.append((label, n_leaf))

    if not records:
        print(f"No usable rows in {csv_path}")
        return

    # ---- exclude the class(es) with the largest number of leaves ----
    max_leaves = max(n for _, n in records)
    filtered = [(lbl, n) for (lbl, n) in records if n != max_leaves]

    if not filtered:
        print(
            f"After excluding max-leaf classes (n={max_leaves}), no records remain; "
            f"skipping plots."
        )
        return

    excluded_labels = [lbl for (lbl, n) in records if n == max_leaves]
    print(
        f"Excluding {len(excluded_labels)} class(es) with max leaves = {max_leaves}: "
        + ", ".join(excluded_labels)
    )

    # ------------------------------------------------------------------
    # 1) Top-N barplot (classes with most leaf subclasses), on FILTERED
    # ------------------------------------------------------------------
    records_sorted = sorted(filtered, key=lambda x: (-x[1], x[0]))[:top_n]
    labels = [lbl for lbl, _ in records_sorted]
    counts = [n for _, n in records_sorted]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(range(len(labels)), counts)
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.ylabel("Number of leaf AAB subclasses")
    plt.title(
        f"Top {len(labels)} AAB classes by number of leaf subclasses\n"
        f"(excluding max-leaf class(es))"
    )

    # grid on y-axis
    ax = plt.gca()
    ax.set_axisbelow(True)
    plt.grid(axis="y", linestyle="--", linewidth=0.5)

    # value labels above each bar (for number of leaves)
    for i, val in enumerate(counts):
        plt.text(
            i,
            val,
            str(val),
            ha="center",
            va="bottom",
            fontsize=8,
        )

    plt.tight_layout()
    plt.savefig(out_dir / f"{base}_top_leaf_aab_classes.png", dpi=150)
    plt.close()

    # ------------------------------------------------------------------
    # 2) Distribution of n_leaf_aab_subclasses, on the SAME FILTERED SET
    # ------------------------------------------------------------------
    all_counts = [n for _, n in filtered]
    max_n = max(all_counts)

    # Define bins: either fine-grained if small max, or coarse if large max
    if max_n <= 20:
        bins = [(i, i) for i in range(0, max_n + 1)]
    else:
        bins = [
            (0, 0),
            (1, 1),
            (2, 5),
            (6, 10),
            (11, 20),
            (21, max_n),
        ]

    bin_labels = []
    bin_values = []
    for lo, hi in bins:
        if lo == hi:
            label = f"{lo}"
        elif hi == max_n and hi > lo:
            label = f"{lo}–{hi}" if hi <= 50 else f"{lo}+"
        else:
            label = f"{lo}–{hi}"
        cnt = sum(1 for v in all_counts if lo <= v <= hi)
        bin_labels.append(label)
        bin_values.append(cnt)

    plt.figure(figsize=(8, 5))
    plt.bar(range(len(bin_labels)), bin_values)
    plt.xticks(range(len(bin_labels)), bin_labels, rotation=0, ha="center")
    plt.ylabel("Number of AAB classes")
    plt.xlabel("Number of leaf AAB subclasses")
    plt.title(
        "Distribution of leaf-subclass counts across AAB classes\n"
        "(excluding max-leaf class(es))"
    )

    ax = plt.gca()
    ax.set_axisbelow(True)
    plt.grid(axis="y", linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.savefig(out_dir / f"{base}_leaf_aab_distribution.png", dpi=150)
    plt.close()



# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def main():
    plot_uniprot_annotation_stats(
        csv_path="../plots_and_stats/makg-core_v1_uniprot_annotations.csv",
        out_prefix="../plots_and_stats/makg-core_v1_uniprot",
        top_n_go_terms=25,
        top_n_keywords=20,
        enrichment_csv_path="../plots_and_stats/makg-core_v1_go_enrichment_gprofiler.csv",
        top_n_go_pvals=25,
    )

    plot_aab_leaf_barplots(
        csv_path="../plots_and_stats/makg-core_v1_aab_leaf_stats.csv",
        out_prefix="../plots_and_stats/makg-core_v1_aab_leaf",
        top_n=30,  # adjust as you like
    )

    write_descendants_csv(
        hp_json=HP_JSON,
        root_id=ROOT_ID,
        out_csv=DESC_CSV,
    )

    desc_ids, mak_ids, _missing = write_not_in_makaao(
        descendants_csv=DESC_CSV,
        makaao_csv=MAKAAO_CSV,
        out_csv=NOT_IN_MAKAAO_CSV,
    )
    print(f"Total descendants: {len(desc_ids)}")
    print(f"Distinct HPO IDs in MAKAAO: {len(mak_ids)}")

    plot_hpo_venn(
        descendants_csv=DESC_CSV,
        makaao_csv=MAKAAO_CSV,
        out_path=VENN_PNG,
    )


if __name__ == "__main__":
    main()
