import sys
import csv
import importlib.util
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class DummyResponse:
    def __init__(self, status_code=404, payload=None):
        self.status_code = status_code
        self._payload = payload or {}

    def json(self):
        return self._payload


def mrconso_row(*, cui, sab, code, label, scui="", sdui="", tty="PT"):
    fields = [""] * 18
    fields[0] = cui
    fields[1] = "ENG"
    fields[2] = "P"
    fields[6] = "Y"
    fields[9] = scui
    fields[10] = sdui
    fields[11] = sab
    fields[12] = tty
    fields[13] = code
    fields[14] = label
    fields[16] = "N"
    return "|".join(fields) + "|\n"


def configure_paths(mod, monkeypatch, data_dir: Path):
    enrichment = data_dir / "enrichment_tables"
    enrichment.mkdir(parents=True)
    assignments = {
        "DATA_DIR": str(data_dir),
        "ENRICH_DIR": str(enrichment),
        "IN_PATH": str(data_dir / "MRCONSO.RRF"),
        "XML_PATH": str(data_dir / "en_product4.xml"),
        "INPUT_CSV_CORE": str(data_dir / "makaao_core.csv"),
        "LOINC_PART_LINK_CSV": str(data_dir / "LoincPartLink_Primary.csv"),
        "OUT_ORPHA_LINKS": str(enrichment / "orphanet_hpo_links.csv"),
        "OUT_ORPHA_UMLS_MAPPINGS": str(enrichment / "orpha_umls_mappings.json"),
        "OUTPUT_CSV_FINAL": str(enrichment / "code_names.csv"),
        "OUT_LOINC_PART_TESTS": str(enrichment / "loinc_part_test_dict.json"),
        "OUT_LOINC_LABELS": str(enrichment / "loinc_labels.json"),
        "REPORT_PATH": str(enrichment / "enrichment_report.md"),
    }
    for name, value in assignments.items():
        monkeypatch.setattr(mod, name, value)
    return enrichment


def test_normalizers_cover_current_identifier_forms():
    mod = load_module(
        Path("scripts/02_create_enrichment_tables.py"), "enrichment_normalizers"
    )

    assert mod.norm_uniprot("UP:P12345") == "P12345"
    assert mod.norm_umls("CUI:C1234567") == "C1234567"
    assert mod.norm_orpha("https://www.orpha.net/en/disease/detail/123") == "123"
    assert mod.norm_chebi("CHEBI:23367") == "23367"
    assert mod.norm_hpo("http://purl.obolibrary.org/obo/HP_0000001") == "HP:0000001"
    assert mod.norm_loinc_part("https://loinc.org/LP12345-6") == "LP12345-6"
    assert mod.norm_loinc_term("LOINC:1234-5") == "1234-5"


@pytest.mark.filterwarnings("ignore:Testing an element.s truth value.*:DeprecationWarning")
def test_02_pipeline_uses_committed_sample_without_network(
    tmp_path, monkeypatch, sample_core_path, sample_core_frame
):
    mod = load_module(Path("scripts/02_create_enrichment_tables.py"), "enrichment_main")
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    enrichment = configure_paths(mod, monkeypatch, data_dir)
    monkeypatch.setattr(mod, "req_get", lambda *args, **kwargs: DummyResponse())

    # Use the repository's committed, license-safe sample as the real core input.
    # Only external/licensed resources are replaced by tiny local fixtures.
    shutil.copyfile(sample_core_path, data_dir / "makaao_core.csv")
    copied_core = pd.read_csv(
        data_dir / "makaao_core.csv",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )
    assert copied_core.shape == sample_core_frame.shape
    assert list(copied_core.columns) == list(sample_core_frame.columns)

    # Select identifiers that are actually requested by the committed sample.
    # Script 02 intentionally ignores unrelated MRCONSO/Orphanet/LOINC records.
    input_cuis, input_orphas = mod.collect_core_orpha_umls_identifiers()
    input_chebis = mod.collect_core_chebi_identifiers()
    input_loinc_parts, _ = mod.collect_core_loinc_parts_and_labels()

    normalized_columns = {
        "".join(str(column).lower().split()): column
        for column in sample_core_frame.columns
    }
    hpo_column = normalized_columns.get("hpo_id")
    input_hpos = set()
    if hpo_column:
        for value in sample_core_frame[hpo_column].tolist():
            input_hpos.update(mod.extract_hpo_ids(str(value)))

    assert input_cuis, "The committed sample must contain at least one UMLS CUI"
    assert input_orphas, "The committed sample must contain at least one ORPHA disease"
    assert input_chebis, "The committed sample must contain at least one ChEBI identifier"
    assert input_hpos, "The committed sample must contain at least one HPO identifier"
    assert input_loinc_parts, "The committed sample must contain at least one LOINC Part"

    cui = sorted(input_cuis)[0]
    orpha = sorted(input_orphas, key=int)[0]
    chebi = sorted(input_chebis, key=int)[0]
    hpo = sorted(input_hpos)[0]
    loinc_part = sorted(input_loinc_parts)[0]
    loinc_term = "1234-5"

    (data_dir / "en_product4.xml").write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<Root>
  <Disorder>
    <OrphaCode>{orpha}</OrphaCode>
    <Name lang="en">Example disease</Name>
    <HPODisorderAssociationList>
      <HPODisorderAssociation>
        <HPO><HPOId>{hpo}</HPOId><HPOTerm>Example phenotype</HPOTerm></HPO>
        <HPOFrequency><Name lang="en">Frequent</Name></HPOFrequency>
      </HPODisorderAssociation>
    </HPODisorderAssociationList>
  </Disorder>
</Root>
""",
        encoding="utf-8",
    )

    with (data_dir / "LoincPartLink_Primary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "LoincNumber",
                "LongCommonName",
                "PartNumber",
                "PartName",
                "PartCodeSystem",
                "LinkTypeName",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "LoincNumber": loinc_term,
                "LongCommonName": "Example laboratory test",
                "PartNumber": loinc_part,
                "PartName": "Example LOINC part",
                "PartCodeSystem": "https://loinc.org",
                "LinkTypeName": "Primary",
            }
        )

    rows = [
        mrconso_row(
            cui=cui,
            sab="ORPHANET",
            code=orpha,
            scui=f"ORPHA:{orpha}",
            label="Example disease",
        ),
        mrconso_row(
            cui="C1000001", sab="HPO", code=hpo, label="Example phenotype"
        ),
        mrconso_row(
            cui="C1000002", sab="CHEBI", code=chebi, label="Example chemical entity"
        ),
        mrconso_row(
            cui="C1000003", sab="LNC", code=loinc_part, label="Example LOINC part"
        ),
        mrconso_row(
            cui="C1000004", sab="LNC", code=loinc_term, label="Example laboratory test"
        ),
    ]
    (data_dir / "MRCONSO.RRF").write_text("".join(rows), encoding="utf-8")

    mod.main()

    mapping = json.loads(Path(mod.OUT_ORPHA_UMLS_MAPPINGS).read_text(encoding="utf-8"))
    assert mapping["orpha_to_umls"].get(orpha) == [cui]
    assert mapping["umls_to_orpha"].get(cui) == [orpha]

    part_tests = json.loads(Path(mod.OUT_LOINC_PART_TESTS).read_text(encoding="utf-8"))
    assert part_tests[loinc_part] == [loinc_term]

    labels = json.loads(Path(mod.OUT_LOINC_LABELS).read_text(encoding="utf-8"))
    assert labels["parts"][loinc_part]
    assert labels["tests"][loinc_term] == "Example laboratory test"

    code_names = pd.read_csv(mod.OUTPUT_CSV_FINAL, keep_default_na=False)
    assert set(code_names.columns) == {"source", "id", "name", "url"}
    actual_sources = set(code_names["source"])
    assert {"UMLS", "HPO", "ORPHA", "ChEBI", "LOINC"} <= actual_sources

    # When the committed sample contains a UniProt identifier, script 02 must
    # still emit UniProt rows even though live API calls are disabled.
    if "uniprot_id" in normalized_columns:
        column = normalized_columns["uniprot_id"]
        if any(str(value).strip() for value in sample_core_frame[column].tolist()):
            assert "UniProt" in actual_sources

    assert Path(mod.OUT_ORPHA_LINKS).is_file()
    assert Path(mod.REPORT_PATH).is_file()
    assert enrichment.is_dir()

