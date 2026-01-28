import importlib.util
from pathlib import Path

import pandas as pd


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


class _DummyResp:
    def __init__(self, status_code=404, payload=None):
        self.status_code = status_code
        self._payload = payload or {}

    def json(self):
        return self._payload


def test_02_pipeline_runs_on_minimal_inputs(tmp_path, monkeypatch):
    """
    Runs the 02 script end-to-end on:
      - a tiny MRCONSO.RRF (local UMLS map)
      - a tiny Orphanet XML (orpha->hpo link table)
      - a tiny makaao_core.csv (IDs to enrich)
    while stubbing network lookups.
    """
    mod = load_module(Path("scripts/02_create_enrichment_tables.py"), "enrich02")

    # Redirect all I/O to tmp_path/data
    data_dir = tmp_path / "data"
    enrich_dir = data_dir / "enrichment_tables"
    data_dir.mkdir()
    enrich_dir.mkdir()

    # Patch module paths
    monkeypatch.setattr(mod, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(mod, "ENRICH_DIR", str(enrich_dir))
    monkeypatch.setattr(mod, "XML_PATH", str(data_dir / "en_product4.xml"))
    monkeypatch.setattr(mod, "OUT_ORPHA_LINKS", str(enrich_dir / "orphanet_hpo_links.csv"))
    monkeypatch.setattr(mod, "INPUT_CSV_CORE", str(data_dir / "makaao_core.csv"))
    monkeypatch.setattr(mod, "OUTPUT_CSV_FINAL", str(enrich_dir / "code_names.csv"))
    monkeypatch.setattr(mod, "IN_PATH", str(data_dir / "MRCONSO.RRF"))

    # Stub network calls (UniProt / OLS / etc.)
    monkeypatch.setattr(mod, "req_get", lambda *a, **k: _DummyResp(status_code=404))

    # --- Minimal MRCONSO.RRF ---
    # Need >= 18 fields, with:
    # parts[1]=="ENG", parts[2]=="P", parts[16]=="N", parts[0]=CUI, parts[14]=name
    parts = [""] * 18
    parts[0] = "C0000005"
    parts[1] = "ENG"
    parts[2] = "P"
    parts[14] = "Test Concept"
    parts[16] = "N"
    (data_dir / "MRCONSO.RRF").write_text("|".join(parts) + "|\n", encoding="utf-8")

    # --- Minimal Orphanet XML (one disorder, one HPO association) ---
    xml = """<?xml version="1.0" encoding="UTF-8"?>
<Root>
  <Disorder>
    <OrphaCode>123</OrphaCode>
    <HPODisorderAssociationList>
      <HPODisorderAssociation>
        <HPO>
          <HPOId>HP:0000001</HPOId>
          <HPOTerm>All</HPOTerm>
        </HPO>
      </HPODisorderAssociation>
    </HPODisorderAssociationList>
  </Disorder>
</Root>
"""
    (data_dir / "en_product4.xml").write_text(xml, encoding="utf-8")

    # --- Minimal makaao_core.csv ---
    # include at least one ID per category
    core = pd.DataFrame(
        [
            {
                "uniprot_id": "UP:P12345",
                "umls_id": "C0000005",
                "chebi_id": "CHEBI:23367",
                "disease_id": "ORPHA:123",
                "loinc_part_id": "LP12345-6",
                "loinc_part": "Example LOINC part",
            }
        ]
    )
    core.to_csv(data_dir / "makaao_core.csv", index=False)

    # Run full pipeline
    mod.main()

    # Assertions: outputs exist and parseable
    out_orpha = Path(mod.OUT_ORPHA_LINKS)
    out_codes = Path(mod.OUTPUT_CSV_FINAL)

    assert out_orpha.exists()
    assert out_codes.exists()

    df_orpha = pd.read_csv(out_orpha)
    # columns produced by script after dropping "rank"
    assert set(df_orpha.columns) >= {"orpha_code", "HPOId", "HPOTerm", "frequency"}

    df_codes = pd.read_csv(out_codes)
    assert set(df_codes.columns) == {"source", "id", "name", "url"}
    assert len(df_codes) >= 1
