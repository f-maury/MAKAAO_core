import importlib.util
from pathlib import Path

def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod

def test_end_to_end_build_sample_kg(tmp_path):
    # 1) process sample core into processed_tables
    p1 = load_module(Path("scripts/01_process_makaao_core_to_tables.py"), "process_core_smoke")
    df = p1.load_core(Path("data/makaao_sample.csv"))

    processed = tmp_path / "processed_tables"
    processed.mkdir()
    # write the full set of processed tables expected by builder
    # (safe to call these even if some columns in sample are sparse)
    p1.write_index_name_en(df, processed)
    p1.write_index_hpo_id(df, processed)
    p1.write_index_parent_index(df, processed)
    p1.write_index_syn_en(df, processed)
    p1.write_index_syn_fr(df, processed)
    p1.write_index_cui_source(df, processed)
    p1.write_index_disease_source(df, processed)
    p1.write_index_uniprot_source(df, processed)
    p1.write_index_chebi_source(df, processed)
    p1.write_index_loinc(df, processed)

    # 2) build KG from processed tables + sample core, with enrichment inputs absent
    p3 = load_module(Path("scripts/03_build_kg_from_tables.py"), "build_kg_smoke")
    p3.BASE_DIR = str(processed) + "/"
    sample_core = tmp_path / "makaao_core.csv"
    sample_core.write_bytes(Path("data/makaao_sample.csv").read_bytes())
    p3.makaao_core_name = str(sample_core)

    out_dir = tmp_path / "kg"
    out_dir.mkdir()
    p3.OUTPUT_OWL_ENRICHED = str(out_dir / "makg-sample.rdf")

    # run main; should tolerate missing enrichment tables
    p3.main()

    out_path = out_dir / "makg-sample.rdf"
    assert out_path.exists()
    assert out_path.stat().st_size > 0
