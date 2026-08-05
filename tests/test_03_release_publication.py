from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import pytest
from rdflib import Graph, OWL, RDF, URIRef


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return load_module(Path("scripts/03_build_kg_from_tables.py"), "build_kg_publication")


def module_records(mod):
    return [
        {
            "name": name,
            "path": f"modules/{name}-module.owl",
        }
        for name in (
            "bao",
            "biolink",
            "chebi",
            "hpo",
            "ordo",
            "prov",
            "sio",
            "skos",
            "uniprot",
        )
    ]


def test_root_catalog_manifest_and_checksums_bind_staged_release(mod, tmp_path):
    stage = tmp_path / "stage"
    reasoning = stage / mod.REASONING_OUTPUT_DIRNAME
    modules = reasoning / "modules"
    modules.mkdir(parents=True)
    records = module_records(mod)
    for record in records:
        (modules / Path(record["path"]).name).write_text(
            "<?xml version='1.0'?><rdf:RDF xmlns:rdf='http://www.w3.org/1999/02/22-rdf-syntax-ns#'/>",
            encoding="utf-8",
        )

    result = {"status": "PASSED_OWL2_DL_AND_HERMIT", "modules": records}
    catalog = mod.write_root_reasoning_catalog(stage, result)
    assert catalog.is_file()
    assert catalog.read_text(encoding="utf-8").count("<uri ") == 9

    graph = Graph()
    graph.add((URIRef("http://example.org/x"), RDF.type, OWL.Class))
    kg = stage / "makaao.rdf"
    tbox = stage / "makaao.owl"
    mod._serialize_rdfxml_verified(graph, kg, "test KG")
    mod._serialize_rdfxml_verified(graph, tbox, "test TBox")

    manifest = {
        "canonical_inputs": {
            "kg": {"sha256": mod._file_sha256(kg)},
            "tbox": {"sha256": mod._file_sha256(tbox)},
        }
    }
    (reasoning / "reasoning-manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )
    (reasoning / "README.txt").write_text("release\n", encoding="utf-8")

    mod.finalize_release_manifest_and_checksums(stage, reasoning, kg, tbox, catalog)
    updated = json.loads((reasoning / "reasoning-manifest.json").read_text())
    assert updated["canonical_inputs"]["root_catalog"]["sha256"] == mod._file_sha256(catalog)
    sums = (reasoning / "SHA256SUMS").read_text(encoding="utf-8")
    assert "../makaao.rdf" in sums
    assert "../makaao.owl" in sums
    assert "../catalog-v001.xml" in sums

    bad_result = {"status": "FAILED", "modules": records}
    with pytest.raises(RuntimeError, match="did not complete"):
        mod.write_root_reasoning_catalog(tmp_path / "bad", bad_result)


def test_transactional_publication_replaces_owned_outputs_and_removes_stale(mod, tmp_path):
    stage = tmp_path / "stage"
    final = tmp_path / "kg"
    (stage / "reasoning" / "reports").mkdir(parents=True)
    (stage / "root.rdf").write_text("new root", encoding="utf-8")
    (stage / "reasoning" / "reports" / "new.txt").write_text("new report", encoding="utf-8")

    (final / "reasoning" / "reports").mkdir(parents=True)
    (final / "root.rdf").write_text("old root", encoding="utf-8")
    (final / "reasoning" / "reports" / "old.txt").write_text("old report", encoding="utf-8")
    (final / "stale.rdf").write_text("stale", encoding="utf-8")
    (final / "unowned.txt").write_text("keep", encoding="utf-8")

    mod.commit_staged_release(stage, final, remove_if_absent=("stale.rdf",))

    assert (final / "root.rdf").read_text() == "new root"
    assert (final / "reasoning" / "reports" / "new.txt").read_text() == "new report"
    assert not (final / "reasoning" / "reports" / "old.txt").exists()
    assert not (final / "stale.rdf").exists()
    assert (final / "unowned.txt").read_text() == "keep"
    assert list(final.glob(".makaao-backup-*")) == []


def test_publication_rolls_back_when_installation_fails(mod, tmp_path, monkeypatch):
    stage = tmp_path / "stage"
    final = tmp_path / "kg"
    stage.mkdir()
    final.mkdir()
    (stage / "a.txt").write_text("new a", encoding="utf-8")
    (stage / "b.txt").write_text("new b", encoding="utf-8")
    (final / "a.txt").write_text("old a", encoding="utf-8")
    (final / "b.txt").write_text("old b", encoding="utf-8")

    real_replace = mod._replace_with_retries
    failed = {"value": False}

    def flaky_replace(source, target, *args, **kwargs):
        if Path(source).name == "b.txt" and Path(source).parent == stage and not failed["value"]:
            failed["value"] = True
            raise PermissionError("simulated installation failure")
        return real_replace(source, target, *args, **kwargs)

    monkeypatch.setattr(mod, "_replace_with_retries", flaky_replace)
    with pytest.raises(PermissionError, match="simulated"):
        mod.commit_staged_release(stage, final)

    assert (final / "a.txt").read_text() == "old a"
    assert (final / "b.txt").read_text() == "old b"
    assert list(final.glob(".makaao-backup-*")) == []


def test_staging_validation_and_java_environment(mod, tmp_path, monkeypatch):
    stage = tmp_path / "stage"
    stage.mkdir()
    (stage / "a").write_text("x", encoding="utf-8")
    mod.validate_staged_root_contents(stage, {"a"})
    with pytest.raises(RuntimeError, match="missing staged outputs"):
        mod.validate_staged_root_contents(stage, {"a", "b"})
    with pytest.raises(RuntimeError, match="unexpected staged outputs"):
        mod.validate_staged_root_contents(stage, set())

    assert mod._strip_controlled_java_options(
        '-Xmx1G -Djava.io.tmpdir=/old -Dkeep="a b"', "JAVA_TOOL_OPTIONS"
    ) == "'-Dkeep=a b'"
    with pytest.raises(ValueError, match="invalid shell-style quoting"):
        mod._strip_controlled_java_options("'unterminated", "JAVA_TOOL_OPTIONS")

    monkeypatch.setenv("JDK_JAVA_OPTIONS", "-Xmx1G -Dkeep=one")
    monkeypatch.setenv("JAVA_TOOL_OPTIONS", "-Djava.io.tmpdir=/old -Dkeep=two")
    monkeypatch.setenv("_JAVA_OPTIONS", "-Xmx2G")
    env = mod.build_reasoner_environment(tmp_path / "java tmp")
    assert "-Xmx1G" not in env.get("JDK_JAVA_OPTIONS", "")
    assert "-Dkeep=one" in env["JDK_JAVA_OPTIONS"]
    assert "-Dkeep=two" in env["JAVA_TOOL_OPTIONS"]
    assert f"-Xmx{mod.JAVA_MAX_HEAP}" in env["JAVA_TOOL_OPTIONS"]
    assert "java.io.tmpdir=" in env["JAVA_TOOL_OPTIONS"]
    assert "_JAVA_OPTIONS" not in env
    assert env["TMPDIR"] == str((tmp_path / "java tmp").resolve())


def test_robot_command_prefers_pinned_jar_then_executable(mod, monkeypatch, tmp_path):
    monkeypatch.setattr(mod, "JAVA_MAX_HEAP", "bad")
    with pytest.raises(ValueError, match="must look like"):
        mod._robot_command()

    # Isolate all lookup locations from the real repository. A developer may
    # legitimately have kg/robot.jar installed, which must not affect this test.
    project_dir = tmp_path / "project"
    kg_dir = tmp_path / "isolated-kg"
    project_dir.mkdir()
    kg_dir.mkdir()
    monkeypatch.setattr(mod, "PROJECT_DIR", project_dir)
    monkeypatch.setattr(mod, "KG_DIR", kg_dir)
    monkeypatch.setattr(mod, "JAVA_MAX_HEAP", "8G")

    pinned_jar = tmp_path / "robot.jar"
    pinned_jar.write_bytes(b"test jar placeholder")
    monkeypatch.setattr(mod, "ROBOT_JAR", str(pinned_jar))

    def with_java_and_robot(name):
        if name == "java":
            return "/usr/bin/java"
        if name == mod.ROBOT_EXECUTABLE:
            return "/usr/bin/robot"
        return None

    monkeypatch.setattr(mod.shutil, "which", with_java_and_robot)
    assert mod._robot_command() == [
        "/usr/bin/java",
        "-Xmx8G",
        "-jar",
        str(pinned_jar),
    ]

    # Once no candidate JAR exists, the executable is the intended fallback.
    pinned_jar.unlink()
    monkeypatch.setattr(mod, "ROBOT_JAR", str(tmp_path / "missing.jar"))
    monkeypatch.setattr(
        mod.shutil,
        "which",
        lambda name: "/usr/bin/robot" if name == mod.ROBOT_EXECUTABLE else None,
    )
    assert mod._robot_command() == ["/usr/bin/robot"]
