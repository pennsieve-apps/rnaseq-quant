"""The synthetic fixture through the processor (both layouts; needs salmon and fastp), and the pieces."""
from __future__ import annotations

import csv
import json
import math
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from processor import aggregate as G
from processor import inputs as I
from processor import main as app
from processor import sparc as S
from tests.make_fixture import ABUNDANCE, TRANSCRIPTS, make_fixture

TOOLS = shutil.which(os.environ.get("SALMON_BIN", "salmon")) and shutil.which(os.environ.get("FASTP_BIN", "fastp"))
needs_tools = pytest.mark.skipif(not TOOLS, reason="salmon and fastp not on PATH (set SALMON_BIN / FASTP_BIN)")


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for name in app.DEFAULTS:
        monkeypatch.delenv(name, raising=False)
    for k, v in params.items():
        monkeypatch.setenv(k, v)
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="session")
def outputs(tmp_path_factory) -> dict[str, Path]:
    if not TOOLS:
        pytest.skip("salmon and fastp not on PATH")
    out = {}
    for layout in ("linked", "flat"):
        fx = make_fixture(tmp_path_factory.mktemp(f"fx-{layout}"), flat=layout == "flat")
        d = tmp_path_factory.mktemp(f"out-{layout}")
        mp = pytest.MonkeyPatch()
        try:
            assert run(fx, d, mp) == 0
        finally:
            mp.undo()
        out[layout] = d
    return out


@needs_tools
@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_run_writes_the_counts_contract(outputs, layout):
    out = outputs[layout]
    c = out / "counts"
    assert {p.name for p in c.iterdir()} == {"gene_counts.csv", "gene_tpm.csv", "samples.csv", "samples_full.csv",
                                             "qc", "summary.json", "provenance.json"}
    assert {p.name for p in (c / "qc").iterdir()} == {f"{s}.{k}" for s in ("A1", "B1", "C1")
                                                       for k in ("fastp.json", "fastp.html", "salmon.json")}
    assert {p.name for p in (out / "transcripts").iterdir()} == {"A1", "B1", "C1", "provenance.json"}
    assert (out / "transcripts" / "B1" / "aux_info" / "meta_info.json").is_file()

    table = {r["sample_id"]: r for r in rows(c / "samples.csv")}
    assert list(next(iter(table.values())))[:7] == app.SAMPLE_COLUMNS
    assert table["A1"]["layout"] == "single" and table["A1"]["files"].count(";") == 1     # two lanes
    assert table["B1"]["layout"] == "paired" and table["B1"]["library_type"] == "IU"
    assert table["C1"]["library_type"] == "U"
    assert all(int(r["reads_in"]) == 3000 and float(r["mapping_rate"]) > 0.95 for r in table.values())
    # SPARC: pools combined, allow-listed columns only in samples.csv, everything in samples_full.csv
    assert table["A1"]["sex"] == "female" and table["B1"]["subject_experimental_group"] == "Tumor"
    assert table["C1"]["sample_experimental_group"] == "tumour"
    assert "date_of_birth" not in table["A1"]
    full = {r["sample_id"]: r for r in rows(c / "samples_full.csv")}
    assert full["A1"]["sparc_sample_id"] == "sam-a1-1;sam-a1-2" and full["A1"]["rin"] == "7.5"
    assert full["A1"]["date_of_birth"] == "2020-01-01;2020-01-02"

    counts = pd.read_csv(c / "gene_counts.csv", index_col="gene_id")
    tpm = pd.read_csv(c / "gene_tpm.csv", index_col="gene_id")
    assert list(counts.columns) == ["gene_name", "A1", "B1", "C1"] and list(counts.index) == ["G1.1", "G2.1", "G3.1", "G4.1"]
    assert counts["C1"]["G2.1"] == 0 and counts["A1"]["G2.1"] > 0                  # Beta absent from C1
    assert counts[["A1", "B1", "C1"]].sum().between(2900, 3004).all()          # rounding
    assert np.allclose(tpm[["A1", "B1", "C1"]].sum(), 1e6, rtol=1e-4)
    # TPM tracks the simulated molecule fractions
    truth = {s: {} for s in ABUNDANCE}
    for s, w in ABUNDANCE.items():
        for (_, g, _, _), a in zip(TRANSCRIPTS, w):
            truth[s][g] = truth[s].get(g, 0) + a / sum(w) * 1e6
    for s in ABUNDANCE:
        assert np.corrcoef([truth[s][g] for g in tpm.index], tpm[s])[0, 1] > 0.95

    summary = json.loads((c / "summary.json").read_text())
    assert summary["n_samples"] == 3 and [x["sample_id"] for x in summary["skipped"]] == ["D1"]
    assert {i["reason"] for i in summary["ignored"]} == {"index read (I1/I2)"}
    prov = json.loads((c / "provenance.json").read_text())
    assert prov["app"] == "rnaseq-quant" and prov["version"] == "0.1.0" and prov["n_ok"] == 3 and prov["n_failed"] == 1
    assert {"app", "version", "tools", "parameters", "inputs", "started", "finished", "samples", "peak_rss_mb",
            "citations"} <= set(prov)
    assert {"10.1038/nmeth.4197", "10.12688/f1000research.7563.2"} <= {c["doi"] for c in prov["citations"]}


@needs_tools
def test_sample_sheet_defines_the_samples(tmp_path, monkeypatch):
    fx = make_fixture(tmp_path / "fx")
    (fx / "reads" / "samples.csv").write_text(
        "sample_id,fastq_1,fastq_2,condition\n"
        "lane-test,A1_S1_L001_R1_001.fastq.gz,,ctrl\n"
        "lane-test,A1_S1_L002_R1_001.fastq.gz,,ctrl\n"
        "escape,../../etc/passwd,,x\n")
    assert run(fx, tmp_path / "out", monkeypatch) == 0
    table = rows(tmp_path / "out" / "counts" / "samples.csv")
    assert [r["sample_id"] for r in table] == ["lane-test"] and table[0]["condition"] == "ctrl"
    summary = json.loads((tmp_path / "out" / "counts" / "summary.json").read_text())
    assert [s["sample_id"] for s in summary["skipped"]] == ["escape"]


@needs_tools
def test_index_of_another_salmon_major_is_refused(tmp_path, monkeypatch):
    fx = make_fixture(tmp_path / "fx")
    info = json.loads((fx / "index" / "index.json").read_text())
    info["salmon_version"] = "1.10.3"
    (fx / "index" / "index.json").write_text(json.dumps(info))
    assert run(fx, tmp_path / "out", monkeypatch) == 1


def test_bad_parameters_exit_2(tmp_path, monkeypatch):
    for params in ({"THREADS": "x"}, {"METADATA_COLUMNS": "Sex"}, {"SAMPLE_SHEET": "../x.csv"}):
        assert run(tmp_path, tmp_path / "out", monkeypatch, **params) == 2


def test_no_index_exit_1(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "x_R1.fastq.gz").write_bytes(b"")
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


@pytest.mark.parametrize("name, expected", [
    ("FN1_S21_L007_R1.fastq.gz", ("FN1", "007", "1")),
    ("s1_S3_L001_R2_001.fastq.gz", ("s1", "001", "2")),
    ("sample_1_R1.fq.gz", ("sample_1", "", "1")),
    ("SRR123_2.fastq", ("SRR123", "", "2")),
    ("liver.fastq.gz", ("liver", "", "1")),
    ("s1_S3_L001_I1_001.fastq.gz", None),
])
def test_parse_name(name, expected):
    assert I.parse_name(Path(name)) == expected


def test_grouping(tmp_path):
    names = ["a/x_S1_L001_R1_001.fastq.gz", "a/x_S1_L002_R1_001.fastq.gz", "a/x_S1_L001_R2_001.fastq.gz",
             "a/x_S1_L002_R2_001.fastq.gz", "b/x_R1.fastq.gz", "c/y_R1.fq.gz", "c/y_R2.fq.gz", "d/z_2.fq.gz"]
    for n in names:
        (tmp_path / n).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / n).write_bytes(b"")
    samples, _ = I.group_fastqs(sorted(tmp_path / n for n in names), tmp_path)
    by = {s.sample_id: s for s in samples}
    assert set(by) == {"a_x", "b_x", "y", "z"}                 # same name in two folders: folder prefix
    assert by["a_x"].paired and len(by["a_x"].r1) == 2 and not by["b_x"].paired and by["y"].paired
    assert by["z"].error and not by["a_x"].error


def test_sparc_subject_clash_and_pools():
    sp = S.Sparc(samples=[{"sample_id": "sam-1", "subject_id": "sub-1", "pool_id": "pool-p", "age": "8w"},
                          {"sample_id": "sam-2", "subject_id": "sub-2", "pool_id": "pool-p", "age": "9w"}],
                 subjects={"sub-1": {"subject_id": "sub-1", "age": "10w", "sex": "male"},
                           "sub-2": {"subject_id": "sub-2", "age": "10w", "sex": "male"}},
                 sample_columns=["sample_id", "subject_id", "pool_id", "age"],
                 subject_columns=["subject_id", "age", "sex"])
    m = S.metadata_for(sp, "primary/pool-p", {"sample_id", "files"})
    assert m == {"sparc_sample_id": "sam-1;sam-2", "subject_id": "sub-1;sub-2", "pool_id": "pool-p",
                 "age": "8w;9w", "subject_age": "10w", "sex": "male"}
    assert S.metadata_for(sp, "primary/sub-2/sam-2", set())["age"] == "9w"
    assert S.metadata_for(sp, "elsewhere", set()) == {}


def _tximport_reference(abundance, counts, length, tx2gene):
    """tximport's summarizeToGene + makeCountsFromAbundance('lengthScaledTPM'), written out loop by loop."""
    genes = list(dict.fromkeys(tx2gene[t] for t in abundance.index))
    samples = list(abundance.columns)
    A = {g: [sum(abundance.loc[t, s] for t in abundance.index if tx2gene[t] == g) for s in samples] for g in genes}
    C = {g: [sum(counts.loc[t, s] for t in abundance.index if tx2gene[t] == g) for s in samples] for g in genes}
    L = {}
    for g in genes:
        txs = [t for t in abundance.index if tx2gene[t] == g]
        row = []
        for j, s in enumerate(samples):
            row.append(sum(abundance.loc[t, s] * length.loc[t, s] for t in txs) / A[g][j] if A[g][j] else math.nan)
        if all(math.isnan(x) for x in row):
            ave = sum(sum(length.loc[t, s] for s in samples) / len(samples) for t in txs) / len(txs)
            row = [ave] * len(samples)
        else:
            ok = [x for x in row if not math.isnan(x)]
            geo = math.exp(sum(math.log(x) for x in ok) / len(ok))
            row = [geo if math.isnan(x) else x for x in row]
        L[g] = row
    new = {g: [A[g][j] * sum(L[g]) / len(samples) for j in range(len(samples))] for g in genes}
    out = {}
    for j, s in enumerate(samples):
        lib = sum(C[g][j] for g in genes)
        tot = sum(new[g][j] for g in genes)
        for g in genes:
            out[(g, s)] = new[g][j] * lib / tot
    return out


def test_length_scaled_tpm_matches_tximport():
    tx = ["t1", "t2", "t3", "t4", "t5", "t6"]
    tx2gene = pd.Series({"t1": "g1", "t2": "g1", "t3": "g2", "t4": "g3", "t5": "g3", "t6": "g4"})
    abundance = pd.DataFrame({"s1": [100, 50, 0, 30, 0, 0], "s2": [20, 80, 10, 0, 0, 0], "s3": [60, 0, 5, 10, 40, 0]},
                             index=tx, dtype=float)
    length = pd.DataFrame({"s1": [900, 400, 1500, 2000, 300, 800], "s2": [950, 420, 1450, 1980, 310, 810],
                           "s3": [880, 390, 1550, 2010, 290, 790]}, index=tx, dtype=float)
    counts = abundance * length / 1000 * 3.0
    got, tpm, missing = G.length_scaled_tpm(abundance, counts, length, tx2gene)
    want = _tximport_reference(abundance, counts, length, tx2gene)
    assert missing == []
    for (g, s), v in want.items():
        assert got.loc[g, s] == pytest.approx(v, rel=1e-12)
    assert tpm.loc["g1", "s1"] == 150 and np.allclose(got.sum(), counts.sum())
    # g2 has TPM 0 in s1 only (geometric-mean length), g4 in every sample (average length)
    assert got.loc["g2", "s1"] == 0 and got.loc["g4"].sum() == 0


def test_transcripts_missing_from_tx2gene_are_left_out():
    ab = pd.DataFrame({"s": [1.0, 2.0]}, index=["a", "b"])
    got, _, missing = G.length_scaled_tpm(ab, ab, ab + 100, pd.Series({"a": "g"}))
    assert missing == ["b"] and list(got.index) == ["g"]
