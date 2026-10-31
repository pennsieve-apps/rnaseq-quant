"""
SPARC dataset metadata (samples.xlsx / subjects.xlsx, or .csv / .tsv) joined onto samples by
their folders, as sc-prep does: a sample whose files are in (or below) a `sam-…` folder gets
that sample row, one in a `pool-…` folder every sample row of that pool, and the subject rows
of those samples' `subject id` (or of a `sub-…` folder above). Pooled rows are combined per
column: one value when they agree, else the distinct values joined with ';'.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from processor.safe_inputs import InputDir, UnsafeInput, assert_defusedxml

MAX_BYTES = 20 * 2**20
MAX_XLSX_UNCOMPRESSED = 200 * 2**20      # zip-bomb guard: sum of the members' declared sizes
SUFFIXES = (".xlsx", ".csv", ".tsv")


def snake(name: str) -> str:
    return re.sub(r"[^0-9a-z]+", "_", str(name).strip().lower()).strip("_")


@dataclass
class Sparc:
    samples: list[dict[str, str]] = field(default_factory=list)
    subjects: dict[str, dict[str, str]] = field(default_factory=dict)
    sample_columns: list[str] = field(default_factory=list)
    subject_columns: list[str] = field(default_factory=list)
    files: list[str] = field(default_factory=list)


def _read_table(path: Path, inp: InputDir) -> list[dict[str, str]]:
    import pandas as pd

    if path.suffix.lower() == ".xlsx":
        assert_defusedxml()                        # openpyxl must parse the (untrusted) XML through defusedxml
        if Path(inp.check(path)).stat().st_size > MAX_BYTES or inp.xlsx_uncompressed_size(path) > MAX_XLSX_UNCOMPRESSED:
            raise UnsafeInput(f"{inp.rel(path)} is larger than {MAX_BYTES // 2**20} MB "
                              f"(or {MAX_XLSX_UNCOMPRESSED // 2**20} MB uncompressed)")
        with inp.open_checked(path) as fh:
            df = pd.read_excel(fh, dtype=str, engine="openpyxl")
    else:
        df = inp.read_csv(path, MAX_BYTES, dtype=str, sep="\t" if path.suffix.lower() == ".tsv" else ",",
                          keep_default_na=False, encoding="utf-8-sig")
    df = df.fillna("")
    df.columns = [snake(c) for c in df.columns]
    rows = [{k: str(v).strip() for k, v in r.items() if k} for r in df.to_dict("records")]
    return [r for r in rows if any(r.values())]


def find(files: list[Path], root: Path, exclude: Path | None) -> dict[str, Path]:
    """The shallowest SPARC samples.* and subjects.* that hold a subject id column."""
    out: dict[str, Path] = {}
    for kind in ("samples", "subjects"):
        cands = [f for f in files if f.stem.lower() == kind and f.suffix.lower() in SUFFIXES
                 and not (exclude and f == exclude)]
        cands.sort(key=lambda p: (len(p.relative_to(root).parts), str(p)))
        if cands:
            out[kind] = cands[0]
    return out


def load(found: dict[str, Path], root: Path) -> Sparc:
    sp = Sparc()
    inp = InputDir(root)
    if "samples" in found:
        rows = _read_table(found["samples"], inp)
        if rows and "subject_id" in rows[0] and "fastq_1" not in rows[0]:
            sp.samples = rows
            sp.sample_columns = list(rows[0])
            sp.files.append(found["samples"].relative_to(root).as_posix())
    if "subjects" in found:
        rows = _read_table(found["subjects"], inp)
        if rows and "subject_id" in rows[0]:
            sp.subjects = {r["subject_id"]: r for r in rows if r.get("subject_id")}
            sp.subject_columns = list(rows[0])
            sp.files.append(found["subjects"].relative_to(root).as_posix())
    return sp


def _combine(rows: list[dict[str, str]], columns: list[str]) -> dict[str, str]:
    out = {}
    for c in columns:
        vals = list(dict.fromkeys(r.get(c, "") for r in rows if r.get(c, "")))
        out[c] = ";".join(vals)
    return out


def metadata_for(sp: Sparc, rel_folder: str, reserved: set[str]) -> dict[str, str]:
    """
    The SPARC columns of a sample whose files are in rel_folder (relative to the input).
    Subject columns that clash with a sample column become subject_<name>; columns that clash
    with this app's own (reserved) become sparc_<name>.
    """
    parts = [p for p in rel_folder.split("/") if p]
    sams = [p for p in parts if p.startswith("sam-")]
    pools = [p for p in parts if p.startswith("pool-")]
    subs = [p for p in parts if p.startswith("sub-")]
    rows = []
    if sams:
        rows = [r for r in sp.samples if r.get("sample_id") == sams[-1]]
    elif pools:
        rows = [r for r in sp.samples if r.get("pool_id") == pools[-1]]
    subject_ids = list(dict.fromkeys(r["subject_id"] for r in rows if r.get("subject_id"))) or subs[-1:]
    subjects = [sp.subjects[s] for s in subject_ids if s in sp.subjects]
    if not rows and not subjects:
        return {}
    out: dict[str, str] = {}
    for k, v in _combine(rows, sp.sample_columns).items():
        out[f"sparc_{k}" if k in reserved else k] = v
    for k, v in _combine(subjects, sp.subject_columns).items():
        if k == "subject_id" and "subject_id" in out:
            continue
        name = f"subject_{k}" if k in sp.sample_columns and not k.startswith("subject_") else k
        out[f"sparc_{name}" if name in reserved else name] = v
    return out
