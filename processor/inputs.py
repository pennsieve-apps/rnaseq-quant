"""
What is under INPUT_DIR: the index (a folder with index.json), the FASTQ files grouped into
samples (by Illumina-style names, or by the optional sample sheet), and the files ignored.
"""
from __future__ import annotations

import csv
import gzip
import json
import re
import zlib
from dataclasses import dataclass, field
from pathlib import Path

FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)
INDEX_READ_RE = re.compile(r"_I[12](_\d{3})?$")
# <sample>[_S<n>][_L<lane>][_R1|_R2|_1|_2][_001]
NAME_RE = re.compile(r"^(?P<stem>.+?)(?:_S\d+)?(?:_L(?P<lane>\d{3}))?(?:[_.]R?(?P<read>[12]))?(?:_\d{3})?$")
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
SHEET_COLUMNS = ("sample_id", "fastq_1", "fastq_2")


def sanitise(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", name).strip("_")
    return re.sub(r"_+", "_", s) or "sample"


def rel(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


# --- index ---------------------------------------------------------------------------------

@dataclass
class Index:
    root: Path
    info: dict

    @property
    def salmon(self) -> Path:
        return self.root / "salmon"


def find_index(files: list[Path], root: Path) -> Index:
    """The one folder holding index.json, salmon/ and tx2gene.csv (rnaseq-index's output)."""
    found = [f.parent for f in files if f.name == "index.json"
             and (f.parent / "salmon").is_dir() and (f.parent / "tx2gene.csv").is_file()]
    if not found:
        raise ValueError("no index under INPUT_DIR (a folder with index.json, salmon/ and tx2gene.csv from rnaseq-index)")
    if len(found) > 1:
        raise ValueError(f"{len(found)} indexes under INPUT_DIR ({', '.join(rel(f, root) or '.' for f in found)}); link one")
    info = json.loads((found[0] / "index.json").read_text())
    return Index(found[0], info)


# --- samples -------------------------------------------------------------------------------

@dataclass
class Sample:
    sample_id: str
    r1: list[Path]
    r2: list[Path]
    folder: Path                                   # where its files are (SPARC sam-/pool-/sub- folders above it)
    extra: dict[str, str] = field(default_factory=dict)   # sample sheet columns
    error: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def paired(self) -> bool:
        return bool(self.r2)

    def files(self, root: Path) -> list[str]:
        return [rel(p, root) for p in self.r1 + self.r2]


def parse_name(path: Path) -> tuple[str, str, str] | None:
    """(stem, lane, read) of a FASTQ file name; None for index reads (I1/I2)."""
    base = FASTQ_RE.sub("", path.name)
    if INDEX_READ_RE.search(base):
        return None
    m = NAME_RE.match(base)
    return m.group("stem"), m.group("lane") or "", m.group("read") or "1"


def first_read_name(path: Path) -> str | None:
    """The name of a FASTQ's first record (before the first space, without /1 or /2); None if unreadable."""
    try:
        with open(path, "rb") as fh:
            gz = fh.read(2) == b"\x1f\x8b"
        with (gzip.open(path, "rb") if gz else open(path, "rb")) as fh:
            line = fh.readline(4096)
    except (OSError, EOFError, zlib.error):
        return None
    words = line[1:].split(None, 1)
    if not line.startswith(b"@") or not words:
        return None
    return re.sub(rb"/[12]$", b"", words[0]).decode("ascii", "replace")


def group_fastqs(fastqs: list[Path], root: Path) -> tuple[list[Sample], list[dict]]:
    """
    Samples by file name, per folder: lanes of one sample together, paired when every R1 has its
    R2 and their first records have the same read name (so single-end files that happen to be
    named _1/_2 are not merged). Sample ids must stay unique after sanitising.
    """
    groups: dict[tuple[Path, str], dict[str, dict[str, Path]]] = {}
    ignored = []
    for f in fastqs:
        parsed = parse_name(f)
        if parsed is None:
            ignored.append({"path": rel(f, root), "reason": "index read (I1/I2)"})
            continue
        stem, lane, read = parsed
        reads = groups.setdefault((f.parent, stem), {"1": {}, "2": {}})
        key = lane or f.name
        if lane in reads[read]:
            # same stem, read and lane twice (e.g. sample.fastq.gz next to sample.fq): keep files apart
            key = f.name
        reads[read][key] = f
    named: list[tuple[str, Sample]] = []          # (name before sanitising, sample)
    for (folder, stem), reads in sorted(groups.items(), key=lambda kv: (str(kv[0][0]), kv[0][1])):
        r1 = [reads["1"][k] for k in sorted(reads["1"])]
        r2 = [reads["2"][k] for k in sorted(reads["2"])]
        s = Sample("", r1, [], folder)
        if not r1:
            s.r1, s.error = r2, "only read 2 files (R2/_2) found"
        elif r2 and sorted(_lanes(r1)) == sorted(_lanes(r2)) and len(r1) == len(r2):
            if all(first_read_name(a) is not None and first_read_name(a) == first_read_name(b) for a, b in zip(r1, r2)):
                s.r2 = r2
            else:
                # read 1 and read 2 by name, but not mates: each file is a single-end sample of its own
                for f in r1 + r2:
                    one = Sample("", [f], [], folder)
                    one.notes.append(f"named like read {parse_name(f)[2]} of {stem}, but its first read name does "
                                     "not match its mate file's: quantified as a single-end sample of its own")
                    named.append((FASTQ_RE.sub("", f.name), one))
                continue
        elif r2:
            s.error = f"{len(r1)} read 1 and {len(r2)} read 2 files that do not pair up by lane"
        named.append((stem, s))
    _assign_ids(named, root)
    return [s for _, s in named], ignored


def _assign_ids(named: list[tuple[str, Sample]], root: Path) -> None:
    """sample_id = the sanitised name; where that clashes, the sanitised relative folder path + name;
    a clash left after that fails every sample involved."""
    def ids(full: set[int]) -> list[str]:
        out = []
        for i, (name, s) in enumerate(named):
            folder = rel(s.folder, root) if s.folder != root else ""
            out.append(sanitise(f"{folder}_{name}" if i in full and folder else name))
        return out
    first = ids(set())
    counts = {x: first.count(x) for x in first}
    final = ids({i for i, x in enumerate(first) if counts[x] > 1})
    for i, (_, s) in enumerate(named):
        s.sample_id = final[i]
        if final.count(final[i]) > 1:
            s.error = (f"sample_id {final[i]} is not unique even with its folder path "
                       f"({', '.join(s.files(root))}); give the samples ids in a sample sheet")


def _lanes(files: list[Path]) -> list[str]:
    return [parse_name(f)[1] for f in files]


# --- sample sheet ----------------------------------------------------------------------------

@dataclass
class SampleSheet:
    path: Path | None = None
    columns: list[str] = field(default_factory=list)
    rows: list[dict[str, str]] = field(default_factory=list)


def read_sheet(path: Path | None) -> SampleSheet:
    """
    sample_id, fastq_1, fastq_2 (optional), then metadata; one row per sample, or per lane. A file
    without a fastq_1 column is not a sample sheet (an empty SampleSheet is returned).
    """
    if path is None:
        return SampleSheet()
    if path.stat().st_size > 10 * 2**20:
        raise ValueError(f"{path.name} is larger than 10 MB")
    with open(path, newline="", encoding="utf-8-sig") as fh:
        head = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(head, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        header = [(h or "").strip() for h in reader.fieldnames or []]
        if "fastq_1" not in header:
            return SampleSheet()           # not a sample sheet (e.g. a SPARC samples.csv)
        if "sample_id" not in header:
            raise ValueError(f"{path.name} has no 'sample_id' column (columns: {', '.join(header)})")
        rows = [{(k or "").strip(): (v or "").strip() for k, v in r.items() if k is not None} for r in reader]
    rows = [r for r in rows if r.get("sample_id") or r.get("fastq_1")]
    return SampleSheet(path, [h for h in header if h and h not in SHEET_COLUMNS], rows)


def _match(value: str, fastqs: list[Path], root: Path) -> list[Path]:
    """FASTQ files whose path relative to the input equals or ends with value (never a path outside it)."""
    want = value.replace("\\", "/").strip().strip("/")
    if not want or Path(want).is_absolute() or ".." in Path(want).parts:
        return []
    return [f for f in fastqs if rel(f, root) == want or rel(f, root).endswith("/" + want)]


def samples_from_sheet(sheet: SampleSheet, fastqs: list[Path], root: Path) -> tuple[list[Sample], list[dict]]:
    """Samples as the sheet defines them (rows sharing a sample_id are lanes); errors per bad row."""
    by_id: dict[str, Sample] = {}
    errors = []
    for row in sheet.rows:
        sid = row.get("sample_id", "")
        if not SAMPLE_ID_RE.match(sid):
            errors.append({"sample_id": sid, "reason": "sample_id must be letters, digits, '_' or '-'"})
            continue
        hits = {}
        for col in ("fastq_1", "fastq_2"):
            if not row.get(col):
                continue
            m = _match(row[col], fastqs, root)
            if len(m) != 1:
                errors.append({"sample_id": sid, "reason": f"{col} {row[col]!r} matches {len(m)} FASTQ files under the input"})
                break
            hits[col] = m[0]
        else:
            if "fastq_1" not in hits:
                errors.append({"sample_id": sid, "reason": "no fastq_1"})
                continue
            s = by_id.setdefault(sid, Sample(sid, [], [], hits["fastq_1"].parent,
                                             {c: row.get(c, "") for c in sheet.columns}))
            s.r1.append(hits["fastq_1"])
            if "fastq_2" in hits:
                s.r2.append(hits["fastq_2"])
    for s in by_id.values():
        if s.r2 and len(s.r2) != len(s.r1):
            s.error = "some rows have fastq_2 and some do not"
    return list(by_id.values()), errors
