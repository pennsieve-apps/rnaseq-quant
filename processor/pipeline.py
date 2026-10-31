"""
One sample: fastp (adapter and quality trimming, QC report), then salmon quant against the index.
Everything happens in a work folder under /tmp that is removed afterwards.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

FASTP = os.environ.get("FASTP_BIN", "fastp")
SALMON = os.environ.get("SALMON_BIN", "salmon")
FASTP_MAX_THREADS = 16          # fastp's own ceiling
LOG_TAIL = 20


class StepError(Exception):
    pass


@dataclass
class Result:
    work: Path
    fastp_json: Path
    fastp_html: Path
    salmon_dir: Path
    reads_in: int
    reads_after_trim: int
    mapping_rate: float
    library_type: str
    seconds_fastp: float
    seconds_salmon: float


def _is_gzip(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _concat(files: list[Path], out: Path) -> Path:
    """Lanes of one read as one file (gzip members concatenate into a valid gzip file)."""
    if len(files) == 1:
        return files[0]
    if len({_is_gzip(f) for f in files}) > 1:
        raise StepError("a sample's lanes mix gzipped and plain FASTQ")
    out = out.with_suffix(".fastq.gz" if _is_gzip(files[0]) else ".fastq")
    with open(out, "wb") as dst:
        for f in files:
            with open(f, "rb") as src:
                shutil.copyfileobj(src, dst, 1 << 22)
    return out


def _run(cmd: list[str], log_path: Path) -> None:
    with open(log_path, "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False,
                            env={**os.environ, "NO_COLOR": "1"}).returncode
    if rc != 0:
        tail = log_path.read_text(errors="replace").splitlines()[-LOG_TAIL:]
        raise StepError(f"{Path(cmd[0]).name} exited with {rc}:\n  " + "\n  ".join(tail))


def tool_version(binary: str) -> str:
    p = subprocess.run([binary, "--version"], capture_output=True, text=True, check=True)
    return (p.stdout + p.stderr).split()[-1]


def run_sample(sample_id: str, r1: list[Path], r2: list[Path], index: Path, threads: int, work_root: Path,
               timer) -> Result:
    work = work_root / sample_id
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    paired = bool(r2)
    in1 = _concat(r1, work / "in_R1")
    in2 = _concat(r2, work / "in_R2") if paired else None

    t0 = timer()
    fastp_json, fastp_html = work / "fastp.json", work / "fastp.html"
    trim1, trim2 = work / "trim_R1.fastq.gz", work / "trim_R2.fastq.gz"
    cmd = [FASTP, "-i", str(in1), "-o", str(trim1), "-j", str(fastp_json), "-h", str(fastp_html),
           "-R", f"{sample_id} fastp report", "-w", str(min(threads, FASTP_MAX_THREADS)), "-z", "1"]
    if paired:
        cmd += ["-I", str(in2), "-O", str(trim2), "--detect_adapter_for_pe"]
    _run(cmd, work / "fastp.log")
    for f in (in1, in2):
        if f is not None and f.parent == work:
            f.unlink()
    seconds_fastp = timer() - t0

    t1 = timer()
    out = work / "salmon"
    cmd = [SALMON, "quant", "-i", str(index), "-l", "A", "--seqBias", "-p", str(threads), "-o", str(out)]
    cmd += (["-1", str(trim1), "-2", str(trim2), "--gcBias"] if paired else ["-r", str(trim1)])
    _run(cmd, work / "salmon.log")
    seconds_salmon = timer() - t1
    for f in (trim1, trim2):
        f.unlink(missing_ok=True)

    fp = json.loads(fastp_json.read_text())["summary"]
    meta = json.loads((out / "aux_info" / "meta_info.json").read_text())
    lib = json.loads((out / "lib_format_counts.json").read_text())
    per_fragment = 2 if paired else 1
    return Result(work, fastp_json, fastp_html, out,
                  reads_in=fp["before_filtering"]["total_reads"] // per_fragment,
                  reads_after_trim=fp["after_filtering"]["total_reads"] // per_fragment,
                  mapping_rate=round(float(meta["percent_mapped"]) / 100, 4),
                  library_type=str(lib.get("expected_format", "")),
                  seconds_fastp=round(seconds_fastp, 1), seconds_salmon=round(seconds_salmon, 1))
