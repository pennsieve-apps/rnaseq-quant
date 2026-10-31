"""
rnaseq-quant: bulk RNA-seq FASTQ to gene counts.

Per sample: fastp (trimming + QC), salmon quant (selective alignment against rnaseq-index's
index, library type detected, sequence bias corrected, fragment GC bias too for paired-end).
Then gene level as tximport's lengthScaledTPM.

  OUTPUT_DIR/counts/gene_counts.csv      gene_id, gene_name, one column per sample (integers)
  OUTPUT_DIR/counts/gene_tpm.csv         same shape, TPM
  OUTPUT_DIR/counts/samples.csv          sample_id, files, layout, reads_in, reads_after_trim,
                                         mapping_rate, library_type, metadata columns
  OUTPUT_DIR/counts/samples_full.csv     sample_id and every metadata column
  OUTPUT_DIR/counts/qc/                  <sample_id>.fastp.json/.html, <sample_id>.salmon.json
  OUTPUT_DIR/counts/summary.json, provenance.json
  OUTPUT_DIR/transcripts/<sample_id>/    salmon's quant.sf, cmd_info.json, lib_format_counts.json,
                                         aux_info/meta_info.json (for tximport/tximeta in R)
  OUTPUT_DIR/transcripts/provenance.json
"""
from __future__ import annotations

import csv
import json
import logging
import os
import re
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

from processor import resources

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("rnaseq-quant")

APP_NAME = "rnaseq-quant"
REPO_ROOT = Path(__file__).resolve().parent.parent
WORK_ROOT = Path("/tmp/rnaseq-quant")
SAMPLE_COLUMNS = ["sample_id", "files", "layout", "reads_in", "reads_after_trim", "mapping_rate", "library_type"]
DEFAULT_METADATA = ("subject_experimental_group,sample_experimental_group,sex,age,species,"
                    "sample_anatomical_location,laterality")
DEFAULTS = {"SAMPLE_SHEET": "samples.csv", "METADATA_COLUMNS": DEFAULT_METADATA, "THREADS": "0"}
CITATIONS = [
    {"doi": "10.1038/nmeth.4197", "pmid": "28263959", "text": "Patro et al. Nature Methods 2017"},
    {"doi": "10.1093/bioinformatics/bty560", "pmid": "30423086", "text": "Chen et al. Bioinformatics 2018"},
    {"doi": "10.1002/imt2.70078", "pmid": "41112039", "text": "Chen. iMeta 2025"},
    {"doi": "10.12688/f1000research.7563.2", "pmid": "26925227", "text": "Soneson et al. F1000Research 2015"},
]
NUM_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) or v).strip() for k, v in DEFAULTS.items()}
    sheet = raw["SAMPLE_SHEET"]
    if "/" in sheet or "\\" in sheet or sheet.startswith(".") or len(sheet) > 200:
        raise ValueError(f"SAMPLE_SHEET must be a file name, got {sheet!r}")
    cols = [c.strip() for c in raw["METADATA_COLUMNS"].split(",") if c.strip()]
    bad = [c for c in cols if not re.fullmatch(r"[a-z0-9_]{1,100}", c)]
    if bad:
        raise ValueError(f"METADATA_COLUMNS must be snake_case column names, got {', '.join(bad)}")
    try:
        threads = int(raw["THREADS"])
    except ValueError as exc:
        raise ValueError(f"THREADS must be an integer, got {raw['THREADS']!r}") from exc
    if threads < 0:
        raise ValueError("THREADS must be >= 0")
    return {"SAMPLE_SHEET": sheet, "METADATA_COLUMNS": cols, "THREADS": threads}


def write_csv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def write_json(path: Path, data) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def major(version: str) -> str:
    return str(version).split(".")[0]


def main() -> int:
    started = datetime.now(timezone.utc)
    t_start = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir or not Path(input_dir).is_dir():
        log.error("INPUT_DIR (an existing folder) and OUTPUT_DIR must be set")
        return 2
    try:
        params = read_parameters()
    except ValueError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    n_threads = resources.threads(params["THREADS"])
    for v in NUM_THREAD_VARS:
        os.environ[v] = str(n_threads)

    import pandas as pd

    from processor import aggregate as G
    from processor import inputs as I
    from processor import pipeline as P
    from processor import sparc as S

    input_dir, output_dir = Path(input_dir), Path(output_dir)

    # 1. what is there
    files = I.walk_files(input_dir)
    try:
        index = I.find_index(files, input_dir)
    except (ValueError, OSError) as exc:
        log.error("%s", exc)
        return 1
    salmon_version = P.tool_version(P.SALMON)
    if major(index.info.get("salmon_version", "")) != major(salmon_version):
        log.error("the index was built with salmon %s; this app runs salmon %s (major versions differ): rebuild it "
                  "with rnaseq-index", index.info.get("salmon_version"), salmon_version)
        return 1
    fastqs = [f for f in files if I.FASTQ_RE.search(f.name) and not f.is_relative_to(index.root)]
    sheet_path = next((f for f in sorted(files, key=lambda p: (len(p.parts), str(p)))
                       if f.name == params["SAMPLE_SHEET"] and not f.is_relative_to(index.root)), None)
    try:
        sheet = I.read_sheet(sheet_path)
    except (ValueError, OSError, UnicodeDecodeError, csv.Error) as exc:
        log.error("cannot use the sample sheet %s: %s", sheet_path, exc)
        return 1
    ignored: list[dict] = []
    if sheet.path:
        samples, sheet_errors = I.samples_from_sheet(sheet, fastqs, input_dir)
        used = {f for s in samples for f in s.r1 + s.r2}
        ignored += [{"path": I.rel(f, input_dir), "reason": "not in the sample sheet"} for f in fastqs if f not in used]
    else:
        samples, ignored = I.group_fastqs(fastqs, input_dir)
        sheet_errors = []
    skipped = [{"sample_id": e["sample_id"], "reason": e["reason"]} for e in sheet_errors]
    for e in sheet_errors:
        log.error("sample sheet row %s: %s", e["sample_id"] or "?", e["reason"])
    try:
        sparc = S.load(S.find(files, input_dir, sheet.path), input_dir)
    except (ValueError, OSError, KeyError) as exc:
        log.warning("SPARC metadata not used: %s", exc)
        sparc = S.Sparc()
    log.info("input: index %s (%s, %s, decoys %s), %d FASTQ files, %d sample(s), sample sheet %s, SPARC %s, "
             "threads %d", I.rel(index.root, input_dir) or ".", index.info.get("organism"), index.info.get("source"),
             index.info.get("decoys"), len(fastqs), len(samples), I.rel(sheet.path, input_dir) if sheet.path else "none",
             ", ".join(sparc.files) or "none", n_threads)
    if not samples:
        log.error("no FASTQ files (.fastq.gz, .fq.gz, .fastq, .fq) under INPUT_DIR")
        return 1

    counts_dir, tx_dir = output_dir / "counts", output_dir / "transcripts"
    qc_dir = counts_dir / "qc"
    for d in (qc_dir, tx_dir):
        d.mkdir(parents=True, exist_ok=True)

    # 2. one sample at a time
    rows, prov_samples, quants = [], [], {}
    reserved = set(SAMPLE_COLUMNS) | set(sheet.columns)
    full_meta: dict[str, dict[str, str]] = {}
    for s in samples:
        t0 = time.monotonic()
        layout = "paired" if s.paired else "single"
        if s.error:
            log.error("%s: skipped: %s", s.sample_id, s.error)
            skipped.append({"sample_id": s.sample_id, "reason": s.error})
            prov_samples.append({"sample_id": s.sample_id, "status": "failed", "error": s.error})
            continue
        try:
            r = P.run_sample(s.sample_id, s.r1, s.r2, index.salmon, n_threads, WORK_ROOT, time.monotonic)
            shutil.copy(r.fastp_json, qc_dir / f"{s.sample_id}.fastp.json")
            shutil.copy(r.fastp_html, qc_dir / f"{s.sample_id}.fastp.html")
            meta = json.loads((r.salmon_dir / "aux_info" / "meta_info.json").read_text())
            lib = json.loads((r.salmon_dir / "lib_format_counts.json").read_text())
            write_json(qc_dir / f"{s.sample_id}.salmon.json", {"meta_info": meta, "lib_format_counts": lib})
            dst = tx_dir / s.sample_id
            (dst / "aux_info").mkdir(parents=True, exist_ok=True)
            for name in ("quant.sf", "cmd_info.json", "lib_format_counts.json", "aux_info/meta_info.json"):
                shutil.copy(r.salmon_dir / name, dst / name)
            shutil.rmtree(r.work, ignore_errors=True)
        except Exception as exc:  # noqa: BLE001 - one sample's failure must not sink the run
            log.error("%s: FAILED: %s\n%s", s.sample_id, exc, traceback.format_exc(limit=3))
            skipped.append({"sample_id": s.sample_id, "reason": str(exc).splitlines()[0]})
            prov_samples.append({"sample_id": s.sample_id, "status": "failed", "error": str(exc)[:2000],
                                 "seconds": round(time.monotonic() - t0, 1)})
            shutil.rmtree(WORK_ROOT / s.sample_id, ignore_errors=True)
            continue
        quants[s.sample_id] = dst / "quant.sf"
        folder = I.rel(s.folder, input_dir) if s.folder != input_dir else ""
        meta_cols = {**S.metadata_for(sparc, folder, reserved), **s.extra}
        full_meta[s.sample_id] = meta_cols
        rows.append({"sample_id": s.sample_id, "files": ";".join(s.files(input_dir)), "layout": layout,
                     "reads_in": r.reads_in, "reads_after_trim": r.reads_after_trim,
                     "mapping_rate": r.mapping_rate, "library_type": r.library_type})
        seconds = round(time.monotonic() - t0, 1)
        prov_samples.append({"sample_id": s.sample_id, "status": "ok", "layout": layout, "files": s.files(input_dir),
                             "reads_in": r.reads_in, "reads_after_trim": r.reads_after_trim,
                             "mapping_rate": r.mapping_rate, "library_type": r.library_type,
                             "seconds_fastp": r.seconds_fastp, "seconds_salmon": r.seconds_salmon, "seconds": seconds})
        log.info("%s: %s, %d file(s), %d reads in, %d after trimming, %.1f%% mapped, library %s, %.0fs "
                 "(fastp %.0fs, salmon %.0fs)", s.sample_id, layout, len(s.r1) + len(s.r2), r.reads_in,
                 r.reads_after_trim, 100 * r.mapping_rate, r.library_type, seconds, r.seconds_fastp, r.seconds_salmon)

    # 3. gene level, tables
    summary: dict = {"n_samples": len(rows), "n_failed": len(skipped)}
    if rows:
        tx2gene = pd.read_csv(index.root / "tx2gene.csv", dtype=str, keep_default_na=False)
        q = G.read_quants(quants)
        counts, tpm, missing = G.length_scaled_tpm(q["TPM"], q["NumReads"], q["EffectiveLength"],
                                                   tx2gene.drop_duplicates("transcript_id").set_index("transcript_id")["gene_id"])
        if missing:
            log.warning("%d quantified transcripts are not in tx2gene.csv and were left out (e.g. %s)",
                        len(missing), ", ".join(missing[:3]))
        names = tx2gene.drop_duplicates("gene_id").set_index("gene_id")["gene_name"]
        order = [g for g in dict.fromkeys(tx2gene["gene_id"]) if g in counts.index]
        for df, path, fmt in ((counts.loc[order].round().astype("int64"), counts_dir / "gene_counts.csv", None),
                              (tpm.loc[order], counts_dir / "gene_tpm.csv", "%.4f")):
            df = df.copy()
            df.insert(0, "gene_name", names.reindex(df.index).to_numpy())
            df.index.name = "gene_id"
            df.to_csv(path, float_format=fmt, lineterminator="\n")
        allowed = [c for c in params["METADATA_COLUMNS"] if any(m.get(c) for m in full_meta.values())]
        extra = sheet.columns + [c for c in allowed if c not in sheet.columns]
        for row in rows:
            row.update({c: full_meta[row["sample_id"]].get(c, "") for c in extra})
        full_cols = list(dict.fromkeys(c for m in full_meta.values() for c, v in m.items() if v or c in sheet.columns))
        write_csv(counts_dir / "samples.csv", SAMPLE_COLUMNS + extra, rows)
        write_csv(counts_dir / "samples_full.csv", ["sample_id"] + full_cols,
                  [{"sample_id": sid, **m} for sid, m in full_meta.items()])
        detected = (counts.loc[order] > 0).sum(axis=0)
        summary.update({
            "n_genes": len(order),
            "median_mapping_rate": float(pd.Series([r["mapping_rate"] for r in rows]).median()),
            "layouts": sorted({r["layout"] for r in rows}),
            "library_types": sorted({r["library_type"] for r in rows}),
            "index": {k: index.info.get(k) for k in ("organism", "source", "assembly", "decoys", "salmon_version")},
            "samples": {r["sample_id"]: {"reads_in": r["reads_in"], "reads_after_trim": r["reads_after_trim"],
                                         "mapping_rate": r["mapping_rate"], "library_type": r["library_type"],
                                         "genes_detected": int(detected[r["sample_id"]]),
                                         "assigned_counts": int(round(counts[r["sample_id"]].sum()))}
                        for r in rows},
        })
    summary["skipped"] = skipped
    summary["ignored"] = ignored
    write_json(counts_dir / "summary.json", summary)

    finished = datetime.now(timezone.utc)
    tools = {"salmon": salmon_version, "fastp": P.tool_version(P.FASTP), "python": sys.version.split()[0]}
    for pkg in ("pandas", "numpy", "openpyxl", "defusedxml"):
        try:
            tools[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    peak = resources.peak_rss_mb()
    provenance = {
        "app": APP_NAME,
        "version": (REPO_ROOT / "VERSION").read_text().strip(),
        "tools": tools,
        "parameters": params,
        "inputs": {"index": I.rel(index.root, input_dir) or ".", "index_info": index.info,
                   "fastq": [I.rel(f, input_dir) for f in fastqs],
                   "sample_sheet": I.rel(sheet.path, input_dir) if sheet.path else None, "sparc": sparc.files},
        "started": started.isoformat(timespec="seconds"),
        "finished": finished.isoformat(timespec="seconds"),
        "samples": prov_samples,
        "skipped": skipped,
        "n_ok": len(rows),
        "n_failed": len(skipped),
        "seconds": round(time.monotonic() - t_start, 1),
        "threads": n_threads,
        "peak_rss_mb": max(peak.values()),
        "peak_rss_mb_by_process": peak,
        "cgroup_peak_mb": resources.cgroup_peak_mb(),
        "citations": CITATIONS,
    }
    for d in (counts_dir, tx_dir):
        write_json(d / "provenance.json", provenance)
    shutil.rmtree(WORK_ROOT, ignore_errors=True)
    log.info("%d sample(s) quantified, %d failed, %.0fs, peak RSS %s MB", len(rows), len(skipped),
             provenance["seconds"], provenance["peak_rss_mb"])
    if not rows:
        log.error("every sample failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
