"""
A tiny synthetic bulk RNA-seq run: six transcripts of four genes, a salmon index built from them
(needs `salmon` on PATH, or SALMON_BIN; or --index copies a prebuilt one), and reads drawn from
known abundances:

  reads/pool-A1/A1_S1_L001_R1_001.fastq.gz, A1_S1_L002_R1_001.fastq.gz   single-end, two lanes
  reads/pool-B1/B1_S2_L001_R1_001.fastq.gz, B1_S2_L001_R2_001.fastq.gz   paired-end
  reads/sam-c1/C1_R1.fastq.gz                                             single-end
  reads/sam-c1/C1_I1.fastq.gz                                             index read (ignored)
  reads/orphan/D1_R2.fastq.gz                                             read 2 only (fails)
  reads/samples.xlsx, reads/subjects.xlsx                                 SPARC metadata
  index/{salmon/, tx2gene.csv, genes.csv, index.json}

  python -m tests.make_fixture data/fixture            # linked: reads/ and index/
  python -m tests.make_fixture data/fixture --flat     # reads at the root next to index/
  python -m tests.make_fixture data/fixture --index tests/fixtures/index   # no salmon needed

tests/fixtures/index is this index as salmon 2.8.0 (the image's) builds it, for a release review,
which runs this generator without salmon.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import random
import shutil
import subprocess
from pathlib import Path

COMP = str.maketrans("ACGT", "TGCA")
# transcript id, gene id, gene name, length
TRANSCRIPTS = [
    ("TX1.1", "G1.1", "Alpha", 1500), ("TX2.1", "G1.1", "Alpha", 900),
    ("TX3.1", "G2.1", "Beta", 1200),
    ("TX4.1", "G3.1", "Gamma", 2000), ("TX5.2", "G3.1", "Gamma", 700),
    ("TX6.1", "G4.1", "Delta", 1000),
]
# relative abundance (molecules) per sample
ABUNDANCE = {
    "A1": [5, 1, 3, 2, 0, 1],
    "B1": [1, 5, 3, 0, 2, 1],
    "C1": [3, 3, 0, 1, 1, 4],
}
READ_LEN, FRAG_MEAN, N_FRAGMENTS = 75, 250, 3000


def sequences(seed: int = 11) -> dict[str, str]:
    rng = random.Random(seed)
    return {t: "".join(rng.choice("ACGT") for _ in range(n)) for t, _, _, n in TRANSCRIPTS}


def build_index(root: Path, seqs: dict[str, str]) -> Path:
    idx = root / "index"
    idx.mkdir(parents=True, exist_ok=True)
    fa = idx / "transcripts.fa"
    fa.write_text("".join(f">{t}\n{s}\n" for t, s in seqs.items()))
    salmon = os.environ.get("SALMON_BIN", "salmon")
    subprocess.run([salmon, "index", "-t", str(fa), "-i", str(idx / "salmon"), "-k", "31", "-p", "2", "-q"],
                   check=True, capture_output=True)
    fa.unlink()
    (idx / "tx2gene.csv").write_text("transcript_id,gene_id,gene_name,gene_type\n" + "".join(
        f"{t},{g},{n},protein_coding\n" for t, g, n, _ in TRANSCRIPTS))
    genes = dict.fromkeys((g, n) for _, g, n, _ in TRANSCRIPTS)
    (idx / "genes.csv").write_text("gene_id,gene_name,gene_type,chrom,length\n" + "".join(
        f"{g},{n},protein_coding,chr1,1000\n" for g, n in genes))
    version = subprocess.run([salmon, "--version"], capture_output=True, text=True, check=True).stdout.split()[-1]
    (idx / "index.json").write_text(json.dumps({
        "organism": "Mus musculus", "source": "synthetic", "assembly": "", "salmon_version": version,
        "kmer_len": 31, "decoys": False, "inputs": {}}, indent=2))
    return idx


def fragments(seqs: dict[str, str], weights: list[int], n: int, rng: random.Random):
    """(fragment sequence) drawn proportional to abundance x length, on either strand."""
    names = [t for t, *_ in TRANSCRIPTS]
    w = [a * len(seqs[t]) for a, t in zip(weights, names)]
    for _ in range(n):
        s = seqs[rng.choices(names, w)[0]]
        flen = min(len(s), max(READ_LEN, int(rng.gauss(FRAG_MEAN, 25))))
        p = rng.randrange(len(s) - flen + 1)
        f = s[p:p + flen]
        yield f if rng.random() < 0.5 else f.translate(COMP)[::-1]


def write_fastq(path: Path, reads: list[str], tag: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for i, r in enumerate(reads):
            fh.write(f"@{tag}:{i}\n{r}\n+\n{'F' * len(r)}\n")


def write_sparc(d: Path) -> None:
    import openpyxl

    sheets = {
        "samples.xlsx": [["sample id", "subject id", "pool id", "sample experimental group", "RIN"],
                         ["sam-a1-1", "sub-1", "pool-A1", "naive", 7.5], ["sam-a1-2", "sub-2", "pool-A1", "naive", 7.5],
                         ["sam-b1-1", "sub-3", "pool-B1", "tumour", 7.1], ["sam-b1-2", "sub-4", "pool-B1", "tumour", 7.1],
                         ["sam-c1", "sub-5", "", "tumour", 6.9]],
        "subjects.xlsx": [["subject id", "subject experimental group", "sex", "species", "date of birth"],
                          ["sub-1", "Control", "female", "Mus musculus", "2020-01-01"],
                          ["sub-2", "Control", "female", "Mus musculus", "2020-01-02"],
                          ["sub-3", "Tumor", "male", "Mus musculus", "2020-01-03"],
                          ["sub-4", "Tumor", "male", "Mus musculus", "2020-01-04"],
                          ["sub-5", "Tumor", "female", "Mus musculus", ""]],
    }
    for name, rows in sheets.items():
        wb = openpyxl.Workbook()
        for r in rows:
            wb.active.append(r)
        wb.save(d / name)


def make_fixture(root: Path, flat: bool = False, seed: int = 42, index: Path | None = None) -> Path:
    seqs = sequences()
    if index:
        shutil.copytree(index, root / "index", dirs_exist_ok=True)
    else:
        build_index(root, seqs)
    reads = root if flat else root / "reads"
    rng = random.Random(seed)
    a1 = [f[:READ_LEN] for f in fragments(seqs, ABUNDANCE["A1"], N_FRAGMENTS, rng)]
    write_fastq(reads / "pool-A1" / "A1_S1_L001_R1_001.fastq.gz", a1[: len(a1) // 2], "A1L1")
    write_fastq(reads / "pool-A1" / "A1_S1_L002_R1_001.fastq.gz", a1[len(a1) // 2:], "A1L2")
    b1 = list(fragments(seqs, ABUNDANCE["B1"], N_FRAGMENTS, rng))
    write_fastq(reads / "pool-B1" / "B1_S2_L001_R1_001.fastq.gz", [f[:READ_LEN] for f in b1], "B1")
    write_fastq(reads / "pool-B1" / "B1_S2_L001_R2_001.fastq.gz", [f.translate(COMP)[::-1][:READ_LEN] for f in b1], "B1")
    c1 = [f[:READ_LEN] for f in fragments(seqs, ABUNDANCE["C1"], N_FRAGMENTS, rng)]
    write_fastq(reads / "sam-c1" / "C1_R1.fastq.gz", c1, "C1")
    write_fastq(reads / "sam-c1" / "C1_I1.fastq.gz", ["ACGTACGT"] * 10, "C1I")
    write_fastq(reads / "orphan" / "D1_R2.fastq.gz", c1[:100], "D1")
    write_sparc(reads)
    return root


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--flat", action="store_true")
    ap.add_argument("--index", type=Path, help="copy this prebuilt index (e.g. tests/fixtures/index) instead of running salmon")
    a = ap.parse_args()
    make_fixture(a.root, a.flat, index=a.index)
