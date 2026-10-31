"""
Gene-level counts and TPM from salmon's transcript-level quant.sf files, as tximport does with
countsFromAbundance = "lengthScaledTPM" (Soneson, Love & Robinson 2015, F1000Research 4:1521,
doi:10.12688/f1000research.7563.2), following tximport's summarizeToGene() and
makeCountsFromAbundance() (thelovelab/tximport, R/summarizeToGene.R and R/helper.R):

  per gene g and sample j, over the transcripts t of g:
    TPM_gj      = sum_t TPM_tj
    counts_gj   = sum_t NumReads_tj
    length_gj   = sum_t TPM_tj * EffectiveLength_tj / TPM_gj
                  (a gene with TPM 0 in sample j: the geometric mean of its length in the other
                  samples; TPM 0 in every sample: the mean over its transcripts of their
                  across-sample mean effective length)
    new_gj      = TPM_gj * mean_j(length_gj)
    lsTPM_gj    = new_gj * sum_g counts_gj / sum_g new_gj

so the counts keep each sample's library size but carry no isoform-length differences between
samples, and DESeq2 / edgeR need no offset.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

QUANT_COLUMNS = ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]


def read_quants(paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    """TPM, NumReads and EffectiveLength as transcripts x samples; every file must list the same transcripts."""
    cols: dict[str, dict[str, np.ndarray]] = {"TPM": {}, "NumReads": {}, "EffectiveLength": {}}
    names = None
    for sid, p in paths.items():
        df = pd.read_csv(p, sep="\t", usecols=QUANT_COLUMNS,
                         dtype={"Name": str, "Length": np.int64, "EffectiveLength": np.float64,
                                "TPM": np.float64, "NumReads": np.float64})
        if names is None:
            names = df["Name"].to_numpy()
        elif not np.array_equal(names, df["Name"].to_numpy()):
            raise ValueError(f"{sid}: quant.sf lists other transcripts than the first sample (another index?)")
        for c in cols:
            cols[c][sid] = df[c].to_numpy()
    return {c: pd.DataFrame(v, index=names) for c, v in cols.items()}


def length_scaled_tpm(abundance: pd.DataFrame, counts: pd.DataFrame, length: pd.DataFrame,
                      tx2gene: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """
    (gene counts, gene TPM, transcripts left out). tx2gene maps transcript_id to gene_id;
    transcripts it does not have are left out, as tximport does.
    """
    keep = abundance.index.isin(tx2gene.index)
    missing = abundance.index[~keep].tolist()
    abundance, counts, length = abundance[keep], counts[keep], length[keep]
    gene = tx2gene.loc[abundance.index].to_numpy()

    abundance_g = abundance.groupby(gene, sort=False).sum()
    counts_g = counts.groupby(gene, sort=False).sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        length_g = (abundance * length).groupby(gene, sort=False).sum() / abundance_g
    ave_len_gene = length.mean(axis=1).groupby(gene, sort=False).mean()
    length_g = _replace_missing_length(length_g, ave_len_gene)

    new = abundance_g.mul(length_g.mean(axis=1), axis=0)
    scaled = new * (counts_g.sum(axis=0) / new.sum(axis=0))
    return scaled, abundance_g, missing


def _replace_missing_length(length_g: pd.DataFrame, ave_len_gene: pd.Series) -> pd.DataFrame:
    """tximport's replaceMissingLength(): NaN lengths (gene TPM 0 in a sample) filled in."""
    m = length_g.to_numpy(copy=True)
    nan = np.isnan(m)
    all_nan = nan.all(axis=1)
    m[all_nan, :] = ave_len_gene.loc[length_g.index[all_nan]].to_numpy()[:, None]
    some = nan.any(axis=1) & ~all_nan
    if some.any():
        with np.errstate(divide="ignore"):
            logm = np.where(nan[some], np.nan, np.log(np.where(nan[some], 1.0, m[some])))
        geo = np.exp(np.nanmean(logm, axis=1))
        rows = m[some]
        rows[nan[some]] = np.broadcast_to(geo[:, None], rows.shape)[nan[some]]
        m[some] = rows
    return pd.DataFrame(m, index=length_g.index, columns=length_g.columns)
