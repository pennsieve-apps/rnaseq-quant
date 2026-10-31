# rnaseq-quant

Turns bulk RNA-seq FASTQ files into a gene count table: the first step of workflow A of the
bulk RNA-seq workflows, `rnaseq-quant (+ index from rnaseq-index) → bulk-de → bulk-report`.

**For a biologist.** Each sample's reads are cleaned (adapter remnants and low-quality ends
removed, with a QC report per sample), then matched to the organism's known transcripts and
counted. Counts are summed to genes in a way that keeps each sample's sequencing depth but
removes differences in which isoform of a gene was expressed, so the table can go straight
into differential expression (DESeq2/PyDESeq2, edgeR). You also get, per sample, how many
reads came in, how many survived trimming, what fraction matched transcripts, and the
library's strandedness as detected from the data.

**For a computational analyst.** Per sample: fastp 1.4.0 (default adapter detection,
`--detect_adapter_for_pe` for paired-end, default quality/length filters; lanes concatenated
first), then `salmon quant -l A --seqBias` (+ `--gcBias` for paired-end) against an
rnaseq-index index (salmon 2.8.0: selective alignment, deterministic two-phase inference).
Gene level is tximport's `countsFromAbundance = "lengthScaledTPM"` (Soneson et al. 2015),
re-implemented in `processor/aggregate.py` from tximport's `summarizeToGene()` and
`makeCountsFromAbundance()` (thelovelab/tximport, devel commit `2c53267`) and checked against
tximport 1.40.0 itself (see "Results"; `tests/golden/` keeps a 30-gene cut of that comparison as a
test).
Counts are rounded to integers. This follows established practice (Conesa et al. 2016:
quality control and trimming, transcript-level quantification; Love et al. 2015 / Soneson et
al. 2015: transcript estimates summarised to genes for count-based DE).

## Ports

| Port | Direction | What |
|---|---|---|
| `reads` | input (subfolder), required | FASTQ (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) anywhere below; optional `samples.csv` sample sheet; optional SPARC `samples.xlsx` / `subjects.xlsx` (or `.csv`/`.tsv`). |
| `index` | input (subfolder), required | rnaseq-index's output (a folder with `index.json`, `salmon/`, `tx2gene.csv`); normally a workspace reference. Its salmon index format (`index_version`, `sshash_format_version` in `salmon/info.json`) must be the one this image's salmon writes (checked at start by indexing two short sequences; salmon 2.8: 2 / 5). salmon 2.9 moved to format 3, so an index from a 2.9-based rnaseq-index is refused here (exit 1), and vice versa. |
| `counts` | output | `OUTPUT_DIR/counts/`: gene tables, sample table, QC; feeds bulk-de. |
| `transcripts` | output | `OUTPUT_DIR/transcripts/`: salmon's per-sample transcript estimates. |

Without links, everything may sit in `INPUT_DIR`: the index is the folder holding
`index.json`; every FASTQ outside it is a read file.

**Samples.** Grouped by file name within a folder: `<sample>_S<n>_L<lane>_R<1|2>_001`,
`<sample>_R1/_R2`, `<sample>_1/_2`, or no read tag (single-end). Lanes of one sample are
concatenated; a sample is paired-end when every R1 has its R2 (same lane) **and** the first
records of R1 and R2 carry the same read name (after the first space, without `/1` `/2`);
otherwise each file is a single-end sample of its own, with a warning (single-end runs named
`SRR…_1` / `SRR…_2` are not merged into a fake pair). A sample with only R2, or with R1/R2 that
do not pair up by lane, fails (logged, skipped). Index reads (`_I1`, `_I2`) are ignored.
Sample ids are the sanitised names; when two clash, each gets its relative folder path as a
prefix, and a clash left after that fails both (give ids in a sample sheet). A sample sheet (`SAMPLE_SHEET`, default `samples.csv`, with
`sample_id, fastq_1, fastq_2`, then any metadata columns; one row per sample or per lane)
replaces the grouping: only its samples are run, its paths are matched against the files
found under the input (equal to, or the end of, their relative path; absolute paths and `..`
never match), and FASTQs it does not list are reported as ignored.

**SPARC metadata**, as sc-prep reads it: a sample whose files are in (or below) a `sam-…`
folder gets that row of `samples.xlsx`; one in a `pool-…` folder gets every row with that
`pool id` (pooled libraries, as in the demo data); the subjects come from those rows'
`subject id` (or a `sub-…` folder above). Column names become snake_case; subject columns that
clash with a sample column become `subject_<name>`; SPARC columns that clash with this app's
own (e.g. SPARC's `sample id`) become `sparc_<name>`. Rows of a pool are combined per column
(one value when they agree, else the distinct values joined by `;`). Only `METADATA_COLUMNS`
go into `samples.csv`; every column goes to `samples_full.csv`. Sample sheet columns go to
both and win over SPARC ones.

## Parameters

| Name | Default | Effect |
|---|---|---|
| `SAMPLE_SHEET` | `samples.csv` | File name of the optional sample sheet (shallowest wins). A file of that name without a `fastq_1` column (e.g. a SPARC `samples.csv`) is not a sheet. |
| `METADATA_COLUMNS` | `subject_id,subject_experimental_group,sample_experimental_group,sex,age,species,sample_anatomical_location,laterality` | SPARC columns copied into `samples.csv` (snake_case). A pooled library gets every pooled value, `;`-joined (e.g. `subject_id` = `sub-M7;sub-M8`). |
| `THREADS` | `0` | Threads for fastp (at most 16) and salmon; 0 = the container's CPU quota. |

Bad values exit 2. No index, an index of another salmon index format, or no FASTQ: exit 1.
A failing sample (bad input) is logged and skipped; the run fails (exit 1) when every sample
failed. When fastp or salmon is **killed** (signal, exit 137: in practice the out-of-memory
killer), the run stops at once with exit 1 and writes no gene tables: the remaining samples
would fail the same way, and a table missing some of them would hide it.
There is no seed: fastp is deterministic and salmon 2.8's default inference is byte-identical
across runs and thread counts.

**Untrusted input.** Files are read only inside the input (`processor/safe_inputs.py`, shared
verbatim by the workflow apps): symlinks are followed while they resolve inside `INPUT_DIR` or
inside the targets of the links the platform places there; others are skipped and logged.
Sample-sheet paths only ever match files found that way. SPARC `.xlsx` files are parsed only
through defusedxml (checked at run time) and refused above 20 MB on disk or 200 MB uncompressed
(zip-bomb guard).

## Outputs

`counts/` (the bulk RNA-seq conventions' Counts contract):

- `gene_counts.csv`: `gene_id, gene_name`, then one column per `sample_id`: lengthScaledTPM
  counts, rounded to integers (so DESeq2 needs no offsets). Genes in `genes.csv`/`tx2gene.csv`
  order.
- `gene_tpm.csv`: same shape, TPM (transcript TPM summed per gene).
- `samples.csv`: `sample_id`, `files` (relative paths, `;`-separated), `layout` (`single`,
  `paired`), `reads_in` and `reads_after_trim` (reads, or read pairs for paired-end, from
  fastp), `mapping_rate` (fraction of trimmed fragments salmon mapped), `library_type`
  (salmon's detected type: `U`, `SF`, `SR`, `IU`, `ISF`, `ISR`, ...), `strand_mapping_bias` and
  `compatible_fragment_ratio` (from salmon's `lib_format_counts.json`), `salmon_warnings`
  (salmon's diagnostics, e.g. a strand bias in an unstranded library; `;`-separated), then
  sample sheet and allow-listed SPARC columns.
- `samples_full.csv`: `sample_id` and every metadata column.
- `qc/<sample_id>.fastp.json`, `qc/<sample_id>.fastp.html` (fastp's reports; the HTML loads
  plotly from `cdn.plot.ly` / `opengene.org` when opened, so its charts need a network
  connection in the browser; the JSON has the same numbers),
  `qc/<sample_id>.salmon.json` (salmon's `meta_info.json` and `lib_format_counts.json`).
- `summary.json`: samples, genes, median mapping rate, layouts, library types, index
  (organism, source, decoys, salmon version and index format), per sample reads, mapping rate,
  library type, strand bias, compatible fragment ratio, salmon's warnings, genes detected,
  assigned counts; `warnings` (each also logged as a WARNING: salmon's diagnostics per sample;
  a mapping rate below 0.8 x the run's median; library types that differ between samples;
  single- and paired-end samples in one run; files split from a fake pair); `skipped`,
  `ignored`, and `aborted` when the run was stopped.
- `provenance.json`: versions (salmon, fastp, Python packages), parameters, inputs found
  (incl. the index's `index.json`), per-sample status and timings, peak memory, `citations`.

`transcripts/` (an addition to the conventions, see below): `<sample_id>/quant.sf`,
`cmd_info.json`, `lib_format_counts.json`, `aux_info/meta_info.json` — the layout tximport,
tximeta and fishpond read, for transcript-level or isoform analyses; `provenance.json`.

## Runtime

`computeTypes: [standard]`, 8 vCPU / 16 GiB, 12 h timeout. One sample at a time; per sample the
work folder in `/tmp` (on the platform the task's ephemeral disk, 20 GiB by default on Fargate)
holds the lanes concatenated (only when there are several), the trimmed reads (gzip level 1:
40 MB per million 50 bp single-end reads) and salmon's intermediate mapping file (RAD, 133 MB per
million such reads, measured), and is removed afterwards: about 1.7 GB per 10 M single-end 50 bp
reads, more for longer or paired-end reads. A 50 M-pair 2x150 bp library needs roughly 20 GB:
raise the ephemeral storage for deep libraries.
Memory is salmon's: the index (1.3 GB on disk for GENCODE M39 transcripts) plus the mapping and
inference state, independent of the number of samples.

Measured with `docker run --platform linux/amd64 --cpus 4` (amd64 emulation on an M-series Mac;
native amd64 is faster) against the transcripts-only GENCODE M39 index, timings and peaks from
`provenance.json`:

| Data | wall time | per sample | peak memory |
|---|---|---|---|
| Discover 427, 12 libraries subsampled to 1 M single-end 50 bp reads (`--memory 6g`) | **337 s** | 23-32 s (fastp 2 s, salmon 22-30 s) | salmon 2.6 GB RSS, Python 0.6 GB; cgroup 4.3 GB incl. page cache |
| Discover 427, full depth (24-29 M reads per library; 11 libraries in one run, `--memory 6g`) | **3,273 s** (55 min) | 4-5 min (fastp 25-31 s, salmon 220-280 s; FT2 385 s) | salmon 3.6 GB RSS, Python 0.6 GB; cgroup 6 GB (the limit, page cache) |
| Synthetic fixture (3 samples x 3,000 reads) | 3 s | | |

In the full-depth run two libraries (T2, and FT2 in a retry at `--memory 5g`) were killed by
the out-of-memory killer while other builders' containers shared the 8 GB Docker VM; T2 and FT2
then ran together on their own at `--memory 6g` (605 s, salmon peak 3.6 GB, cgroup 4.4 GB). Give
the container at least twice salmon's peak (the trimmed reads and salmon's intermediate mapping
file are written while it runs); 16 GiB is ample for a transcripts-only index.

So a typical bulk library (25 M reads) takes about 5 minutes on 4 emulated vCPUs; 100 such
samples fit in the 12 h timeout on 8 vCPUs. A decoy-aware index is much larger (2.2 billion
k-mers for mouse vs 179 M, see rnaseq-index) and was not measured here: it raises salmon's
memory, which is why the default is 16 GiB rather than 8. Image: 344 MB (fastp built from source, statically linked).

## Results on the demo data

Discover 427 (RE-JOIN; 12 pooled single-end 50 bp libraries of mouse trigeminal ganglion
neurons, female/male x naive/tongue tumour), GENCODE M39 transcripts-only index from
rnaseq-index, defaults. Outputs in the workflow repository: `data/work/tg-tumour/quant/`
(1 M-read subsample, all 12 libraries: the input for bulk-de), `quant-full-11-samples/` and
`quant-full-T2-FT2/` (full depth), `quant-vs-authors/` (the comparisons below and the script).

- **Every library is unstranded (`U`)**, single-end, 87-89 % of reads mapped, except **FT2:
  64 % mapped**, with salmon's warning of a possible strand bias (0.47) in an unstranded
  library. The authors' own counts for FT2 are also the lowest (4.1 M vs 6.4-8.0 M), so FT2 is
  the weakest library in both pipelines; bulk-de should look at it in the PCA.
- fastp removed almost nothing (2,000-4,000 of ~25 M reads per library): the reads were already
  adapter-trimmed and of high quality.
- 77,296 genes in the table (78,374 in the index: 1,078 genes have only transcripts whose
  sequences duplicate another gene's, which salmon collapses); at full depth 30,000-34,000 genes
  have a count > 0 per library.
- **tximport check:** `tximport(type = "salmon", tx2gene, countsFromAbundance =
  "lengthScaledTPM")` (tximport 1.40.0, R 4.6.1) on the same `quant.sf` files gives the same
  rounded counts in all 927,552 gene x sample cells (largest unrounded difference 0.49999) and the
  same TPM up to the 4 decimals written.
- **Against the authors' counts** (`<id>_Genes_ReadCount.txt`: TopHat2 on UCSC mm9 and HTSeq, per
  Ibrahim et al. 2023; 22,614 gene symbols), joined on gene name (our counts summed per
  `gene_name`; 18,253 names in common), Spearman correlation per library:

  | Reads | all common genes | genes with a count in either | genes with ≥ 10 authors' counts |
  |---|---|---|---|
  | Full depth (24-29 M reads; T2 and FT2 from their own run) | median **0.957** (0.945-0.959; FT2 0.910) | median 0.940 (FT2 0.860) | median 0.931 (FT2 0.831) |
  | 1 M-read subsample | median 0.931 (0.913-0.933; FT2 0.877) | median 0.903 | median 0.891 |

  Our totals over the common genes are ~2.5x the authors' (e.g. FN1 17.6 M vs 7.2 M): HTSeq on
  mm9 drops multi-mapping and ambiguous reads and has an older, smaller annotation, while salmon
  divides such reads between transcripts. The agreement is what two different reference builds
and counting models (alignment to mm9 + HTSeq vs salmon on GENCODE M39) give, and it does not
improve when low-count genes are left out (≥ 10 authors' counts: median 0.931). Each library
matches its own authors' counts best (full depth, 11 x 11 comparison), though only by a small
margin (0.957 vs ≤ 0.952 for another library of the same kind), as expected for replicate
pools of the same tissue.

## Limitations

- Counts are estimates, not read tallies: reads shared between genes are divided by salmon's
  likelihood model, so they do not match featureCounts/HTSeq counts exactly (see Results).
- Gene-level counts hide isoform switches; use `transcripts/` for those.
- Only what is in the index's annotation can be counted; novel genes and isoforms cannot.
- `--gcBias` needs fragment lengths, so it is applied to paired-end data only; single-end
  libraries get sequence-bias correction only.
- No UMI deduplication or rRNA removal: the app is for standard poly-A or ribo-depleted bulk
  libraries. 3'-tag libraries (QuantSeq, BRB-seq) need different counting.
- Library type is detected per sample; a mix of strandedness across samples is warned about in
  the log and `summary.json` but not acted on.
- Single- and paired-end libraries in one run are quantified alike, but their counts differ
  systematically (fragment length and GC correction, mappability): include `layout` (in
  `samples.csv`) as a covariate in the DE design, or analyse them separately.

## Deviations from the brief

- `--validateMappings` is not passed: since salmon 2.0, selective alignment (what that flag
  turned on in salmon ≤ 1.x) is the default and only mapping mode, and salmon 2.8 accepts the
  flag only to ignore it with a warning (salmon `MIGRATION.md`).
- A second output port, `transcripts`, keeps salmon's per-sample `quant.sf` (the conventions
  list only `counts`); throwing away transcript-level estimates would prevent isoform
  analyses and re-aggregation with tximport/tximeta in R.

## Licences

Our code: Apache-2.0 (`LICENSE`). salmon 2.x is BSD-3-Clause (C++ salmon ≤ 1.12 was GPL-3.0),
fastp MIT (built with ISA-L BSD-3-Clause, libdeflate MIT, Highway Apache-2.0/BSD-3-Clause),
pandas / NumPy BSD-3-Clause, openpyxl MIT, defusedxml PSF-2.0; nothing in the image
is GPL. The lengthScaledTPM formula is implemented from the paper and tximport's (LGPL) source;
no tximport code is included. See `THIRD_PARTY.md`.

## Credits

- **salmon** (Rob Patro and the COMBINE-lab; Patro, Duggal, Love, Irizarry, Kingsford)
  quantifies each sample's transcripts.
- **fastp** (Shifu Chen and OpenGene; Chen, Zhou, Chen, Gu) trims reads and writes the QC
  reports; built here from source with Intel's **ISA-L**, Eric Biggers' **libdeflate** and
  Google's **Highway**.
- **Selective alignment with decoys** (Srivastava, Malik, Sarkar, Zakeri, Almodaresi, Soneson,
  Love, Kingsford, Patro) is cited when the index has decoys.
- **tximport** (Soneson, Love, Robinson) defines the lengthScaledTPM gene-level summary this
  app re-implements.
- **GENCODE** (the GENCODE consortium; Mudge, Frankish et al.) provides the demo annotation,
  mouse release M39.
- **Discover dataset 427** (Tumanov, Shein, Lai, Chen, Ruparel; RE-JOIN, CC BY 4.0), described
  by Ibrahim et al. 2023, is the demo data.

## Citations

Papers checked against Crossref (title, authors and their order, journal, year, volume, issue,
pages) and PubMed (E-utilities `esearch` on the DOI, then `esummary`); the dataset against
DataCite and the Discover API. Full author lists are in `CITATION.cff`.

- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware
  quantification of transcript expression. *Nat Methods* 14(4):417-419 (2017).
  doi:10.1038/nmeth.4197. PMID 28263959.
- Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor.
  *Bioinformatics* 34(17):i884-i890 (2018). doi:10.1093/bioinformatics/bty560. PMID 30423086.
- Chen S. fastp 1.0: An ultra-fast all-round tool for FASTQ data quality control and
  preprocessing. *iMeta* 4(5):e70078 (2025). doi:10.1002/imt2.70078. PMID 41112039. (the
  fastp 1.x release used here; the 2023 iMeta paper on deduplication is not cited because
  deduplication is not used)
- Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level
  estimates improve gene-level inferences. *F1000Research* 4:1521 (2015).
  doi:10.12688/f1000research.7563.2. PMID 26925227. (version 2, the one PubMed indexes;
  Crossref dates version 1, doi:10.12688/f1000research.7563.1, 2015-12-30 and version 2
  2016-02-29)
- Srivastava A, Malik L, Sarkar H, Zakeri M, Almodaresi F, Soneson C, Love MI, Kingsford C,
  Patro R. Alignment and mapping methodology influence transcript abundance estimation.
  *Genome Biol* 21(1):239 (2020). doi:10.1186/s13059-020-02151-8. PMID 32894187. (cited in
  provenance when the index has decoys)
- Mudge JM, Carbonell-Sala S, Diekhans M, et al. (56 authors). GENCODE 2025: reference gene
  annotation for human and mouse. *Nucleic Acids Res* 53(D1):D966-D975 (2025).
  doi:10.1093/nar/gkae1078. PMID 39565199. (cited in provenance for a GENCODE index)
- Conesa A, Madrigal P, Tarazona S, Gomez-Cabrero D, Cervera A, McPherson A, Szcześniak MW,
  Gaffney DJ, Elo LL, Zhang X, Mortazavi A. A survey of best practices for RNA-seq data
  analysis. *Genome Biol* 17:13 (2016). doi:10.1186/s13059-016-0881-8. PMID 26813401.
  (practice followed; not a dependency)
- Love MI, Anders S, Kim V, Huber W. RNA-Seq workflow: gene-level exploratory analysis and
  differential expression. *F1000Research* 4:1070 (2015). doi:10.12688/f1000research.7035.1.
  PMID 26674615. (practice followed; not a dependency)

Demo data: Tumanov AV, Shein SA, Lai Z, Chen Y, Ruparel SB. Sex-dependent differences in the
genomic profile of lingual sensory neurons in naive and tongue-tumor bearing mice. SPARC
Portal, Pennsieve Discover dataset 427 v1 (2025). doi:10.26275/gaap-zlpj. CC BY 4.0.
Described by: Ibrahim T, Wu P, Wang LJ, Fang-Mei C, Murillo J, Merlo J, Shein SS, Tumanov AV,
Lai Z, Weldon K, Chen Y, Ruparel S. Sex-dependent differences in the genomic profile of
lingual sensory neurons in naïve and tongue-tumor bearing mice. *Sci Rep* 13(1):13117 (2023).
doi:10.1038/s41598-023-40380-6. PMID 37573456.
