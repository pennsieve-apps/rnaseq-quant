# Third-party software and data

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components; their licences apply to the image. Nothing in it is GPL.

| Component | What | Licence | Source |
|---|---|---|---|
| salmon 2.8.0 | `salmon quant`: selective alignment and transcript abundance estimation. The release's prebuilt `salmon-cli-x86_64-unknown-linux-gnu.tar.xz`, checked against the sha256 published with the release (`10a7f799…1efa6`). Licence file at `/usr/local/share/licenses/salmon/LICENSE`. salmon 2.x is a Rust rewrite; the binary statically includes Rust crates under MIT / Apache-2.0 / BSD licences (listed in salmon's `Cargo.lock` at tag `v2.8.0`). | BSD-3-Clause (salmon ≥ 2.0; C++ salmon ≤ 1.12 was GPL-3.0) | https://github.com/COMBINE-lab/salmon/releases/tag/v2.8.0 |
| fastp 1.4.0 | Adapter and quality trimming, QC report. The statically linked Linux build the fastp README links (`https://opengene.org/fastp/fastp.1.4.0`); upstream publishes no checksum, so the sha256 recorded at first download (`e56c3218…90844`) is pinned in the `Dockerfile`. Its licence file (tag `v1.4.0`) is at `/usr/local/share/licenses/fastp/LICENSE`. The static binary includes the libraries fastp's README lists as dependencies: ISA-L (BSD-3-Clause), libdeflate (MIT) and Google Highway (Apache-2.0 or BSD-3-Clause). | MIT | https://github.com/OpenGene/fastp |
| pandas 3.1.0 | Tables, gene-level summing. | BSD-3-Clause | https://pandas.pydata.org |
| NumPy 2.5.4 | Arrays (wheel bundles OpenBLAS, BSD-3-Clause, and the GCC Fortran runtime under the GCC Runtime Library Exception). | BSD-3-Clause | https://numpy.org |
| openpyxl 3.1.5, et_xmlfile 2.0.0 | Reading SPARC `samples.xlsx` / `subjects.xlsx`. | MIT | https://foss.heptapod.net/openpyxl/openpyxl |
| defusedxml 0.7.1 | openpyxl parses the (untrusted) xlsx XML through it, which blocks entity-expansion and external-entity attacks. | PSF-2.0 | https://github.com/tiran/defusedxml |
| python-dateutil 3.0.0, six 1.17.0, tzdata 2026.6 | pandas dependencies. | Apache-2.0 / BSD-3-Clause; MIT; Apache-2.0 | PyPI |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12), pinned by digest. | PSF-2.0 (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

The gene-level summarisation re-implements the published lengthScaledTPM formula of tximport
(LGPL ≥ 2); no tximport code is included.

## Data used in testing (not in the image)

| Data | Licence | Source |
|---|---|---|
| Discover dataset 427 v1 (Tumanov, Shein, Lai, Chen, Ruparel; RE-JOIN): 12 single-end mouse trigeminal ganglion libraries with the authors' gene counts | CC BY 4.0 | https://doi.org/10.26275/gaap-zlpj |
| GENCODE mouse release M39 (index built by rnaseq-index) | Public data; EMBL-EBI terms of use (https://www.ebi.ac.uk/about/terms-of-use); cite GENCODE 2025 | https://www.gencodegenes.org/mouse/release_M39.html |
