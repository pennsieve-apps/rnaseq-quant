# Third-party software and data

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components; their licences apply to the image. Nothing in it is GPL apart from
the C/C++ runtimes linked into fastp (see their row).

| Component | What | Licence | Source |
|---|---|---|---|
| salmon 2.8.0 | `salmon quant`: selective alignment and transcript abundance estimation. The release's prebuilt `salmon-cli-x86_64-unknown-linux-gnu.tar.xz`, checked against the sha256 published with the release (`10a7f799…1efa6`). Licence file at `/usr/local/share/licenses/salmon/LICENSE`. salmon 2.x is a Rust rewrite; the binary statically includes Rust crates under MIT / Apache-2.0 / BSD licences (listed in salmon's `Cargo.lock` at tag `v2.8.0`). | BSD-3-Clause (salmon ≥ 2.0; C++ salmon ≤ 1.12 was GPL-3.0) | https://github.com/COMBINE-lab/salmon/releases/tag/v2.8.0 |
| fastp 1.4.0 | Adapter and quality trimming, QC report. Built from source in a builder stage of the `Dockerfile` as a static binary: GitHub archive of tag `v1.4.0` by commit `ee125a52…1351`, checked against its sha256. Licence file at `/usr/local/share/licenses/fastp/LICENSE`. | MIT | https://github.com/OpenGene/fastp |
| ISA-L v2.31.0 | Fast gzip in fastp (statically linked); commit `bd226375…cc0a`, archive sha256 checked. `/usr/local/share/licenses/isa-l/LICENSE`. | BSD-3-Clause | https://github.com/intel/isa-l |
| libdeflate v1.26 | Fast deflate in fastp (statically linked); commit `92e6a0db…3dda`, archive sha256 checked. `/usr/local/share/licenses/libdeflate/COPYING`. | MIT | https://github.com/ebiggers/libdeflate |
| Google Highway 1.3.0 | SIMD in fastp (statically linked); commit `ac0d5d29…a93f`, archive sha256 checked. `/usr/local/share/licenses/highway/LICENSE`. | Apache-2.0 or BSD-3-Clause | https://github.com/google/highway |
| glibc / libstdc++ (Debian 12) | Statically linked into fastp by the build. | glibc LGPL-2.1+ (static linking: the `Dockerfile` and the pinned sources let anyone relink fastp against another glibc, which is what the LGPL asks); libstdc++ GPL-3.0 with the GCC Runtime Library Exception (no copyleft on the linked program) | https://www.debian.org |
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
