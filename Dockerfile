# rnaseq-quant: fastp + salmon (official release binaries) + pandas.
# linux/amd64; no network at run time; runs as uid 1000.
FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337 AS fetch

# salmon 2.8.0 (BSD-3-Clause), the release's x86_64 Linux build, checked against the sha256 the
# release publishes (salmon-cli-x86_64-unknown-linux-gnu.tar.xz.sha256).
ARG SALMON_VERSION=2.8.0
ARG SALMON_SHA256=10a7f7994df7ef0592ac191cc79a1abd41316463000b307647a9dcab21e1efa6
ADD --checksum=sha256:${SALMON_SHA256} \
    https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-cli-x86_64-unknown-linux-gnu.tar.xz /tmp/salmon.tar.xz
# fastp 1.4.0 (MIT), the statically linked Linux build the fastp README points to. Upstream
# publishes no checksum; this sha256 was recorded when the file was first downloaded (2026-10-31).
ARG FASTP_VERSION=1.4.0
ARG FASTP_SHA256=e56c3218b76cc1d6d8a7a32315fa6a41407f5d271fe561fc62cee15941190844
ADD --checksum=sha256:${FASTP_SHA256} https://opengene.org/fastp/fastp.${FASTP_VERSION} /tmp/fastp
ADD --checksum=sha256:52e18866a994ee9839aa6d15696bc1fe425e6466631d3a57a75c089dd84ae88e \
    https://raw.githubusercontent.com/OpenGene/fastp/v${FASTP_VERSION}/LICENSE /tmp/fastp.LICENSE
RUN apt-get update && apt-get install -y --no-install-recommends xz-utils && rm -rf /var/lib/apt/lists/* \
    && tar -xJf /tmp/salmon.tar.xz -C /tmp \
    && install -m 0755 /tmp/salmon-cli-x86_64-unknown-linux-gnu/salmon /usr/local/bin/salmon \
    && install -D -m 0644 /tmp/salmon-cli-x86_64-unknown-linux-gnu/LICENSE /usr/local/share/licenses/salmon/LICENSE \
    && install -m 0755 /tmp/fastp /usr/local/bin/fastp \
    && install -D -m 0644 /tmp/fastp.LICENSE /usr/local/share/licenses/fastp/LICENSE

FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

COPY --from=fetch /usr/local/bin/salmon /usr/local/bin/fastp /usr/local/bin/
COPY --from=fetch /usr/local/share/licenses /usr/local/share/licenses
RUN salmon --version && fastp --version

WORKDIR /processor
# Pinned Python dependencies (binary wheels only).
COPY requirements.txt /processor/requirements.txt
RUN pip install --no-cache-dir --only-binary=:all: -r /processor/requirements.txt \
    && find /usr/local/lib/python3.12 -depth -type d -name __pycache__ -exec rm -rf {} + \
    && find /usr/local/lib/python3.12/site-packages -depth -type d \( -name tests -o -name test \) -exec rm -rf {} +

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH=/processor \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    HOME=/tmp

# The platform runs the task as uid 1000; nothing is written outside OUTPUT_DIR and /tmp.
USER 1000

ENTRYPOINT ["python3", "-m", "processor.main"]
