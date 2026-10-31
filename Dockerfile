# rnaseq-quant: fastp (built from source) + salmon (official release binary) + pandas.
# linux/amd64; no network at run time; runs as uid 1000.
FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337 AS fetch

# salmon 2.8.0 (BSD-3-Clause), the release's x86_64 Linux build, checked against the sha256 the
# release publishes (salmon-cli-x86_64-unknown-linux-gnu.tar.xz.sha256).
ARG SALMON_VERSION=2.8.0
ARG SALMON_SHA256=10a7f7994df7ef0592ac191cc79a1abd41316463000b307647a9dcab21e1efa6
ADD --checksum=sha256:${SALMON_SHA256} \
    https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-cli-x86_64-unknown-linux-gnu.tar.xz /tmp/salmon.tar.xz
RUN apt-get update && apt-get install -y --no-install-recommends xz-utils && rm -rf /var/lib/apt/lists/* \
    && tar -xJf /tmp/salmon.tar.xz -C /tmp \
    && install -m 0755 /tmp/salmon-cli-x86_64-unknown-linux-gnu/salmon /usr/local/bin/salmon \
    && install -D -m 0644 /tmp/salmon-cli-x86_64-unknown-linux-gnu/LICENSE /usr/local/share/licenses/salmon/LICENSE

# fastp 1.4.0 (MIT) built from source as a static binary, with the libraries its README lists
# (ISA-L, libdeflate, Google Highway), each fetched by commit SHA and checked against the sha256
# of GitHub's source archive. The build toolchain stays in this stage.
FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337 AS fastp
RUN apt-get update && apt-get install -y --no-install-recommends build-essential cmake nasm autoconf automake \
        libtool pkg-config && rm -rf /var/lib/apt/lists/*
# fastp v1.4.0, ISA-L v2.31.0, libdeflate v1.26, Highway 1.3.0
ARG FASTP_COMMIT=ee125a52fddae9d8189e13cd8a33ad04988a1351
ARG ISAL_COMMIT=bd226375027899087bd48f3e59b910430615cc0a
ARG LIBDEFLATE_COMMIT=92e6a0db9fa848d742f9eb286c92afc60f2c3dda
ARG HWY_COMMIT=ac0d5d297b13ab1b89f48484fc7911082d76a93f
ADD --checksum=sha256:fe59a5f585a3a7a944a5462a882dac7460f94dab3f648d52640a0d33662d6ea0 \
    https://github.com/OpenGene/fastp/archive/${FASTP_COMMIT}.tar.gz /src/fastp.tar.gz
ADD --checksum=sha256:41ac176bb363d52c626f9061e043c7e0b80dcf0ac0a65e91d09506e38ff49652 \
    https://github.com/intel/isa-l/archive/${ISAL_COMMIT}.tar.gz /src/isa-l.tar.gz
ADD --checksum=sha256:af56498149827f5a5e3f1ffb653c79c8c75557f3dc2caf6002cac38a59b844b7 \
    https://github.com/ebiggers/libdeflate/archive/${LIBDEFLATE_COMMIT}.tar.gz /src/libdeflate.tar.gz
ADD --checksum=sha256:a7a816f4b62a0414ff0d39c0a8875847468dd6f7c9ad71781a24469a756cdea7 \
    https://github.com/google/highway/archive/${HWY_COMMIT}.tar.gz /src/highway.tar.gz
WORKDIR /src
RUN for f in *.tar.gz; do tar -xzf "$f"; done \
    && cd /src/isa-l-${ISAL_COMMIT} && ./autogen.sh && ./configure --prefix=/usr/local --libdir=/usr/local/lib \
        --enable-static --disable-shared && make -j"$(nproc)" && make install \
    && cd /src/libdeflate-${LIBDEFLATE_COMMIT} && cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_LIBDIR=lib \
        -DLIBDEFLATE_BUILD_SHARED_LIB=OFF -DLIBDEFLATE_BUILD_GZIP=OFF && cmake --build build -j"$(nproc)" \
        && cmake --install build \
    && cd /src/highway-${HWY_COMMIT} && cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_LIBDIR=lib \
        -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF -DHWY_ENABLE_TESTS=OFF -DHWY_ENABLE_EXAMPLES=OFF \
        -DHWY_ENABLE_CONTRIB=OFF && cmake --build build -j"$(nproc)" && cmake --install build \
    && cd /src/fastp-${FASTP_COMMIT} && make -j"$(nproc)" LIBRARY_DIRS=/usr/local/lib \
    && ./fastp --version 2>&1 | grep -q "fastp 1.4.0" \
    && ! ldd ./fastp \
    && install -m 0755 fastp /usr/local/bin/fastp \
    && install -D -m 0644 LICENSE /usr/local/share/licenses/fastp/LICENSE \
    && install -D -m 0644 /src/isa-l-${ISAL_COMMIT}/LICENSE /usr/local/share/licenses/isa-l/LICENSE \
    && install -D -m 0644 /src/libdeflate-${LIBDEFLATE_COMMIT}/COPYING /usr/local/share/licenses/libdeflate/COPYING \
    && install -D -m 0644 /src/highway-${HWY_COMMIT}/LICENSE /usr/local/share/licenses/highway/LICENSE

FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

COPY --from=fetch /usr/local/bin/salmon /usr/local/bin/
COPY --from=fastp /usr/local/bin/fastp /usr/local/bin/
COPY --from=fetch /usr/local/share/licenses /usr/local/share/licenses
COPY --from=fastp /usr/local/share/licenses /usr/local/share/licenses
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
