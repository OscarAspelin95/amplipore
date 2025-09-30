FROM rust:1.90.0-bullseye AS builder

ENV SINTAX_RS_VERSION="0.0.1-alpha" \
    FASTQ_RS_VERSION="0.0.2-alpha"

# Compile Rust dependencies.
RUN apt-get update && apt-get install -y wget \
    && mkdir -p /usr/src/deps \
    # sintax_rs
    && cd /usr/src/deps \
    && wget https://github.com/OscarAspelin95/sintax_rs/archive/refs/tags/v"${SINTAX_RS_VERSION}".tar.gz \
    && tar -xf v"${SINTAX_RS_VERSION}".tar.gz && cd sintax_rs-"${SINTAX_RS_VERSION}" \
    && cargo build --release && cp ./target/release/sintax_rs /usr/local/bin/sintax_rs \
    # fastq_rs
    && cd /usr/src/deps \
    && wget https://github.com/OscarAspelin95/fastq_rs/archive/refs/tags/v"${FASTQ_RS_VERSION}".tar.gz \
    && tar -xf v"${FASTQ_RS_VERSION}".tar.gz && cd fastq_rs-"${FASTQ_RS_VERSION}" \
    && cargo build --release && cp ./target/release/fastq_rs /usr/local/bin/fastq_rs


FROM python:3.13.7-bookworm

ENV USEARCH_VERSION="12.0-beta" \
    DASH_PORT="8000"

WORKDIR /usr/src/app

EXPOSE ${DASH_PORT}

COPY --from=builder /usr/local/bin/sintax_rs /usr/local/bin/sintax_rs
COPY --from=builder /usr/local/bin/fastq_rs /usr/local/bin/fastq_rs

COPY ./requirements.txt requirements.txt
RUN apt-get update && apt-get install -y curl && pip install -r requirements.txt \
    && wget https://github.com/rcedgar/usearch12/releases/download/v"${USEARCH_VERSION}"1/usearch_linux_x86_"${USEARCH_VERSION}" \
    && chmod +x ./usearch_linux_x86_"${USEARCH_VERSION}" && mv ./usearch_linux_x86_"${USEARCH_VERSION}" /usr/local/bin/usearch