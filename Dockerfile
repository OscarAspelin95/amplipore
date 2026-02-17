FROM rust:1.90.0-bullseye AS builder

ENV SINTAX_RS_VERSION="v0.0.2" \
	FASTQ_RS_VERSION="v0.0.2"

# Compile Rust dependencies from source with SIMD acceleration.
RUN apt-get update && apt-get install -y git \
	# sintax_rs
	&& git clone --branch ${SINTAX_RS_VERSION} --depth 1 https://github.com/OscarAspelin95/sintax_rs.git /usr/src/sintax_rs \
	&& cd /usr/src/sintax_rs \
	&& RUSTFLAGS="-C target-cpu=native" cargo build --release \
	&& cp target/release/sintax_rs /usr/local/bin/sintax_rs \
	# fastq_rs
	&& git clone --branch ${FASTQ_RS_VERSION} --depth 1 https://github.com/OscarAspelin95/fastq_rs.git /usr/src/fastq_rs \
	&& cd /usr/src/fastq_rs \
	&& RUSTFLAGS="-C target-cpu=native" cargo build --release \
	&& cp target/release/fastq_rs /usr/local/bin/fastq_rs \
	&& rm -rf /usr/src/sintax_rs /usr/src/fastq_rs


FROM python:3.13.7-bookworm

ENV USEARCH_VERSION="12.0-beta" \
	DASH_PORT="8000" \
	BLAST_VERSION="2.17.0"

WORKDIR /usr/src/app

EXPOSE ${DASH_PORT}

COPY --from=builder /usr/local/bin/sintax_rs /usr/local/bin/sintax_rs
COPY --from=builder /usr/local/bin/fastq_rs /usr/local/bin/fastq_rs

COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

ENV	DB_PATH="/usr/src/db/db.fasta" \
	UV_COMPILE_BYTECODE=1 \
	UV_LINK_MODE=copy \
	PATH="/usr/src/app/.venv/bin:$PATH"

COPY ./pyproject.toml ./uv.lock ./
RUN uv sync --frozen --no-dev --no-install-project \
	&& apt-get update && apt-get install -y curl wget \
	&& mkdir -p /usr/src/deps \
	# Usearch
	&& cd /usr/src/deps \
	&& wget https://github.com/rcedgar/usearch12/releases/download/v"${USEARCH_VERSION}"1/usearch_linux_x86_"${USEARCH_VERSION}" \
	&& chmod +x ./usearch_linux_x86_"${USEARCH_VERSION}" && mv ./usearch_linux_x86_"${USEARCH_VERSION}" /usr/local/bin/usearch \
	# BLAST
	&& wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/"${BLAST_VERSION}"/ncbi-blast-"${BLAST_VERSION}"+-x64-linux.tar.gz \
	&& tar -xf ncbi-blast-"${BLAST_VERSION}"+-x64-linux.tar.gz \
	&& cp ncbi-blast-"${BLAST_VERSION}"+/bin/* /usr/local/bin/ \
	# Remove dependencies
	&& rm -r /usr/src/deps


COPY ./app /usr/src/app/

# Download database
RUN mkdir -p /usr/src/db \
	&& mkdir -p /tmp/db \
	&& uv run python3 database.py --outdir /tmp/db \
	&& cp /tmp/db/db.fasta ${DB_PATH} \
	&& rm -r /tmp/db
