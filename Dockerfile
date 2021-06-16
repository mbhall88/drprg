FROM rust:1.52 AS builder

COPY . /drprg

WORKDIR /drprg

ARG MAFFT_URL="https://mafft.cbrc.jp/alignment/software/mafft_7.475-1_amd64.deb"
ARG DEB="/mafft.deb"

RUN cargo build --release \
    && make EXTDIR=/usr/bin pandora \
    && make EXTDIR=/usr/bin makeprg \
    && make EXTDIR=/usr/bin bcftools \
    && wget -O "$DEB" "$MAFFT_URL"


FROM ubuntu:focal

ARG DEB="/mafft.deb"

COPY --from=builder /drprg/target/release/drprg /usr/bin/pandora /usr/bin/make_prg /usr/bin/
COPY --from=builder "$DEB" "$DEB"

ARG BCFTOOLS_URL="https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2"
ARG PKGS="ca-certificates wget autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libperl-dev libgsl0-dev"

RUN apt update \
    && apt-get install --no-install-recommends -y $PKGS \
    && update-ca-certificates -f \
    && dpkg -i "$DEB" \
    && rm -f "$DEB"

WORKDIR /bcftools
RUN ( wget -O - "$BCFTOOLS_URL" | tar -xjf - ) \
	&& cd bcftools-1.12 \
	&& ./configure --prefix=/usr \
	&& make \
	&& make install \
	&& apt remove -y $PKGS \
	&& rm -rf /var/lib/apt/lists/* \
	&& rm -rf /bcftools

RUN drprg --help && pandora --help && make_prg --help && mafft --version && bcftools --version