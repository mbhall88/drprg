FROM rust:1.51 AS builder

COPY . /drprg

WORKDIR /drprg

ARG MAFFT_URL="https://mafft.cbrc.jp/alignment/software/mafft_7.475-1_amd64.deb"
ARG DEB="/mafft.deb"

RUN cargo build --release \
    && make EXTDIR=/usr/bin pandora \
    && make EXTDIR=/usr/bin makeprg \
    && wget -O "$DEB" "$MAFFT_URL"


FROM ubuntu:bionic

ARG DEB="/mafft.deb"

COPY --from=builder /drprg/target/release/drprg /usr/bin/pandora /usr/bin/make_prg /usr/bin/
COPY --from=builder "$DEB" "$DEB"

RUN dpkg -i "$DEB" \
    && rm -f "$DEB"

RUN drprg --help && pandora --help && make_prg --help && mafft --version