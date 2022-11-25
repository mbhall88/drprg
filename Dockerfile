FROM rust:1.65 AS builder

COPY . /drprg

ARG JUST_URL="https://just.systems/install.sh"
ARG V_JUST="1.2.0"
ARG MAFFT_URL="https://mafft.cbrc.jp/alignment/software/mafft_7.505-1_amd64.deb"
ARG DEB="/mafft.deb"

RUN apt update \
    && apt install -y cmake \
    && cargo install --locked --path /drprg --root /usr \
    && (wget -O - "$JUST_URL" | bash -s -- --to /usr/bin --tag "$V_JUST") \
    && just -f /drprg/justfile EXTDIR=/usr/bin pandora \
    && just -f /drprg/justfile EXTDIR=/usr/bin bcftools \
    && just -f /drprg/justfile EXTDIR=/usr/bin makeprg \
    && wget -O "$DEB" "$MAFFT_URL" \
    && dpkg -i "$DEB" \
    && apt remove -y cmake \
    && rm -rf /drprg "$DEB" /var/lib/apt/lists/*

RUN drprg --help && pandora --help && make_prg --help && mafft --version && bcftools --version

# These two are required for the paper pipeline
RUN cargo install --root /usr --version 0.7.0 rasusa && rasusa --version

ARG SEQFU_URL="https://github.com/telatin/seqfu2/releases/download/v1.16.0/SeqFu-v1.16.0-Linux-x86_64.zip"
RUN wget -qO /seqfu.zip $SEQFU_URL && unzip -d /usr/bin /seqfu.zip && rm /seqfu.zip && seqfu version
