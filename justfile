PROJECT := "drprg"
EXTDIR := absolute_path("src/ext")

make_extdir:
    mkdir -p {{ EXTDIR }}

MAKE_PRG_VERSION := "0.4.0"
MAKEPRG_URL := "https://github.com/iqbal-lab-org/make_prg/releases/download/" + MAKE_PRG_VERSION + "/make_prg_" + MAKE_PRG_VERSION
MAKEPRG := join(EXTDIR, "make_prg")

# download make_prg binary
makeprg: make_extdir
    wget "{{ MAKEPRG_URL }}" -O "{{ MAKEPRG }}"
    chmod +x "{{ MAKEPRG }}"

PANDORA_VERSION := "0.10.0-alpha.0.1"
PANDORA_URL := "https://github.com/mbhall88/pandora/releases/download/" + PANDORA_VERSION + "/pandora-linux-precompiled-v" + PANDORA_VERSION
PANDORA := join(EXTDIR, "pandora")

# download pandora binary
pandora: make_extdir
    wget "{{ PANDORA_URL }}" -O "{{ PANDORA }}"
    chmod +x "{{ PANDORA }}"

BCFTOOLS_VERSION := "1.15.1"
BCFTOOLS_URL := "https://github.com/samtools/bcftools/releases/download/" + BCFTOOLS_VERSION + "/bcftools-" + BCFTOOLS_VERSION + ".tar.bz2"
BCFTOOLS := join(EXTDIR, "bcftools")

# download and build bcftools binary
bcftools: make_extdir
    #!/usr/bin/env bash
    set -euxo pipefail
    cd {{ EXTDIR }}
    wget {{BCFTOOLS_URL}} -O - | tar -xjf -
    cd bcftools-{{ BCFTOOLS_VERSION }}
    ./configure --prefix=$PWD
    make
    make install
    cd {{ invocation_directory() }}
    ln -sf bcftools-{{ BCFTOOLS_VERSION }}/bin/bcftools {{ BCFTOOLS }}

MAFFT_VERSION := "7.505"
MAFFT_URL := "https://mafft.cbrc.jp/alignment/software/mafft-" + MAFFT_VERSION + "-without-extensions-src.tgz"
MAFFT_SRC := join(EXTDIR, "mafft")
MAFFT_BIN := join(MAFFT_SRC, "bin/mafft")

# download and build mafft
mafft: make_extdir
    #!/usr/bin/env bash
    set -euxo pipefail
    mkdir -p {{MAFFT_SRC}}
    wget -O - {{MAFFT_URL}} | tar xzf - -C {{MAFFT_SRC}} --strip-components=1
    sed -i "1s?/usr/local?{{MAFFT_SRC}}?" {{MAFFT_SRC}}/core/Makefile
    make -C {{MAFFT_SRC}}/core install PREFIX={{MAFFT_SRC}}

# download and build all dependencies
deps: makeprg pandora bcftools mafft

# run clippy to check for linting issues
lint:
    cargo clippy --all-features --all-targets -- -D warnings

# run all tests
test:
    cargo test -v --all-targets --no-fail-fast

# get coverage with tarpaulin
coverage:
    cargo tarpaulin -t 300 -- --test-threads 1