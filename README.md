> ⚠️**In early development stage - do not use**⚠️

👩‍⚕Dr. PRG - **D**rug **R**esistance **P**rediction with **R**eference **G**raphs️👨‍⚕️

[![codecov](https://codecov.io/gh/mbhall88/drprg/branch/main/graph/badge.svg?token=2WAA6MZRKK)](https://codecov.io/gh/mbhall88/drprg)
[![Rust](https://github.com/mbhall88/drprg/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/mbhall88/drprg/actions/workflows/rust.yml)


[TOC]:#

## Table of Contents
- [Install](#install)
  - [Container](#container)
  - [Conda](#conda)
  - [Local](#local)
- [Usage](#usage)
  - [Build](#build)


## Install

### Container

[![Docker Repository on Quay](https://quay.io/repository/mbhall88/drprg/status "Docker Repository on Quay")](https://quay.io/repository/mbhall88/drprg)

A Docker container is available for all commits/branches/versions. To view the available
tags, visit <https://quay.io/repository/mbhall88/drprg?tab=tags>

For example, to use the latest commit on the `main` branch, the URI is

```
$ TAG="latest"
$ URI="quay.io/mbhall88/drprg:$TAG"
```

#### Docker

To run `drprg` using the above container with Docker

```
$ docker pull "$URI"
$ docker run -it "$URI" drprg --help
```

#### Singularity

To run `drprg` using the above container with [Singularity]

```
$ singularity exec "docker://$URI" drprg --help
```

### Conda

<!--todo-->

### Local

```
$ cargo build --release
$ target/release/drprg -h
```

#### Dependencies

`drprg` relies on [`pandora`][pandora] for all subcommands. Additionally, if you want to
build a custom panel, you will need [`make_prg`][makeprg] (prototype) and
[`mafft`][mafft]. To download dependencies and place them in their default location, run

```shell script
# all dependencies
$ make deps
# pandora only
$ make pandora
# make_prg only
$ make makeprg
# mafft only
$ make mafft
```

By default, the external dependencies will be downloaded to `src/ext`. This can be
changed by specifying a path to `EXTDIR` when installing the external dependencies.

```shell script
$ make deps EXTDIR="some/other/dir"
```

## Usage

### Build

Build a `drprg` panel from a [`mykrobe`][mykrobe]-style panel.

```
$ drprg build -a annotation.gff3 -i panel.tsv -f ref.fa -o outdir
```

<!--todo: document the input and output files for build-->

#### Full usage

```
$ drprg build --help
Build a drprg panel from a mykrobe-style panel

USAGE:
    drprg build [FLAGS] [OPTIONS] --gff <gff-file> --panel <panel-file> --fasta <reference-file>

FLAGS:
    -F, --force      Force overwriting existing files. Use this if you want to build from scratch
    -h, --help       Prints help information
        --keep       Keep all temporary files that would otherwise be deleted
    -V, --version    Prints version information
    -v, --verbose    Use verbose output

OPTIONS:
    -a, --gff <gff-file>            Annotation file that will be used to gather information about genes in panel
    -M, --mafft <mafft-exec>        Path to MAFFT executable. Will try in src/ext or $PATH if not given
    -m, --makeprg <makeprg-exec>    Path to make_prg executable. Will try in src/ext or $PATH if not given
    -l, --match-len <match-len>     Minimum number of consecutive characters which must be identical for a match in
                                    make_prg [default: 7]
    -o, --outdir <outdir>           Directory to place output [default: .]
    -P, --padding <padding>         Number of bases of padding to add to start and end of each gene [default: 100]
    -p, --pandora <pandora-exec>    Path to pandora executable. Will try in src/ext or $PATH if not given
    -i, --panel <panel-file>        Panel to build index from
    -f, --fasta <reference-file>    Reference genome in FASTA format (must be indexed with samtools faidx)
    -t, --threads <threads>         Number of threads to use. Use 0 to select the number automatically [default: 1]
```

[pandora]: https://github.com/rmcolq/pandora
[mafft]: https://mafft.cbrc.jp/alignment/software/
[makeprg]: https://github.com/leoisl/make_prg/releases/tag/v0.2.0_prototype
[mykrobe]: https://github.com/Mykrobe-tools/mykrobe
[Singularity]: https://sylabs.io/