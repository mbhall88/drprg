# Installation

## Conda

[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/drprg)](https://anaconda.org/bioconda/drprg)
[![bioconda version](https://anaconda.org/bioconda/drprg/badges/platforms.svg)](https://anaconda.org/bioconda/drprg)
![Conda](https://img.shields.io/conda/dn/bioconda/drprg)

```
conda install drprg
```

## Container

[![Docker Repository on Quay](https://quay.io/repository/mbhall88/drprg/status "Docker Repository on Quay")](https://quay.io/repository/mbhall88/drprg)

A Docker container is available for all commits/branches/versions. To view the available
tags, visit <https://quay.io/repository/mbhall88/drprg?tab=tags>

For example, to use the latest commit on the `main` branch, the URI is

```
$ TAG="latest"
$ URI="quay.io/mbhall88/drprg:$TAG"
```

### Docker

To run `drprg` using the above container with Docker

```
$ docker pull "$URI"
$ docker run -it "$URI" drprg --help
```

### Singularity

To run `drprg` using the above container with [Singularity]

```
$ singularity exec "docker://$URI" drprg --help
```

## Prebuilt binary

If you use the prebuilt binary, you must have the [external dependecies](#dependencies) installed separately.

```shell
curl -sSL drprg.mbh.sh | sh
# or with wget
wget -nv -O - drprg.mbh.sh | sh
```

You can also pass options to the script like so

```
$ curl -sSL drprg.mbh.sh | sh -s -- --help
install.sh [option]

Fetch and install the latest version of drprg, if drprg is already
installed it will be updated to the latest version.

Options
        -V, --verbose
                Enable verbose output for the installer

        -f, -y, --force, --yes
                Skip the confirmation prompt during installation

        -p, --platform
                Override the platform identified by the installer

        -b, --bin-dir
                Override the bin installation directory [default: /usr/local/bin]

        -a, --arch
                Override the architecture identified by the installer [default: x86_64]

        -B, --base-url
                Override the base URL used for downloading releases [default: https://github.com/mbhall88/drprg/releases]

        -h, --help
                Display this help message
```


## Cargo

[![Crates.io](https://img.shields.io/crates/v/drprg.svg)](https://crates.io/crates/drprg)

If installing via cargo, you must have the [external dependecies](#dependencies) installed separately.

```
$ cargo install drprg
```

## Local

**Minimum supported Rust version**: `1.65.0`

```
$ cargo build --release
$ target/release/drprg -h
```

### Dependencies

`drprg` relies on:
- [`pandora`][pandora]
- [`bcftools`][bcftools]
- [`make_prg`][makeprg]
- [`mafft`][mafft]

You can install the dependencies using the provided [`justfile`][just]

```shell script
# all dependencies
$ just deps

# pandora only
$ just pandora

# make_prg only
$ just makeprg

# mafft only
$ just mafft

# bcftools only
$ just bcftools
```

By default, the external dependencies will be downloaded to `src/ext`. This can be
changed by specifying a path to `EXTDIR` when installing the external dependencies.

```shell script
$ just deps EXTDIR="some/other/dir"
```

[pandora]: https://github.com/rmcolq/pandora
[mafft]: https://mafft.cbrc.jp/alignment/software/
[makeprg]: https://github.com/leoisl/make_prg/
[mykrobe]: https://github.com/Mykrobe-tools/mykrobe
[Singularity]: https://sylabs.io/
[bcftools]: https://samtools.github.io/bcftools/bcftools.html
[just]: https://github.com/casey/just
