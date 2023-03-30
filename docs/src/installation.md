# Installation

## Conda

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