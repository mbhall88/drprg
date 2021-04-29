> ⚠️**In early development stage - do not use**⚠️

👩‍⚕Dr. PRG - **D**rug **R**esistance **P**rediction with **R**eference **G**raphs️👨‍⚕️

## Install

```
$ cargo build --release
$ target/release/drprg -h
```

### Dependencies

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

[pandora]: https://github.com/rmcolq/pandora
[mafft]: https://mafft.cbrc.jp/alignment/software/
[makeprg]: https://github.com/leoisl/make_prg/releases/tag/v0.2.0_prototype