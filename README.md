👩‍⚕Dr. PRG - **D**rug **r**esistance **P**rediction with **R**eference **G**raphs️👨‍⚕️

[![codecov](https://codecov.io/gh/mbhall88/drprg/branch/main/graph/badge.svg?token=2WAA6MZRKK)](https://codecov.io/gh/mbhall88/drprg)
[![Rust](https://github.com/mbhall88/drprg/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/mbhall88/drprg/actions/workflows/rust.yml)
[![github release version](https://img.shields.io/github/v/release/mbhall88/drprg)](https://github.com/mbhall88/drprg/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI:10.1101/2023.05.04.539481](http://img.shields.io/badge/DOI-10.1101/2023.05.04.539481-B31B1B.svg)](https://doi.org/10.1101/2023.05.04.539481)

Full documentation: <https://mbh.sh/drprg/>

As the name suggests, DrPRG is a tool for predicting drug resistance from sequencing
data. It can be used for any species, provided an index is available for that species.
The [documentation][docs] outlines which species have prebuilt indices and also a guide
for how to create your own.

## Quick usage

Download the latest *M. tuberculosis* prebuilt index

```
drprg index --download mtb
```

Predict resistance from an Illumina fastq

```
drprg predict -x mtb -i reads.fq --illumina -o outdir/
```

### Help

```
$ drprg -h
Drug Resistance Prediction with Reference Graphs

Usage: drprg [OPTIONS] <COMMAND>

Commands:
  build    Build an index to predict resistance from
  predict  Predict drug resistance
  index    Download and interact with indices
  help     Print this message or the help of the given subcommand(s)

Options:
  -v, --verbose        Use verbose output
  -t, --threads <INT>  Maximum number of threads to use [default: 1]
  -h, --help           Print help (see more with '--help')
  -V, --version        Print version
```

## Citation

> Michael B Hall, Leandro Lima, Lachlan J M Coin, Zamin Iqbal. Drug resistance prediction for Mycobacterium tuberculosis with reference graphs. bioRxiv 2023.05.04.539481 [Preprint]. 2023. doi: [10.1101/2023.05.04.539481](https://doi.org/10.1101/2023.05.04.539481)


```bib
@article {Hall2023.05.04.539481,
        author = {Michael B Hall and Leandro Lima and Lachlan J M Coin and Zamin Iqbal},
        title = {Drug resistance prediction for Mycobacterium tuberculosis with reference graphs},
        elocation-id = {2023.05.04.539481},
        year = {2023},
        doi = {10.1101/2023.05.04.539481},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2023/05/04/2023.05.04.539481},
        eprint = {https://www.biorxiv.org/content/early/2023/05/04/2023.05.04.539481.full.pdf},
        journal = {bioRxiv}
}
```

[docs]: https://mbh.sh/drprg/

[guide]: https://mbh.sh/drprg/guide/build/build.html