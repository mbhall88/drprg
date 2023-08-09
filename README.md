👩‍⚕Dr. PRG - **D**rug **r**esistance **P**rediction with **R**eference **G**raphs️👨‍⚕️

[![codecov](https://codecov.io/gh/mbhall88/drprg/branch/main/graph/badge.svg?token=2WAA6MZRKK)](https://codecov.io/gh/mbhall88/drprg)
[![Rust](https://github.com/mbhall88/drprg/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/mbhall88/drprg/actions/workflows/rust.yml)
[![github release version](https://img.shields.io/github/v/release/mbhall88/drprg)](https://github.com/mbhall88/drprg/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![10.1099/mgen.0.001081](https://img.shields.io/badge/doi-10.1099%2Fmgen.0.001081-blue)](https://doi.org/10.1099/mgen.0.001081)

Full documentation: <https://mbh.sh/drprg/>

As the name suggests, Dr. PRG (pronounced "Doctor P-R-G") is a tool for predicting drug resistance from sequencing
data. It can be used for any species, provided an index is available for that species.
The [documentation][docs] outlines which species have prebuilt indices and also a guide
for how to create your own.

## Quick Installation

```
conda install -c bioconda drprg
```

Linux is currently the only supported platform; however, there is a Docker container that can be used on other platforms. 

See the [installation guide][install] for more options.

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

> Hall MB, Lima L, Coin LJM, Iqbal Z (2023) Drug resistance prediction for Mycobacterium tuberculosis with reference graphs. Microbial Genomics 9:001081. doi: [10.1099/mgen.0.001081](https://doi.org/10.1099/mgen.0.001081) 


```bib
@article{hall_drug_2023,
	title = {Drug resistance prediction for {Mycobacterium} tuberculosis with reference graphs},
	volume = {9},
	copyright = {All rights reserved},
	issn = {2057-5858},
	url = {https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001081},
	doi = {10.1099/mgen.0.001081},
	number = {8},
	journal = {Microbial Genomics},
	author = {Hall, Michael B. and Lima, Leandro and Coin, Lachlan J. M. and Iqbal, Zamin},
	year = {2023},
	pages = {001081},
}
```

[docs]: https://mbh.sh/drprg/

[guide]: https://mbh.sh/drprg/guide/build/build.html
[install]: https://mbh.sh/drprg/installation.html
