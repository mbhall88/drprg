👩‍⚕Dr. PRG - **D**rug **r**esistance **P**rediction with **R**eference **G**raphs️👨‍⚕️

[![codecov](https://codecov.io/gh/mbhall88/drprg/branch/main/graph/badge.svg?token=2WAA6MZRKK)](https://codecov.io/gh/mbhall88/drprg)
[![Rust](https://github.com/mbhall88/drprg/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/mbhall88/drprg/actions/workflows/rust.yml)

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

## Citation

Preprint coming soon...

[docs]: https://mbh.sh/drprg/

[guide]: https://mbh.sh/drprg/guide/build/build.html