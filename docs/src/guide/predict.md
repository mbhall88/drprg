# Predict

The `predict` subcommand is used to predict resistance for a sample from an index.

At its simplest

```
drprg predict -i reads.fq -x mtb -o outdir
```

`drprg` is a bit "new-age" in that it assume the reads are Nanopore. If they're
Illumina, use the `-I/--illumina` option.

See [Predict Output](./predict-output.md) documentation for a detailed description of
what results/output files and formats to expect.

## Required

### Index

The index is provided via the `-x/--index` option. It can either be a path to an index,
or the name of a [downloaded index](./download.md). As with the [`index`](./download.md)
subcommmand, you can specify a version if you don't want to use the latest.

### Input reads

A fastq (or fasta) file of the reads you want to predict resistance from - provided via
the `-i/--input` option. If you have paired reads in two files, simply combine them and
pass the combined file - interleave order doesn't matter. For example

```
cat r1.fq r2.fq > combined.fq
drprg predict -i combined.fq ...
```

`gzip`-compressed files are also accepted.

## Optional

### Sample name

Identifier to use for your output files. By default, it will be set to the file name
prefix (e.g. `name` for a fastq named `name.fq.gz`). Provided via the `-s/--sample`
option.

### Minimum allele frequency

Provided via the `-f/--maf` option. If an alternate allele has at least this fraction of
the depth, a minor resistance ("r") prediction is made. By default, this is set to `1.0`
for Nanopore data (i.e. minor allele detection is off) and `0.1` when using
the `--illumina` option. For example, if a variant is called as the reference allele for
Illumina reads, but an alternate allele has more than 10% of the depth on that position,
a minor resistance call is made for the alternate allele.

### Ignore synonymous

Using the `-S/--ignore-synonymous` option will prevent synonymous mutations from
appearing as unknown resistance calls. However, any synonymous mutations in the
catalogue will still be considered.

## Quick usage

```
$ drprg predict -h
Predict drug resistance

Usage: drprg predict [OPTIONS] --index <DIR> --input <FILE>

Options:
  -v, --verbose        Use verbose output
  -t, --threads <INT>  Maximum number of threads to use [default: 1]
  -h, --help           Print help (see more with '--help')

Input/Output:
  -x, --index <DIR>      Name of a downloaded index or path to an index
  -i, --input <FILE>     Reads to predict resistance from
  -o, --outdir <DIR>     Directory to place output [default: .]
  -s, --sample <SAMPLE>  Identifier to use for the sample
  -I, --illumina         Sample reads are from Illumina sequencing

Filter:
  -S, --ignore-synonymous     Ignore unknown (off-catalogue) variants that cause a synonymous substitution
  -f, --maf <FLOAT[0.0-1.0]>  Minimum allele frequency to call variants [default: 1]
```

## Full usage

```
$ drprg predict --help
Predict drug resistance

Usage: drprg predict [OPTIONS] --index <DIR> --input <FILE>

Options:
  -p, --pandora <FILE>
          Path to pandora executable. Will try in src/ext or $PATH if not given

  -v, --verbose
          Use verbose output

  -m, --makeprg <FILE>
          Path to make_prg executable. Will try in src/ext or $PATH if not given

  -t, --threads <INT>
          Maximum number of threads to use

          Use 0 to select the number automatically

          [default: 1]

  -M, --mafft <FILE>
          Path to MAFFT executable. Will try in src/ext or $PATH if not given

  -h, --help
          Print help (see a summary with '-h')

Input/Output:
  -x, --index <DIR>
          Name of a downloaded index or path to an index

  -i, --input <FILE>
          Reads to predict resistance from

          Both fasta and fastq are accepted, along with compressed or uncompressed.

  -o, --outdir <DIR>
          Directory to place output

          [default: .]

  -s, --sample <SAMPLE>
          Identifier to use for the sample

          If not provided, this will be set to the input reads file path prefix

  -I, --illumina
          Sample reads are from Illumina sequencing

Filter:
  -S, --ignore-synonymous
          Ignore unknown (off-catalogue) variants that cause a synonymous substitution

  -f, --maf <FLOAT[0.0-1.0]>
          Minimum allele frequency to call variants

          If an alternate allele has at least this fraction of the depth, a minor resistance ("r") prediction is made. Set to 1 to disable. If --illumina is passed, the default is 0.1

          [default: 1]

      --debug
          Output debugging files. Mostly for development purposes

  -d, --min-covg <INT>
          Minimum depth of coverage allowed on variants

          [default: 3]

  -D, --max-covg <INT>
          Maximum depth of coverage allowed on variants

          [default: 2147483647]

  -b, --min-strand-bias <FLOAT>
          Minimum strand bias ratio allowed on variants

          For example, setting to 0.25 requires >=25% of total (allele) coverage on both strands for an allele.

          [default: 0.01]

  -g, --min-gt-conf <FLOAT>
          Minimum genotype confidence (GT_CONF) score allow on variants

          [default: 0]

  -L, --max-indel <INT>
          Maximum (absolute) length of insertions/deletions allowed

  -K, --min-frs <FLOAT>
          Minimum fraction of read support

          For example, setting to 0.9 requires >=90% of coverage for the variant to be on the called allele

          [default: 0]
```