# Build an index

The `build` subcommand is used to build a custom index to [`predict`](../predict.md)
from.

```
$ drprg build -a annotation.gff3 -i catalogue.tsv -f ref.fa -o outdir
```

## Required

### Catalogue

Also referred to as a "panel" sometimes. This is the catalogue of mutations that confer
resistance (or susceptibility). See [the Catalogue docs](./catalogue.md) for more
details of what format this file must take. Provided via the `-i/--panel` option.

### Reference genome

A FASTA file of the reference genome that will act as the "skeleton" for the reference
graph. Provided via the `-f/--fasta` option.

### Annotation

The annotation is a [GFF3][gff] file for the species you are building the index for.
Provided via the `-a/--gff` option.

### VCF

A VCF to build the reference graph from. See [VCF](./vcf.md) for a detailed description
of this file. Provided via the `-b/--vcf` option. This option is mutually exclusive
with [`--prebuilt-prg`](#prebuilt-reference-graph).

## Optional

### Expert rules

These describe a set of rules - rather than a single mutation - that cause resistance.
See [Expert rules](./expert.md) for a detailed description of the rules. Provided via
the `-r/--rules` option.

### Padding

Number of bases to add to the start and end of each gene (locus) in the reference graph.
Set to 100 by default. Provided via the `-P/--padding` option.

### Prebuilt reference graph

**Advanced users only.** For those who know how to build their own reference graphs, and
would prefer to do that themselves by hand, rather than let DrPRG do it for them. The
value passed to this option must be a directory containing:

- A reference graph file named `dr.prg`
- A directory of multiple sequence alignments (`/msas/`) that the reference graph was
  built from. There must be a fasta file alignment for each gene in this
  directory `<gene>.fa`.
- An optional [`pandora`][pandora] index. If not present, one will be created.

The reference graph is expected to contain the reference sequence (
with [padding](#padding)) for each gene according to the [annotation](#annotation)
and [reference sequence](#reference-genome) provided.

Provided with the `-d/--prebuilt-prg` option.

### Version

The version name for the index. Set to the current date (`YYYYMMDD`) by default.
Provided via the `--version` option.

## Quick usage

```
$ drprg build -h
Build an index to predict resistance from

Usage: drprg build [OPTIONS] --gff <FILE> --panel <FILE> --fasta <FILE>

Options:
  -v, --verbose             Use verbose output
  -t, --threads <INT>       Maximum number of threads to use [default: 1]
  -a, --gff <FILE>          Annotation file that will be used to gather information about genes in catalogue
  -i, --panel <FILE>        Panel/catalogue to build index for
  -f, --fasta <FILE>        Reference genome in FASTA format (must be indexed with samtools faidx)
  -b, --vcf <FILE>          An indexed VCF to build the index PRG from. If not provided, then a prebuilt PRG must be given. See `--prebuilt-prg`
  -P, --padding <INT>       Number of bases of padding to add to start and end of each gene [default: 100]
  -o, --outdir <DIR>        Directory to place output [default: .]
  -d, --prebuilt-prg <DIR>  A prebuilt PRG to use
  -r, --rules <FILE>        "Expert rules" to be applied in addition to the catalogue
  -I, --no-fai              Don't index --fasta if an index doesn't exist
  -C, --no-csi              Don't index --vcf if an index doesn't exist
      --version <VERSION>   Version to use for the index [default: 20230404]
  -h, --help                Print help (see more with '--help')
```

## Full usage

```
$ drprg build --help
Build an index to predict resistance from

Usage: drprg build [OPTIONS] --gff <FILE> --panel <FILE> --fasta <FILE>

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

  -B, --bcftools <FILE>
          Path to bcftools executable. Will try in src/ext or $PATH if not given

  -a, --gff <FILE>
          Annotation file that will be used to gather information about genes in catalogue

  -i, --panel <FILE>
          Panel/catalogue to build index for

  -f, --fasta <FILE>
          Reference genome in FASTA format (must be indexed with samtools faidx)

  -b, --vcf <FILE>
          An indexed VCF to build the index PRG from. If not provided, then a prebuilt PRG must be given. See `--prebuilt-prg`

  -P, --padding <INT>
          Number of bases of padding to add to start and end of each gene

          [default: 100]

  -o, --outdir <DIR>
          Directory to place output

          [default: .]

  -d, --prebuilt-prg <DIR>
          A prebuilt PRG to use.

          Only build the panel VCF and reference sequences - not the PRG. This directory MUST contain a PRG file named `dr.prg`, along, with a directory called `msas/` that contains an MSA fasta file for each gene `<gene>.fa`. There can optionally also be a pandora index file, but if not, the indexing will be performed by drprg. Note: the PRG is expected to contain the reference sequence for each gene according to the annotation and reference genome given (along with padding) and must be in the forward strand orientation.

  -r, --rules <FILE>
          "Expert rules" to be applied in addition to the catalogue.

          CSV file with blanket rules that describe resistance (or susceptibility). The columns are <variant type>,<gene>,<start>,<end>,<drug(s)>. See the docs for a detailed explanation.

  -l, --match-len <INT>
          Minimum number of consecutive characters which must be identical for a match in make_prg

          [default: 5]

  -N, --max-nesting <INT>
          Maximum nesting level when constructing the reference graph with make_prg

          [default: 5]

  -k, --pandora-k <INT>
          Kmer size to use for pandora

          [default: 15]

  -w, --pandora-w <INT>
          Window size to use for pandora

          [default: 11]

  -I, --no-fai
          Don't index --fasta if an index doesn't exist

  -C, --no-csi
          Don't index --vcf if an index doesn't exist

      --version <VERSION>
          Version to use for the index

          [default: 20230404]

  -h, --help
          Print help (see a summary with '-h')
```

[gff]: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

[pandora]: https://github.com/rmcolq/pandora