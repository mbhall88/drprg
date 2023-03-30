# Build an index


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

