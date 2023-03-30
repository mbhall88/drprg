# Expert rules

These are blanket rules that describe resistance (or susceptibility). The file is a CSV
with each row representing a rule and is passed to `drprg build` via the `--rules`
option. The format of each row is

```csv
vartype,gene,start,end,drug
```

1. `vartype`: the variant type of the rule. Supported types are:
    * [`frameshift`][frameshift] - Any insertion or deletion whose length is not a
      multiple of three
    * [`missense`][missense] - A DNA change that results in a different amino acid
    * [`nonsense`][nonsense] - A DNA change that results in a stop codon instead of an
      amino acid
    * `absence` - Gene is absent
2. `gene`: the name of the gene the rule applies to
3. `start`: An optional start position for the rule to apply from. The position is in
   codon coordinates where the rule applies to amino acid changes and is 1-based
   inclusive. If not provided, the start of the gene is inferred. If you want to include
   the upstream (promoter) region of the gene, use negative coordinates.
4. `end`: An optional end position for the rule to apply to. The position is in codon
   coordinates where the rule applies to amino acid changes and is 1-based inclusive. If
   not provided, the end of the gene is inferred.
5. `drug`: A semi-colon-delimited (`;`) list of drugs the rule impacts. If the rule
   confers susceptibility, use `NONE` for this column.

If there are certain rules you need for your
species-of-interest, [raise an issue][issue], and we can look at implementing it.

### Example

This is an example of the *M. tuberculosis* expert rules file used in our paper.

```csv
missense,rpoB,426,452,Rifampicin
nonsense,rpoB,426,452,Rifampicin
frameshift,rpoB,1276,1356,Rifampicin
nonsense,katG,,,Isoniazid
frameshift,katG,,,Isoniazid
absence,katG,,,Isoniazid
nonsense,ethA,,,Ethionamide
frameshift,ethA,,,Ethionamide
absence,ethA,,,Ethionamide
nonsense,gid,,,Streptomycin
frameshift,gid,,,Streptomycin
absence,gid,,,Streptomycin
nonsense,pncA,,,Pyrazinamide
frameshift,pncA,,,Pyrazinamide
absence,pncA,,,Pyrazinamide
missense,katG,315,315,Isoniazid
missense,gid,125,125,Streptomycin
missense,rpoB,425,425,Rifampicin
missense,gid,136,136,Streptomycin
```

The row

```csv
frameshift,pncA,,,Pyrazinamide
```

says that a frameshift *anywhere* within the *pncA* gene will cause resistance to
Pyrazinamide

```csv
nonsense,rpoB,426,452,Rifampicin
frameshift,rpoB,1276,1356,Rifampicin
```

these two rules illustrate the context of the start and end coordinates. In the first
row, we say that any nonsense mutation between 426 and 452 in *rpoB* causes resistance
to Rifampicin. As nonsense mutations only apply to amino acid changes, the coordinates
are in codon-space. Whereas the second row describes a frameshift, which only applies to
nucleotides; therefore, 1276 and 1356 are in bases-space (i.e. the 1276th
nucleotide/base). (As an aside, these two rules both apply to the same region -
the [RRDR])

```csv
missense,katG,315,315,Isoniazid
```

describes any missense mutation at position 315 in *katG* causing isoniazid resistance.


[frameshift]: https://www.genome.gov/genetics-glossary/Frameshift-Mutation

[missense]: https://www.genome.gov/genetics-glossary/Missense-Mutation

[nonsense]: https://www.genome.gov/genetics-glossary/Nonsense-Mutation

[issue]: https://github.com/mbhall88/drprg/issues/new/choose

[RRDR]: https://doi.org/10.1016/j.cmi.2016.09.006
