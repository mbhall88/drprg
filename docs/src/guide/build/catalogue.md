# Catalogue of mutations

Also referred to as a "panel" sometimes. This is the catalogue of mutations that confer
resistance (or susceptibility). Provided via the `-i/--panel` option.

This file is a tab-delimited (TSV) file with four columns

1. The gene name of the mutation
2. The mutation in the form `<ref><pos><alt>`, where `<ref>` is the reference nucleotide
   or amino acid, `<pos>` is the 1-based position in the gene/protein, and `<alt>` is
   the nucleotide or amino acid the reference is changed to. `<pos>` is with respect to
   the type of mutation - i.e. if the mutation is a amino acid change, the `<pos>` must
   be position within the protein (codon position). Promoter mutations can be specied
   with a negative symbol - e.g. -10 means 10 nucleotide before the first position in
   the gene. One-letter codes are to be used for amino acids.
3. The residue type - `DNA` (nucleotide change) or `PROT` (amino acid change).
4. Comma-separated list of the drug(s) the mutation confers resistance to. Use `NONE` if
   the mutation is not associated with resistance. 

## Example

```
pncA    TCG196TAG   DNA Pyrazinamide
pncA    T142R   PROT    Pyrazinamide
tlyA    S159*   PROT    Capreomycin
embB    M306N   PROT    Ethambutol
rpoB    C1275CCA    DNA Rifampicin
gid R137P   PROT    Streptomycin
gid R118S   PROT    Streptomycin
tlyA    G123GC  DNA Capreomycin
pncA    G97D    PROT    Pyrazinamide
fabG1   C-15X   DNA Ethionamide,Isoniazid
```