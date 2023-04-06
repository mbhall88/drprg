# Prediction Output

Inside the `predict` output directory (`-o/--outdir`) you'll find a collection of files
and directories. We outline the files that are likely to be of interest to most users.

## Prediction JSON

This is `<sample>.drprg.json`. The value given to the `-s/--sample`
option dictates the name - i.e. `<sample>.drprg.json`.

This is **the only file most users will want/need to interact with**. It contains the
resistance prediction for each drug in the index's catalogue.

#### Example

This is a trimmmed (toy) example JSON output for a sample

```json
{
  "genes": {
    "absent": [
      "ahpC"
    ],
    "present": [
      "embA",
      "embB",
      "ethA",
      "fabG1",
      "gid",
      "gyrA",
      "gyrB",
      "inhA",
      "katG",
      "pncA",
      "rpoB",
      "rrs"
    ]
  },
  "sample": "toy",
  "susceptibility": {
    "Amikacin": {
      "evidence": [
        {
          "gene": "rrs",
          "residue": "DNA",
          "variant": "A1401X",
          "vcfid": "b815ed3f"
        }
      ],
      "predict": "F"
    },
    "Ethambutol": {
      "evidence": [
        {
          "gene": "embB",
          "residue": "PROT",
          "variant": "M306I",
          "vcfid": "a290b118"
        }
      ],
      "predict": "r"
    },
    "Ethionamide": {
      "evidence": [
        {
          "gene": "ethA",
          "residue": "PROT",
          "variant": "A381P",
          "vcfid": "169f75d4"
        }
      ],
      "predict": "U"
    },
    "Isoniazid": {
      "evidence": [
        {
          "gene": "fabG1",
          "residue": "DNA",
          "variant": "G-17T",
          "vcfid": "de9b689e"
        },
        {
          "gene": "katG",
          "residue": "PROT",
          "variant": "S315T",
          "vcfid": "acaa8ca2"
        }
      ],
      "predict": "R"
    },
    "Levofloxacin": {
      "evidence": [],
      "predict": "S"
    }
  },
  "version": {
    "drprg": "0.1.0",
    "index": "20230308"
  }
}
```

The keys of the JSON are

- `genes`: This contains a list of genes in the index reference graph which are present
  and absent
- `sample`: The value passed to the `-s/--sample` option
- `susceptibility`: The keys of this entry are the drugs in the index catalogue. Each
  drug's entry contains [`evidence`](#evidence) supporting the value in
  the [`predict`](#predict) section.

### Predict

The `predict` entry for a drug is the resistance prediction for the sample. Possible
values are

- `S`: susceptible. This is the "default" prediction. If no mutations are detected for
  the sample, it is assumed to be susceptible
- `F`: failed. Genotyping failed for one or more mutations for this drug. See
  the [prediction VCF](#prediction-vcf) for more information
- `U`: unknown. One or more mutations that are not present in the index catalogue were
  detected in a gene associated with this drug
- `R`: resistant. One or more mutations from the index catalogue that confer resistance
  were detected
- `u` or `r`: The same as the uppercase versions, but the mutation(s) were detected in
  a [minor allele](./predict.md#minimum-allele-frequency).

### Evidence

This is a list of the mutations supporting the prediction. The `residue` is one of `DNA`
or `PROT` indicating whether the mutation describes a nucleotide or amino acid change,
respectively.

The `variant` is of the form `<ref><pos><alt>`; where `<ref>` is the reference sequence
at position `<pos>` and `<alt>` is the nucleotide/amino acid the reference is changed
to. See the [catalogue docs](./build/catalogue.md) for more information.

The `vcfid` is the value in [the VCF](#prediction-vcf) `ID` column for this mutation,
making it easier to find a mutation in [the VCF](#prediction-vcf).

## Prediction VCF

This is `<sample>.drprg.bcf`. As this file is a [BCF], you will need to
use [`bcftools`][bcftools] to view it - e.g. `bcftools view sample.drprg.bcf`.

You should only need to interact with this file if you want further information about
the exact evidence supporting a mutation being called, why a mutation was called as
failed (`F`), or why a mutation *wasn't* called. For those mutations in the JSON, you
can easily look them up using the `vcfid`, which can be found in the `ID` (third)
column of the BCF. Or you can just use [`grep`][grep]. For example, to look up the
Isoniazid S315T mutation
in *katG* from [the example](#example) you can use

```
$ bcftools view sample.drprg.bcf | grep acaa8ca2
katG    1044    acaa8ca2        GC      AC,CA,CC        .       PASS    VC=PH_SNPs;GRAPHTYPE=SIMPLE;PDP=0,0.0123457,0,0.987654;VARID=katG_S315T;PREDICT=R       GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF      3:0,1,0,42:0,0,0,38:0,1,0,42:0,0,0,38:0,1,0,127:0,0,0,116:1,1,1,0:-523.019,-514.096,-523.019,-7.87925:506.217
```

All `INFO` and `FORMAT` fields are defined in the header of the BCF file. We recommend
reading the [VCF/BCF][BCF] specifications for help with how to interpret data in a VCF
file.

[BCF]: https://samtools.github.io/hts-specs/VCFv4.4.pdf

[bcftools]: https://github.com/samtools/bcftools/

[grep]: https://www.man7.org/linux/man-pages/man1/grep.1.html