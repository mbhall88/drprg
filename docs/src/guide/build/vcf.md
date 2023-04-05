# VCF to build reference graph

The input VCF (`-b/--vcf`) for `build` has some very particular requirements. If you
aren't very familiar with the [formal VCF specifications][vcf], please read them
first.

This VCF file is what `drprg` uses to construct the reference graph from. Variants from
it are applied to the [skeleton reference genome](./build.md#reference-genome). It
represents the population variation you would like to build into the reference graph.
For a greater understanding of these concepts, we recommend
reading [this paper][pandora].

To add variants from multiple samples, simply use a multi-sample VCF.

The `CHROM` field should be the gene name. As such, the `POS` is with respect to the
gene. It's important to ensure the `POS` takes into account
the [padding](./build.md#padding) that will be used when running `build`. A helper
script for converting a VCF from reference genome coordinates to gene coordinates can be
found [here][pos-script]. In addition, `drprg` will raise an error if the coordinates
are incorrect with respect to the annotation and reference genome you provide.

[vcf]: https://samtools.github.io/hts-specs/VCFv4.3.pdf

[pandora]: https://github.com/rmcolq/pandora

[pos-script]: https://github.com/mbhall88/drprg/blob/main/scripts/extract_panel_genes_from_vcf.py