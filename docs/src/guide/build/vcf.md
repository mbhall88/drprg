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


## Example from *Mycobacterium* tuberculosis

We provide an example of how the VCF for the [work on *M. tuberculosis*][paper] (MTB) was generated.

### Call variants for population

The first step is to select the population of samples you want variation from. For the MTB reference graph, we took a subset of the [15,211 global MTB isolates published by the CRyPTIC consortium][cryptic]. Briefly, the 15,211 samples were individually variant-called using [`clockwork`][clockwork] and then joint-genotyped using [`minos`][minos]. The "wide" VCF obtained at the end contains genotypes for all 15,211 samples against all variants discovered across the population. 

**Note**: the reference genome you call these variants against must be the same the [skeleton reference](./build.md#reference-genome).

### (Optional) subset population

In our running example, we have 15,211 samples in our VCF. This is far too much to be useful for building a reference graph. [This paper][forge] illustrates the problem nicely. In addition, Figure 3.2 and Sections 3.6.4 and 3.6.5 in [this thesis][thesis] show that using too many samples does not provide any benefit to variant precision or recall, but costs increased memory usage and CPU time.

To create a representative subset, we split the VCF file into separate lineage VCFs, where all samples in a VCF are of the same lineage. We randomly chose 20 samples from each lineage 1 through 4, as well as 20 samples from all other lineages combined. In addition, we included [17 clinical samples][gramtools] representing MTB global diversity (lineages 1-6) to give a total of 117 samples.

We generated the subsets by listing all samples in the VCF file

```
bcftools query --list-samples > samples.txt
```

and then splitting this list by lineage based on lineage classifications. If we have a list of Lineage 1 samples in a text file, we can get 20 random samples from it with

```
shuf -n 20 L1.txt > L1.subsampled.txt
```

we then combine all of these text files of sample names into a single text file, extract the 117 samples, and filter the remaining variants with

```
bcftools view -S samples.subsampled.txt wide.vcf.gz |
  bcftools view -a -U -c=1:nref -o MTB.vcf.gz
```

`-a` trims ALT alleles not seen in the genotype fields, `-U` excludes sites without a called genotype, and `-c=1:nref` removes positions there is no non-reference alleles.

You can find the resulting VCF (`MTB.vcf.gz`) we generated [here][sparse].

### (Optional) Manually add "orphan" mutations

We also manually added ~50 mutations that were commonly found in minor frequencies but were not contained in the population VCF.

This is done by creating a file of mutations of the form `<gene>_<mutation>` - this is effectively column 1 and 2 from [the catalogue](./catalogue.md) separated by an underscore. For example

```
rpoB_I491F
rplC_C154R
ddn_L49P
gid_Q125*
rrs_g1484t
embB_L74R
```

the difference between this and the catalogue though is that DNA mutations should have the nucleotides in lowercase (e.g. `rrs_g1484t`). A VCF can be generated for this list of mutations using [this script][orphan_script].

```
python create_orphan_mutations.py -r ref.fa -m orphan_mutations.txt \
  -g annotation.gff -o orphan.vcf
```

We then merge this VCF of "orphan" mutations into the population VCF with

```
bcftools merge MTB.vcf.gz orphan.vcf |
  bcftools norm -c e -f reference.fa -o merged.vcf
```

`-c e` just tells `bcftools norm` to error if the REF alleles in the VCFs do not match the reference.

### Extract catalogue genes from VCF

DrPRG expects the VCF `CHROM` field to be the name of the gene, rather than the reference contig/chromosome name - which is what we currently have. As such, the `POS` must be with respect to the gene. It's important to ensure the `POS` takes into account the [padding](./build.md#padding) that will be used when running `build`. 

A helper script for converting a VCF from reference genome coordinates to gene coordinates can be found [here][pos-script]. In addition, `drprg` will raise an error if the coordinates are incorrect with respect to the annotation and reference genome you provide. By providing this script with the [catalogue/panel](./catalogue.md) you will use in [`build`](./build.md), it will extract only those variants that fall within the genes in your catalogue.

```
python extract_panel_genes_from_vcf.py --padding 100 \
  -i catalogue.tsv -g annotation.gff --vcf merged.vcf -o final.vcf
```

`final.vcf` can now be used as the [input VCF for `build`](./build.md#vcf).


[vcf]: https://samtools.github.io/hts-specs/VCFv4.3.pdf

[pandora]: https://github.com/rmcolq/pandora

[paper]: https://doi.org/10.1101/2023.05.04.539481

[pos-script]: https://github.com/mbhall88/drprg/blob/main/scripts/extract_panel_genes_from_vcf.py

[cryptic]: https://doi.org/10.1371/journal.pbio.3001721

[clockwork]: https://github.com/iqbal-lab-org/clockwork

[minos]: https://github.com/iqbal-lab-org/minos

[forge]: https://doi.org/10.1186/s13059-018-1595-x

[thesis]: https://doi.org/10.17863/CAM.81350

[gramtools]: https://doi.org/10.1186/s13059-021-02474-0

[sparse]: https://doi.org/10.6084/m9.figshare.23625327

[orphan_script]: https://github.com/mbhall88/drprg/blob/main/scripts/create_orphan_mutations.py