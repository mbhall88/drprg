rule download_panel:
    output:
        panel=RESOURCES / "panel.original.tsv",
        var2drug=RESOURCES / "var2drug.original.tsv",
    log:
        LOGS / "download_panel.log",
    params:
        url=config["panel"]["url"],
        md5=config["panel"]["md5"],
        panel_fname_pattern="panel.dedup",
        var2drug_fname_pattern="var2drug.dedup",
    container:
        CONTAINERS["base"]
    shadow:
        "shallow"
    resources:
        mem_mb=int(0.5 * GB),
    shell:
        """
        PANEL="panel.tar.gz"
        wget {params.url} -O "$PANEL" 2> {log}
        md5hash=$(md5sum $PANEL | cut -d ' ' -f 1) 2>> {log}
        if [ "$md5hash" != {params.md5} ]; then
            echo "ERROR: md5 hash of downloaded panel ($md5hash) does not match the expected value ({params.md5})"
            exit 1
        fi &>> {log}
        FNAME=$(tar -tzf $PANEL | grep {params.panel_fname_pattern}) 2>> {log}
        tar -zxvOf "$PANEL" "$FNAME" > {output.panel} 2>> {log}
        FNAME=$(tar -tzf $PANEL | grep {params.var2drug_fname_pattern}) 2>> {log}
        tar -zxvOf "$PANEL" "$FNAME" > {output.var2drug} 2>> {log}
        """


rule convert_mutations:
    """There is an inhA promotor mutation which is actually a fabG1 synonymous mutation.
    For drprg it is probably better represented as the synonymous mutation.
    In addition, we need to create a panel that has a column for drug."""
    input:
        panel=rules.download_panel.output.panel,
        var2drug=rules.download_panel.output.var2drug,
    output:
        panel=RESOURCES / "panel.altered.tsv",
    log:
        LOGS / "alter_inhA_mutation.log",
    container:
        CONTAINERS["python"]
    resources:
        mem_mb=int(0.5 * GB),
    params:
        old="inhA_G-154A",  # fabG1_L203L fabG1_ctg606cta fabG1_g608a
        new="fabG1_CTG606CTA",
        new_alphabet="DNA",
    script:
        str(SCRIPTS / "convert_mutations.py")


rule extract_panel_genes_from_popn_vcf:
    input:
        annotation=RESOURCES / "NC_000962.3.gff3",
        vcf=config["population_vcf"],
        index=config["population_vcf"] + ".csi",
        panel=rules.convert_mutations.output.panel,
    output:
        vcf=RESULTS / "drprg/popn_prg/popn.bcf",
    log:
        LOGS / "extract_panel_genes_from_popn_vcf.log",
    params:
        padding=PADDING,
    conda:
        str(ENVS / "extract_panel_genes_from_vcf.yaml")
    script:
        str(SCRIPTS / "extract_panel_genes_from_vcf.py")


rule index_popn_vcf:
    input:
        rules.extract_panel_genes_from_popn_vcf.output.vcf,
    output:
        RESULTS / "drprg/popn_prg/popn.bcf.csi",
    params:
        extra="-f",
    wrapper:
        "0.76.0/bio/bcftools/index"


rule create_references:
    input:
        genome=RESOURCES / "NC_000962.3.fa",
        faidx=RESOURCES / "NC_000962.fa.fai",
        annotation=rules.extract_panel_genes_from_popn_vcf.input.annotation,
        panel=rules.extract_panel_genes_from_popn_vcf.input.panel,
    output:
        fasta=RESULTS / "drprg/popn_prg/genes.fa",
        faidx=RESULTS / "drprg/popn_prg/genes.fa.fai",
    log:
        LOGS / "create_references.log",
    params:
        padding=rules.extract_panel_genes_from_popn_vcf.params.padding,
    conda:
        str(ENVS / "create_references.yaml")
    script:
        str(SCRIPTS / "create_references.py")

rule create_popn_pre_msas:
    input:
        vcf=rules.extract_panel_genes_from_popn_vcf.output,
        vcfidx=rules.index_popn_vcf.output[0],
        references=rules.create_references.output.fasta,
    output:
        directory(RESULTS / "drprg/popn_prg/pre_msas"),
    log:
        LOGS / "create_popn_pre_msas.log",
    conda:
        str(ENVS / "create_popn_pre_msas.yaml")
    script:
        str(SCRIPTS / "create_popn_pre_msas.py")