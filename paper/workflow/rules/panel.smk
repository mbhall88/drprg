rule download_panel:
    output:
        panel=RESOURCES / "panel.original.tsv",
        var2drug=RESOURCES / "var2drug.original.json",
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
        tar -zxOf "$PANEL" "$FNAME" > {output.panel} 2>> {log}
        FNAME=$(tar -tzf $PANEL | grep {params.var2drug_fname_pattern}) 2>> {log}
        tar -zxOf "$PANEL" "$FNAME" > {output.var2drug} 2>> {log}
        """


rule convert_mutations:
    """There is an inhA promotor mutation which is actually a fabG1 synonymous mutation.
    For drprg it is probably better represented as the synonymous mutation.
    In addition, we need to create a panel that has a column for drug."""
    input:
        panel=rules.download_panel.output.panel,
        var2drug=rules.download_panel.output.var2drug,
    output:
        panel=RESOURCES / "panel.converted.tsv",
    log:
        LOGS / "convert_mutations.log",
    container:
        CONTAINERS["python"]
    resources:
        mem_mb=int(0.5 * GB),
    params:
        remove_mutations={"inhA_G-154A"},  # is actually fabG1_L203L
        update_drugs={
            "fabG1_CTG607CTC": ["Isoniazid", "Ethionamide"],
            "fabG1_CTG607CTT": ["Isoniazid", "Ethionamide"],
            "fabG1_CTG607TTA": ["Isoniazid", "Ethionamide"],
            "fabG1_CTG607TTG": ["Isoniazid", "Ethionamide"],
            "fabG1_CTG607CTA": ["Isoniazid", "Ethionamide"],
        },
    script:
        str(SCRIPTS / "convert_mutations.py")


rule create_orphan_mutations:
    input:
        mutations=RESOURCES / "target.mutations.tsv",
        reference=RESOURCES / "NC_000962.3.fa",
        annotation=RESOURCES / "NC_000962.3.gff3",
        gene_reference=WORKFLOW.parent.parent / "tests/cases/predict/genes.fa",
    output:
        vcf=RESULTS / "drprg/popn_prg/common_mutations.bcf",
        vcfidx=RESULTS / "drprg/popn_prg/common_mutations.bcf.csi",
    log:
        LOGS / "create_orphan_mutations.log",
    resources:
        mem_mb=int(0.5 * GB),
    shadow:
        "shallow"
    params:
        offset=config["padding"],
    conda:
        ENVS / "create_orphan_mutations.yaml"
    script:
        SCRIPTS / "create_orphan_mutations.py"


rule extract_panel_genes_from_popn_vcf:
    input:
        annotation=RESOURCES / "NC_000962.3.gff3",
        vcf=config["population_vcf"],
        panel=rules.convert_mutations.output.panel,
    output:
        vcf=RESULTS / "drprg/popn_prg/full.merged.bcf",
    log:
        LOGS / "extract_panel_genes_from_popn_vcf.log",
    params:
        padding=PADDING,
    conda:
        str(ENVS / "extract_panel_genes_from_vcf.yaml")
    script:
        str(SCRIPTS / "extract_panel_genes_from_vcf.py")


rule merge_reference_vcfs:
    input:
        popn_vcf=rules.extract_panel_genes_from_popn_vcf.output.vcf,
        mutations_vcf=rules.create_orphan_mutations.output.vcf,
        reference=rules.create_orphan_mutations.input.gene_reference,
    output:
        vcf=RESULTS / "drprg/popn_prg/final.bcf",
        vcfidx=RESULTS / "drprg/popn_prg/final.bcf.csi",
    log:
        LOGS / "merge_reference_vcfs.log",
    shadow:
        "shallow"
    resources:
        mem_mb=int(0.5 * GB),
    container:
        CONTAINERS["bcftools"]
    shell:
        """
        (bcftools merge {input.mutations_vcf} {input.popn_vcf} \
            | bcftools norm -f {input.reference} -c e -o {output.vcf} -) 2> {log}
        bcftools index -f {output.vcf} 2>> {log}
        """


rule download_who_panel:
    output:
        panel=RESOURCES / "who.panel.tsv",
    log:
        LOGS / "download_who_panel.log",
    container:
        CONTAINERS["base"]
    params:
        url=config["who_panel_url"],
    shell:
        "wget {params.url} -O {output.panel} 2> {log}"


rule add_non_resistance_mutations:
    input:
        panel=rules.convert_mutations.output.panel,
        known=rules.download_who_panel.output.panel,
        features=rules.extract_panel_genes_from_popn_vcf.input.annotation,
        reference=RESOURCES / "NC_000962.3.fa",
    output:
        panel=RESOURCES / "panel.with_susceptible_mutations.tsv",
    log:
        LOGS / "add_non_resistance_mutations.log",
    container:
        CONTAINERS["python"]
    resources:
        mem_mb=int(0.5 * GB),
    params:
        keep_grades=(4, 5),
        no_drug="NONE",
        exclude=["inhA_T4I"],  # incorrect variant that should be fabG1_T4I which is already in catalogue
    script:
        str(SCRIPTS / "add_non_resistance_mutations.py")


rule filter_panel_for_expert_rules:
    input:
        panel=rules.add_non_resistance_mutations.output.panel,
        rules=RESOURCES / "rules.csv",
    output:
        panel=RESOURCES / "panel.filtered.tsv",
        rules=RESOURCES / "rules.extra.csv",
    log:
        LOGS / "filter_panel_for_expert_rules.log",
    container:
        CONTAINERS["python"]
    params:
        script=SCRIPTS / "filter_panel_for_expert_rules.py",
    shell:
        """
        python {params.script} {input.panel} {input.rules} {output.panel} {output.rules} 2> {log}
        """


rule drprg_build:
    input:
        panel=rules.filter_panel_for_expert_rules.output.panel,
        ref=rules.add_non_resistance_mutations.input.reference,
        annotation=rules.extract_panel_genes_from_popn_vcf.input.annotation,
        vcf=rules.merge_reference_vcfs.output.vcf,
        vcfidx=rules.merge_reference_vcfs.output.vcfidx,
        rules=rules.filter_panel_for_expert_rules.output.rules,
    output:
        outdir=directory(RESULTS / "drprg/index/w{w}/k{k}"),
        prg=RESULTS / "drprg/index/w{w}/k{k}/dr.prg",
        vcf=RESULTS / "drprg/index/w{w}/k{k}/panel.bcf",
        vcf_idx=RESULTS / "drprg/index/w{w}/k{k}/panel.bcf.csi",
        ref=RESULTS / "drprg/index/w{w}/k{k}/genes.fa",
    log:
        LOGS / "drprg_build/w{w}/k{k}.log",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    threads: 2
    container:
        CONTAINERS["drprg"]
    params:
        options="-v -w {w} -k {k}",
        match_len=config["match_len"],
        padding=config["padding"],
    shell:
        """
        drprg build {params.options} -l {params.match_len} -P {params.padding} \
            -a {input.annotation} -o {output.outdir} -i {input.panel} \
            -f {input.ref} -t {threads} -b {input.vcf} -r {input.rules} 2> {log}
        """


rule download_tbprofiler_db:
    output:
        db=directory(RESULTS / "tbprofiler/tbdb"),
    container:
        CONTAINERS["base"]
    log:
        LOGS / "download_tbprofiler_db.log",
    params:
        url=config["tbdb_url"],
        outdir=lambda wildcards, output: Path(output.db).parent,
    shadow:
        "shallow"
    shell:
        """
        wget {params.url} -O tbdb.zip 2> {log}
        unzip -d {params.outdir} tbdb.zip 2>> {log}
        mv {params.outdir}/tbdb* {output.db} 2>> {log}
        # clear the other annotations file
        > {output.db}/tbdb.other_annotations.csv
        # clear the watchlist as it can cause error
        > {output.db}/tbdb.watchlist.csv
        """


rule mykrobe_to_hgvs:
    input:
        panel=rules.convert_mutations.output.panel,
        gff=RESOURCES / "NC_000962.3.gff3",
    output:
        panel=RESOURCES / "mykrobe_to_hgvs.csv",
    log:
        LOGS / "mykrobe_to_hgvs.log",
    container:
        CONTAINERS["python"]
    params:
        script=SCRIPTS / "mykrobe_to_hgvs.py",
        opts="-v",
        expert_rules=config["expert_rules"],
    shell:
        """
        python {params.script} {params.opts} -E {params.expert_rules} \
            -i {input.panel} -g {input.gff} -o {output.panel} 2> {log}
        """


rule create_tbprofiler_db:
    input:
        panel=rules.mykrobe_to_hgvs.output.panel,
        db=rules.download_tbprofiler_db.output.db,
    output:
        touch(RESULTS / "tbprofiler/.db.built"),
    log:
        LOGS / "create_tbprofiler_db.log",
    conda:
        str(ENVS / "tbprofiler.yaml")
    params:
        opts="--load --custom --include_original_mutation",
    shell:
        """
        cp {input.panel} {input.db}/tbdb.csv 2> {log}
        cd {input.db} || exit 1
        tb-profiler create_db {params.opts} 2>> {log}
        """
