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
        old="inhA_G-154A",  # fabG1_L203L fabG1_ctg607cta fabG1_g609a
        new="fabG1_CTG607CTA",
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
    log:
        LOGS / "index_popn_vcf.log",
    params:
        extra="-f",
    container:
        CONTAINERS["bcftools"]
    shell:
        "bcftools index {params.extra} {input} 2> {log}"


rule create_references:
    input:
        genome=RESOURCES / "NC_000962.3.fa",
        faidx=rules.index_reference.output.index,
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


rule create_popn_msas:
    input:
        pre_msas=rules.create_popn_pre_msas.output[0],
    output:
        directory(RESULTS / "drprg/popn_prg/msas"),
    log:
        LOGS / "create_popn_msas.log",
    container:
        CONTAINERS["mafft"]
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    shell:
        """
        mkdir -p {output[0]} 2> {log}
        for f in {input.pre_msas}/*.fa
        do
            outname={output[0]}/$(basename "$f")
            mafft --auto --thread {threads} "$f" > "$outname"
        done 2>> {log}
        """


rule make_popn_prgs:
    input:
        msas=rules.create_popn_msas.output[0],
    output:
        prg=RESULTS / "drprg/popn_prg/prgs/dr.prg",
        update_ds=RESULTS / "drprg/popn_prg/prgs/dr.update_DS",
        prgs=directory(RESULTS / "drprg/popn_prg/prgs/dr_prgs"),
    log:
        LOGS / "make_popn_prgs.log",
    shadow:
        "shallow"
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(2 * GB),
    container:
        CONTAINERS["drprg"]
    params:
        match_len=config["match_len"],
        prefix=lambda wildcards, output: Path(output.prg).with_suffix(""),
    shell:
        """
        make_prg from_msa -t {threads} --min_match_len {params.match_len} \
          -o {params.prefix} -i {input.msas} > {log} 2>&1 
        mv {output.prg}.fa {output.prg} 2>> {log}
        """


rule index_popn_prg:
    input:
        prg=rules.make_popn_prgs.output.prg,
        update_ds=rules.make_popn_prgs.output.update_ds,
        prgs=rules.make_popn_prgs.output.prgs,
    output:
        index=RESULTS / "drprg/popn_prg/prgs/w{w}/k{k}/dr.prg.k{k}.w{w}.idx",
        kmer_prgs=directory(
            RESULTS / "drprg/popn_prg/prgs/w{w}/k{k}/kmer_prgs",
        ),
        prg=RESULTS / "drprg/popn_prg/prgs/w{w}/k{k}/dr.prg",
        update_ds=RESULTS / "drprg/popn_prg/prgs/w{w}/k{k}/dr.update_DS",
        prgs=directory(RESULTS / "drprg/popn_prg/prgs/w{w}/k{k}/dr_prgs"),
    log:
        LOGS / "index_popn_prg/w{w}/k{k}.log",
    params:
        options="-v",
        outdir=lambda wildcards, output: Path(output.prg).parent,
    threads: 4
    container:
        CONTAINERS["drprg"]
    shell:
        """
        cp -R {input.prg} {output.prg} 2> {log}
        cp -R {input.update_ds} {output.update_ds} 2>> {log}
        cp -R {input.prgs} {output.prgs} 2>> {log}
        pandora index {params.options} -t {threads} -w {wildcards.w} -k {wildcards.k} \
            {output.prg} &>> {log}
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
        features=rules.create_references.input.annotation,
        reference=rules.create_references.input.genome,
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


rule drprg_build:
    input:
        panel=rules.add_non_resistance_mutations.output.panel,
        ref=rules.create_references.input.genome,
        annotation=rules.create_references.input.annotation,
        prg_index=rules.index_popn_prg.output.index,
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
        padding=rules.create_references.params.padding,
        prebuilt_dir=lambda wildcards, input: Path(input.prg_index).parent,
    shell:
        """
        drprg build {params.options} -l {params.match_len} -P {params.padding} \
            -a {input.annotation} -o {output.outdir} -i {input.panel} \
            -f {input.ref} -t {threads} -d {params.prebuilt_dir} 2> {log}
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
    shell:
        "python {params.script} {params.opts} -i {input.panel} -g {input.gff} -o {output.panel} 2> {log}"


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
