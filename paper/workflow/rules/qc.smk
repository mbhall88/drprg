rule preprocessing:
    input:
        run_dir=rules.download_data.output.outdir,
        run_info=rules.aggregate_run_info.output.run_info,
    output:
        fastq=RESULTS / "preprocessing/{tech}/{proj}/{sample}/{run}.fq.gz",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
    log:
        LOGS / "preprocessing/{tech}/{proj}/{sample}/{run}.log",
    conda:
        str(ENVS / "preprocessing.yaml")
    params:
        opts=lambda wildcards: "-I -l 30 --cut_tail --dedup --stdout"
        if wildcards.tech == "illumina"
        else "-q 7",
        script=SCRIPTS / "preprocessing.sh",
    shadow:
        "shallow"
    shell:
        """
        bash {params.script} -r {wildcards.run} -i {input.run_info} \
            -o {output.fastq} -t {threads} {params.opts} 2> {log}
        """


rule build_decontamination_db:
    output:
        fasta=RESOURCES / "decontamination/remove_contam.fa.gz",
        metadata=RESOURCES / "decontamination/remove_contam.tsv",
    params:
        script=SCRIPTS / "download_tb_reference_files.pl",
        outdir=lambda wildcards, output: Path(output.fasta).parent,
    container:
        CONTAINERS["clockwork"]
    log:
        LOGS / "build_decontamination_db.log",
    shadow:
        "shallow"
    shell:
        "perl {params.script} {params.outdir} &> {log}"


BWA_EXTNS = [".amb", ".ann", ".bwt", ".pac", ".sa"]


rule index_decontam_db_with_bwa:
    input:
        fasta=rules.build_decontamination_db.output.fasta,
    output:
        index=multiext(
            str(RESOURCES / "decontamination/remove_contam.fa.gz"), *BWA_EXTNS
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
    log:
        LOGS / "index_decontam_db_with_bwa.log",
    container:
        CONTAINERS["bwa"]
    shell:
        "bwa index {input.fasta} 2> {log}"


rule index_decontam_db_with_minimap2:
    input:
        fasta=rules.build_decontamination_db.output.fasta,
    output:
        index=RESOURCES / "decontamination/remove_contam.fa.gz.map-ont.mmi",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
    log:
        LOGS / "index_decontam_db_with_minimap2.log",
    threads: 4
    container:
        CONTAINERS["minimap2"]
    params:
        opts="-x map-ont -I 12G",
    shell:
        "minimap2 {params.opts} -t {threads} -d {output.index} {input.fasta} 2> {log}"


rule map_to_decontam_db:
    input:
        bwa_index=rules.index_decontam_db_with_bwa.output.index,
        mm2_index=rules.index_decontam_db_with_minimap2.output.index,
        ref=rules.index_decontam_db_with_bwa.input.fasta,
        reads=rules.preprocessing.output.fastq,
        run_info=rules.aggregate_run_info.output.run_info,
    output:
        bam=RESULTS / "mapped/{tech}/{proj}/{sample}/{run}.sorted.bam",
        index=RESULTS / "mapped/{tech}/{proj}/{sample}/{run}.sorted.bam.bai",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
    params:
        opts=lambda wildcards: "-a -x map-ont"
        if wildcards.tech == "nanopore"
        else "-I -M",
        script=SCRIPTS / "map_to_decontam_db.sh",
        ref=lambda wildcards, input: input.ref if wildcards.tech == "illumina" else input.mm2_index,
    conda:
        str(ENVS / "aln_tools.yaml")
    log:
        LOGS / "map_to_decontam_db/{tech}/{proj}/{sample}/{run}.log",
    shell:
        """
        bash {params.script} -r {wildcards.run} -i {input.run_info} -R {input.reads} \
            -o {output.bam} -d {params.ref} -t {threads} {params.opts} 2> {log}
        """


rule filter_contamination:
    input:
        bam=rules.map_to_decontam_db.output.bam,
        metadata=rules.build_decontamination_db.output.metadata,
    output:
        keep_ids=RESULTS / "filtered/{tech}/{proj}/{sample}/{run}/keep.reads",
        contam_ids=RESULTS / "filtered/{tech}/{proj}/{sample}/{run}/contaminant.reads",
        unmapped_ids=RESULTS / "filtered/{tech}/{proj}/{sample}/{run}/unmapped.reads",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    conda:
        str(ENVS / "filter.yaml")
    params:
        script=SCRIPTS / "filter_contamination.py",
        extra="--ignore-secondary",
        outdir=lambda wildcards, output: Path(output.keep_ids).parent,
    log:
        LOGS / "filter_contamination/{tech}/{proj}/{sample}/{run}.log",
    shell:
        """
        python {params.script} {params.extra} \
            -i {input.bam} \
            -m {input.metadata} \
            -o {params.outdir} 2> {log}
        """


rule extract_decontaminated_reads:
    input:
        reads=rules.map_to_decontam_db.input.reads,
        read_ids=rules.filter_contamination.output.keep_ids,
    output:
        reads=RESULTS / "filtered/{tech}/{proj}/{sample}/{run}/{run}.filtered.fq.gz",
        stats=RESULTS / "filtered/{tech}/{proj}/{sample}/{run}/{run}.filtered.stats.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(8 * GB) * attempt,
    log:
        LOGS / "extract_decontaminated_reads/{tech}/{proj}/{sample}/{run}.log",
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        seqkit grep -o {output.reads} -f {input.read_ids} {input.reads} 2> {log}
        seqkit stats -a -T {output.reads} > {output.stats} 2>> {log}
        """

def infer_stats_files(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for run, row in df.iterrows():
        proj = row["bioproject"]
        sample = row["biosample"]
        p = (
            RESULTS
            /  f"filtered/{wildcards.tech}/{proj}/{sample}/{run}/{run}.filtered.stats.tsv",
        )
        files.append(p)

    return files

def infer_keep_files(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for run, row in df.iterrows():
        proj = row["bioproject"]
        sample = row["biosample"]
        p = (
            RESULTS
            /  f"filtered/{wildcards.tech}/{proj}/{sample}/{run}/keep.reads",
        )
        files.append(p)

    return files

def infer_contam_files(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for run, row in df.iterrows():
        proj = row["bioproject"]
        sample = row["biosample"]
        p = (
            RESULTS
            /  f"filtered/{wildcards.tech}/{proj}/{sample}/{run}/contaminant.reads",
        )
        files.append(p)

    return files

def infer_unmapped_files(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for run, row in df.iterrows():
        proj = row["bioproject"]
        sample = row["biosample"]
        p = (
            RESULTS
            /  f"filtered/{wildcards.tech}/{proj}/{sample}/{run}/unmapped.reads",
        )
        files.append(p)

    return files

rule qc_summary:
    input:
        stats=infer_stats_files,
        keep_ids=infer_keep_files,
        contam_ids=infer_contam_files,
        unmapped_ids=infer_unmapped_files,
    output:
        summary=RESULTS / "QC/{tech}.qc.csv",
    log:
        LOGS / "qc_summary/{tech}.log",
    params:
        genome_size=4411532,
    script:
        str(SCRIPTS / "qc_summary.py")
