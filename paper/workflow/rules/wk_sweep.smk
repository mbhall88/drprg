"""These rules are for sweeping through pandora values for w and k for each technology
to see what works best for drprg.
"""
WK_SWEEP = RESULTS / "drprg/wk_sweep"


rule seqtk_mergepe:
    input:
        reads=lambda wildcards: infer_h2h_reads(wildcards, tech="illumina"),
    output:
        merged=WK_SWEEP / "reads/illumina/{sample}.fq.gz",
    params:
        compress_lvl=9,
    log:
        LOGS / "seqtk_mergepe/{sample}.log",
    threads: 2
    wrapper:
        "0.75.0-7-g05edf56/bio/seqtk/mergepe"


rule copy_nanopore_reads:
    input:
        reads=lambda wildcards: infer_h2h_reads(wildcards, tech="nanopore"),
    output:
        reads=WK_SWEEP / "reads/nanopore/{sample}.fq.gz",
    log:
        LOGS / "copy_nanopore_reads/{sample}.log",
    resources:
        mem_mb=int(0.3 * GB),
    shell:
        "cp {input.reads} {output.reads} 2> {log}"


rule drprg_predict_sweep:
    input:
        index=rules.drprg_build.output.outdir,
        reads=WK_SWEEP / "reads/{tech}/{sample}.fq.gz",
    output:
        outdir=directory(WK_SWEEP / "predict/w{w}/k{k}/{tech}/{sample}"),
        report=WK_SWEEP / "predict/w{w}/k{k}/{tech}/{sample}/{sample}.drprg.json",
        vcf=WK_SWEEP / "predict/w{w}/k{k}/{tech}/{sample}/{sample}.drprg.bcf",
    log:
        LOGS / "drprg_predict_sweep/w{w}/k{k}/{tech}/{sample}.log",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    container:
        CONTAINERS["drprg"]
    params:
        opts=" ".join(["-v", "-s {sample}", "-u", "--failed"]),
        filters=lambda wildcards: drprg_filter_args(wildcards),
        tech_flag=lambda wildcards: "-I" if wildcards.tech == "illumina" else "",
    shell:
        """
        drprg predict {params.opts} {params.filters} {params.tech_flag} \
          -o {output.outdir} -i {input.reads} -x {input.index} -t {threads} 2> {log}
        """