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
