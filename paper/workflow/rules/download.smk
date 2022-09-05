rule download_data:
    output:
        outdir=directory(RESULTS / "download/{tech}/{proj}/{sample}/{run}"),
        run_info=RESULTS / "download/{tech}/{proj}/{sample}/{run}/fastq-run-info.json",
    log:
        LOGS / "download_data/{tech}/{proj}/{sample}/{run}.log",
    container:
        CONTAINERS["fastq_dl"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(0.5 * GB),
    params:
        db="sra",
    shell:
        "fastq-dl --outdir {output.outdir} {wildcards.run} {params.db} > {log} 2>&1"




rule aggregate_run_info:
    input:
        dirs=infer_download_dirs,
    output:
        run_info=RESULTS / "download/{tech}/run_info.tsv",
    log:
        LOGS / "aggregate_run_info/{tech}.log",
    params:
        delim="\t",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "aggregate_run_info.py")
