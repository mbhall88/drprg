rule download_data:
    output:
        outdir=directory(RESULTS / "download/{tech}/{proj}/{sample}/{run}"),
        run_info=RESULTS / "download/{tech}/{proj}/{sample}/{run}/fastq-run-info.json",
    log:
        LOGS / "download_data/{tech}/{proj}/{sample}/{run}.log",
    container:
        CONTAINERS["fastq_dl"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(2 * GB),
    params:
        db="sra",
    shadow:
        "shallow"
    shell:
        "fastq-dl --outdir {output.outdir} {wildcards.run} {params.db} > {log} 2>&1"


rule validate_run_info:
    input:
        info=rules.download_data.output.run_info,
    output:
        run_info=RESULTS / "validate/{tech}/{proj}/{sample}/{run}/run_info.tsv",
    log:
        LOGS / "validate_run_info/{tech}/{proj}/{sample}/{run}.log",
    params:
        delim="\t",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "validate_run_info.py")
