rule mykrobe_depth:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
        run_info=rules.validate_run_info.output.run_info,
    output:
        report=RESULTS
        / "depth/mykrobe/{depth}/{tech}/{proj}/{sample}/{run}.mykrobe.json.gz",
    log:
        LOGS / "mykrobe_depth/{depth}/{tech}/{proj}/{sample}/{run}.log",
    shadow:
        "shallow"
    resources:
        mem_mb=int(4 * GB),
    params:
        seed=88,
        genome_size=config["genome_size"],
        mykrobe_opts=rules.mykrobe_predict.params.opts,
        tech_opts=infer_mykrobe_tech_opts,
    conda:
        ENVS / "mykrobe.yaml"
    script:
        SCRIPTS / "mykrobe_depth.sh"


rule combine_mykrobe_depth_reports:
    input:
        reports=infer_mykrobe_depth_reports,
    output:
        report=RESULTS / "depth/mykrobe/{tech}.summary.csv",
    log:
        LOGS / "combine_mykrobe_reports/{tech}.log",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_mykrobe_depth_reports.py")


rule tbprofiler_depth:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
        run_info=rules.validate_run_info.output.run_info,
        db=rules.create_tbprofiler_db.output[0],
    output:
        report=RESULTS
        / "depth/tbprofiler/{depth}/{tech}/{proj}/{sample}/{run}/results/{run}.results.json",
        vcf=RESULTS
        / "depth/tbprofiler/{depth}/{tech}/{proj}/{sample}/{run}/vcf/{run}.targets.csq.vcf.gz",
    log:
        LOGS / "tbprofiler_depth/{depth}/{tech}/{proj}/{sample}/{run}.log",
    shadow:
        "shallow"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(6 * GB),
    conda:
        str(ENVS / "tbprofiler.yaml")
    params:
        opts="--txt --no_trim -p {run}",
        min_depth=lambda wildcards: 10 if wildcards.depth > 20 else 1,
        outdir=lambda wildcards, output: Path(output.report).parent.parent,
        seed=rules.mykrobe_depth.params.seed,
        genome_size=rules.mykrobe_depth.params.genome_size,
    script:
        SCRIPTS / "tbprofiler_depth.sh"


rule combine_tbprofiler_depth_reports:
    input:
        reports=infer_tbprofiler_depth_reports,
    output:
        report=RESULTS / "depth/tbprofiler/{tech}.summary.csv",
    log:
        LOGS / "combine_tbprofiler_reports/{tech}.log",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_tbprofiler_depth_reports.py")


rule drprg_depth:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
        index=RESULTS / f"drprg/index/w{W}/k{K}",
        run_info=rules.validate_run_info.output.run_info,
    output:
        report=RESULTS
        / "depth/drprg/{depth}/{tech}/{proj}/{sample}/{run}/{run}.drprg.json",
        vcf=RESULTS / "depth/drprg/{depth}/{tech}/{proj}/{sample}/{run}/{run}.drprg.bcf",
    shadow:
        "shallow"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    container:
        CONTAINERS["drprg"]
    log:
        LOGS / "drprg_depth/{depth}/{tech}/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--sample {run}",
                "--verbose",
                "--failed",
                "--ignore-synonymous",
            ]
        ),
        tech_opts=infer_drprg_tech_opts,
        filters=lambda wildcards: drprg_filter_args(wildcards, None if wildcards.depth > 20 else 1),
        outdir=lambda wildcards, output: Path(output.report).parent,
        seed=rules.mykrobe_depth.params.seed,
        genome_size=rules.mykrobe_depth.params.genome_size,
    threads: 2
    script:
        SCRIPTS / "drprg_depth.sh"


rule combine_drprg_depth_reports:
    input:
        reports=infer_drprg_depth_reports,
    output:
        report=RESULTS / "depth/drprg/{tech}.summary.csv",
    log:
        LOGS / "combine_drprg_depth_reports/{tech}.log",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_drprg_depth_reports.py")
