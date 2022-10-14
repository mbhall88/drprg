rule mykrobe_depth:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
        run_info=rules.validate_run_info.output.run_info,
    output:
        report=RESULTS / "depth/mykrobe/{depth}/{tech}/{proj}/{sample}/{run}.mykrobe.json.gz",
    log:
        LOGS / "mykrobe_depth/{depth}/{tech}/{proj}/{sample}/{run}.log",
    shadow:
        "shallow"
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
