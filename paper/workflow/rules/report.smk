rule compare_sn_and_sp:
    input:
        summary_files=expand(
            str(RESULTS / "amr_predictions/{tool}/{{tech}}/summary.csv"), tool=TOOLS
        ),
        phenotypes=lambda wildcards: config[f"{wildcards.tech}_samplesheet"],
    output:
        plots=report(
            multiext(str(PLOTS / "sn_sp/{tech}"), ".png", ".svg"),
            category="Sn/Sp",
            subcategory="Figure",
            labels={"Technology": "{tech}"},
        ),
        table=report(
            TABLES / "sn_sp/summary.{tech}.csv",
            category="Sn/Sp",
            subcategory="Tables",
            labels={"Technology": "{tech}", "Table": "Summary"},
        ),
        classification=report(
            TABLES / "sn_sp/classifications.{tech}.csv",
            category="Sn/Sp",
            subcategory="Tables",
            labels={"Technology": "{tech}", "Table": "Classifications"},
        ),
    log:
        LOGS / "compare_sn_and_sp/{tech}.log",
    params:
        minor_is_susceptible=False,
        unknown_is_resistant=False,
        failed_is_resistant=False,
        figsize=(13, 8),
        dpi=300,
        sn_marker="^",
        sp_marker="o",
        ignore_drugs=("ciprofloxacin", "all"),
        min_num_phenotypes=10,
    conda:
        str(ENVS / "compare_sn_and_sp.yaml")
    script:
        str(SCRIPTS / "compare_sn_and_sp.py")


rule aggregate_predict_benchmarks:
    input:
        bench=infer_benchmark_reports,
    output:
        summary=BENCH / "predict/{tech}.summary.csv",
    log:
        LOGS / "aggregate_predict_benchmarks/{tech}.log",
    container:
        CONTAINERS["python"]
    params:
        delim=",",
    script:
        SCRIPTS / "aggregate_predict_benchmarks.py"
