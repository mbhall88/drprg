rule compare_sn_and_sp:
    input:
        summary_files=expand(
            str(RESULTS / "amr_predictions/{tool}/{{tech}}/summary.csv"), tool=TOOLS
        ),
        phenotypes=lambda wildcards: config[f"{wildcards.tech}_samplesheet"],
    output:
        plots=multiext(str(PLOTS / "sn_sp/{tech}"), ".png", ".svg"),
        table=TABLES / "sn_sp/summary.{tech}.csv",
        classification=TABLES / "sn_sp/classifications.{tech}.csv",
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
        ignore_drugs={"ciprofloxacin"},
    conda:
        str(ENVS / "compare_sn_and_sp.yaml")
    script:
        str(SCRIPTS / "compare_sn_and_sp.py")
