def infer_mykrobe_tech_opts(wildcards):
    if wildcards.tech == "illumina":
        return "-e 0.001 --ploidy haploid"
    else:
        return "--ont"


rule mykrobe_predict:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
    output:
        report=RESULTS
        / "amr_predictions/mykrobe/{tech}/{proj}/{sample}/{run}.mykrobe.json.gz",
    shadow:
        "shallow"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    container:
        CONTAINERS["mykrobe"]
    log:
        LOGS / "mykrobe_predict/{tech}/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-A",
                "-O json",
                "-D 0.20",
                "--species tb",
                "--sample {run}",
            ]
        ),
        tech_opts=infer_mykrobe_tech_opts,
        base_json=lambda wildcards, output: Path(output.report).with_suffix("")
    threads: 2
    shell:
        """
        mykrobe predict {params.tech_opts} {params.opts} -o {params.base_json} \
            -i {input.reads} -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        gzip {params.base_json} 2>> {log}
        """


def infer_mykrobe_reports(wildcards):
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
            / f"amr_predictions/mykrobe/{wildcards.tech}/{proj}/{sample}/{run}.mykrobe.json.gz"
        )
        files.append(p)

    return files


rule combine_mykrobe_reports:
    input:
        reports=infer_mykrobe_reports,
    output:
        report=RESULTS / "amr_predictions/mykrobe/{tech}/summary.csv",
    log:
        LOGS / "combine_mykrobe_reports/{tech}.log",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_mykrobe_reports.py")
