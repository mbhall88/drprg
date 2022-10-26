

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
    benchmark:
        BENCH / "predict/mykrobe/{tech}/{proj}/{sample}/{run}.tsv"
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
        base_json=lambda wildcards, output: Path(output.report).with_suffix(""),
    threads: 2
    shell:
        """
        mykrobe predict {params.tech_opts} {params.opts} -o {params.base_json} \
            -i {input.reads} -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        gzip {params.base_json} 2>> {log}
        """


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


# ==========
# DRPRG
# ==========
rule drprg_predict:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
        index=RESULTS / f"drprg/index/w{W}/k{K}",
    output:
        report=RESULTS / "amr_predictions/drprg/{tech}/{proj}/{sample}/{run}/{run}.drprg.json",
        vcf=RESULTS / "amr_predictions/drprg/{tech}/{proj}/{sample}/{run}/{run}.drprg.bcf",
    shadow:
        "shallow"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    container:
        CONTAINERS["drprg"]
    benchmark:
        BENCH / "predict/drprg/{tech}/{proj}/{sample}/{run}.tsv"
    log:
        LOGS / "drprg_predict/{tech}/{proj}/{sample}/{run}.log",
    params:
        opts=" ".join(
            [
                "--sample {run}",
                "--verbose",
                "--ignore-synonymous",
            ]
        ),
        tech_opts=infer_drprg_tech_opts,
        filters=drprg_filter_args,
        outdir=lambda wildcards, output: Path(output.report).parent,
    threads: 2
    shell:
        """
        drprg predict {params.opts} {params.tech_opts} {params.filters} \
            -i {input.reads} -o {params.outdir} -x {input.index} -t {threads} 2> {log}
        """


rule combine_drprg_reports:
    input:
        reports=infer_drprg_reports,
    output:
        report=RESULTS / "amr_predictions/drprg/{tech}/summary.csv",
    log:
        LOGS / "combine_drprg_reports/{tech}.log",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_drprg_reports.py")


# ==========
# TBProfiler
# ==========
rule tbprofiler_predict:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
        run_info=rules.validate_run_info.output.run_info,
        db=rules.create_tbprofiler_db.output[0],
    output:
        report=RESULTS
        / "amr_predictions/tbprofiler/{tech}/{proj}/{sample}/{run}/results/{run}.results.json",
        vcf=RESULTS
        / "amr_predictions/tbprofiler/{tech}/{proj}/{sample}/{run}/vcf/{run}.targets.csq.vcf.gz",
    log:
        LOGS / "tbprofiler_predict/{tech}/{proj}/{sample}/{run}.log",
    shadow:
        "shallow"
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(6 * GB),
    conda:
        str(ENVS / "tbprofiler.yaml")
    params:
        opts="--txt --no_trim -p {run}",
        outdir=lambda wildcards, output: Path(output.report).parent.parent,
    benchmark:
        BENCH / "predict/tbprofiler/{tech}/{proj}/{sample}/{run}.tsv"
    shell:
        """
        exec 2> {log} # send all stderr from this script to the log file

        reads={input.reads}
        run_acc={wildcards.run}
        run_info={input.run_info}

        files_str=$(grep "$run_acc" "$run_info" | cut -f2)
        n_files=$(awk -F \; '{{print NF}}' <<< "$files_str")

        if [ "$n_files" -eq 2 ]; then
            tmpout=$(mktemp -d)
            prefix="${{tmpout}}/${{run_acc}}"
            # we need to deinterleave the fastq file
            seqfu deinterleave -o "$prefix" --check "$reads"
            input_arg=("-1" "${{prefix}}_R1.fq" "-2" "${{prefix}}_R2.fq")
        else
            input_arg=("-1" "$reads")
        fi

        tb-profiler profile "${{input_arg[@]}}" {params.opts} -t {threads} -d {params.outdir} --platform {wildcards.tech}

        """


rule combine_tbprofiler_reports:
    input:
        reports=infer_tbprofiler_reports,
    output:
        report=RESULTS / "amr_predictions/tbprofiler/{tech}/summary.csv",
    log:
        LOGS / "combine_tbprofiler_reports/{tech}.log",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_tbprofiler_reports.py")
