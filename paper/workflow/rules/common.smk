BWA_EXTNS = [".amb", ".ann", ".bwt", ".pac", ".sa"]


def infer_h2h_reads(wildcards, tech=None):
    sample = wildcards.sample
    tech = wildcards.tech if tech is None else tech
    if sample.startswith("mada"):
        site = "madagascar"
    elif sample.startswith("R"):
        site = "south_africa"
    else:
        site = "birmingham"

    site_dir = QC_DIR / f"subsampled/{site}/{tech}/{sample}"
    if tech == "nanopore":
        return site_dir / f"{sample}.subsampled.fastq.gz"
    elif tech == "illumina":
        return [site_dir / f"{sample}.subsampled.R{i}.fastq.gz" for i in [1, 2]]
    else:
        raise ValueError(f"Got unknown tech {tech}")


def infer_mykrobe_tech_opts(wildcards):
    return {"illumina": "-e 0.001 --ploidy haploid", "nanopore": "--ont"}[
        wildcards.tech
    ]


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


def infer_mykrobe_depth_reports(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for dp in config["depths"]:
        for run, row in df.iterrows():
            proj = row["bioproject"]
            sample = row["biosample"]
            p = (
                RESULTS
                / f"depth/mykrobe/{dp}/{wildcards.tech}/{proj}/{sample}/{run}.mykrobe.json.gz"
            )
            files.append(p)

    return files


def infer_drprg_tech_opts(wildcards) -> str:
    return {"illumina": "-I -C 4"}.get(wildcards.tech, "")


def infer_preprocessing_tech_opts(wildcards) -> str:
    return {"illumina": "-I -l 30 --cut_tail --dedup --stdout", "nanopore": "-q 7"}[
        wildcards.tech
    ]


def infer_map_opts(wildcards):
    return {"illumina": "-I -M", "nanopore": "-a -x map-ont"}[wildcards.tech]


def infer_map_ref_index(wildcards, input) -> str:
    return {"illumina": input.ref, "nanopore": input.mm2_index}[wildcards.tech]


def infer_download_dirs(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    dirs = []
    for run, row in df.iterrows():
        p = (
            RESULTS
            / f"download/{wildcards.tech}/{row.bioproject}/{row.biosample}/{run}"
        )
        dirs.append(p)

    return dirs


def infer_stats_files(wildcards):
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
            / f"filtered/{wildcards.tech}/{proj}/{sample}/{run}/{run}.filtered.stats.tsv"
        )
        files.append(str(p))

    return files


def infer_keep_files(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for run, row in df.iterrows():
        proj = row["bioproject"]
        sample = row["biosample"]
        p = RESULTS / f"filtered/{wildcards.tech}/{proj}/{sample}/{run}/keep.reads"
        files.append(str(p))

    return files


def infer_contam_files(wildcards):
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
            / f"filtered/{wildcards.tech}/{proj}/{sample}/{run}/contaminant.reads"
        )
        files.append(str(p))

    return files


def infer_unmapped_files(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for run, row in df.iterrows():
        proj = row["bioproject"]
        sample = row["biosample"]
        p = RESULTS / f"filtered/{wildcards.tech}/{proj}/{sample}/{run}/unmapped.reads"
        files.append(str(p))

    return files


def drprg_filter_args(wildcards, override_depth: int=None) -> str:
    """Generate CLI args for drprg filters"""
    filters = config.get("filters", {})

    if override_depth is not None:
        filters["min_covg"] = override_depth

    args = ""
    for (flag, key) in [
        ("-d", "min_covg"),
        ("-b", "min_strand_bias"),
        ("-g", "min_gt_conf"),
        ("-L", "max_indel"),
        ("-K", "min_frs"),
    ]:
        if key in sorted(filters):
            if key == "min_gt_conf":
                val = filters[key][wildcards.tech]
            else:
                val = filters[key]
            args += f"{flag} {val} "

    return args


def infer_drprg_reports(wildcards):
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
            / f"amr_predictions/drprg/{wildcards.tech}/{proj}/{sample}/{run}/{run}.drprg.json"
        )
        files.append(p)

    return files


def infer_drprg_depth_reports(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for dp in config["depths"]:
        for run, row in df.iterrows():
            proj = row["bioproject"]
            sample = row["biosample"]
            p = (
                RESULTS
                / f"depth/drprg/{dp}/{wildcards.tech}/{proj}/{sample}/{run}/{run}.drprg.json"
            )
            files.append(p)

    return files


def infer_tbprofiler_reports(wildcards):
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
            / f"amr_predictions/tbprofiler/{wildcards.tech}/{proj}/{sample}/{run}/results/{run}.results.json"
        )
        files.append(p)

    return files


def infer_tbprofiler_depth_reports(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for dp in config["depths"]:
        for run, row in df.iterrows():
            proj = row["bioproject"]
            sample = row["biosample"]
            p = (
                RESULTS
                / f"depth/tbprofiler/{dp}/{wildcards.tech}/{proj}/{sample}/{run}/results/{run}.results.json"
            )
            files.append(p)

    return files


def infer_benchmark_reports(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for run, row in df.iterrows():
        proj = row["bioproject"]
        sample = row["biosample"]
        for tool in TOOLS:
            p = BENCH / f"predict/{tool}/{wildcards.tech}/{proj}/{sample}/{run}.tsv"
            files.append(p)

    return files
