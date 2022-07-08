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


def drprg_filter_args(wildcards) -> str:
    """Generate CLI args for drprg filters"""
    filters = config.get("filters", {})
    args = ""
    for (flag, key) in [
        ("-d", "min_covg"),
        ("-b", "min_strand_bias"),
        ("-g", "min_gt_conf"),
        ("-L", "max_indel"),
        ("-K", "min_frs"),
    ]:
        if key in filters:
            if key == "min_gt_conf":
                val = filters[key][wildcards.tech]
            else:
                val = filters[key]
            args += f"{flag} {val} "

    return args