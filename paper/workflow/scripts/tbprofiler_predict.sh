#!/usr/bin/env bash
set -euxo pipefail

# shellcheck disable=SC2154
exec 2> "${snakemake_log[0]}" # send all stderr from this script to the log file

# shellcheck disable=SC2154
reads="${snakemake_input[reads]}"
# shellcheck disable=SC2154
run_acc="${snakemake_wildcards[run]}"
run_info="${snakemake_input[run_info]}"

files_str=$(grep "$run_acc" "$run_info" | cut -f2)
IFS=';' read -r -a files <<< "$files_str"
n_files="${#files[@]}"

tmpout=$(mktemp -d)
trap 'rm -rf -- "$tmpout"' EXIT

if [ "$n_files" -eq 2 ]; then
    prefix="${tmpout}/${run_acc}"
    # we need to deinterleave the fastq file
    seqfu deinterleave -o "$prefix" --check "$reads"
    input_arg=("-1" "${prefix}_R1.fq" "-2" "${prefix}_R2.fq")
else
    input_arg=("-1" "$reads")
fi

# shellcheck disable=SC2154
IFS=" " read -r -a opts <<< "${snakemake_params[opts]}"

if [ "${snakemake_wildcards[tech]}" = "nanopore" ]; then
    opts+=("--platform" "nanopore")
fi

# shellcheck disable=SC2154
tb-profiler profile "${input_arg[@]}" "${opts[@]}" -t "${snakemake[threads]}" -d "${snakemake_params[outdir]}"

rm -rf "$tmpout"