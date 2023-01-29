#!/usr/bin/env bash

set -euxo pipefail

exec 2> "${snakemake_log[0]}" # send all stderr from this script to the log file

reads="${snakemake_input[reads]}"
run_acc="${snakemake_wildcards[run]}"
run_info="${snakemake_input[run_info]}"
depth="${snakemake_wildcards[depth]}"
seed="${snakemake_params[seed]}"
genome_size="${snakemake_params[genome_size]}"

files_str=$(grep "$run_acc" "$run_info" | cut -f2)
n_files=$(awk -F \; '{print NF}' <<< "$files_str")
tmpout=$(mktemp -d)
trap 'rm -rf -- "$tmpout"' EXIT
prefix="${tmpout}/${run_acc}"
subreads="${prefix}_${depth}.fq"

if [ "$n_files" -eq 2 ]; then
    # we need to deinterleave the fastq file
    seqfu deinterleave -o "$prefix" --check "$reads"

    # now we subsample
    in_r1="${prefix}_R1.fq"
    out_r1="${prefix}_${depth}_R1.fq"
    in_r2="${prefix}_R2.fq"
    out_r2="${prefix}_${depth}_R2.fq"

    rasusa -i "$in_r1" "$in_r2" -o "$out_r1" "$out_r2" -s $seed -c $depth -g $genome_size

    # re-interleave the subsampled reads
    seqfu interleave --check -1 "$out_r1" -2 "$out_r2" > "$subreads"

else
    # no need to deinterleave, we can just subsample the input fastq
    rasusa -i "$reads" -o "$subreads" -s $seed -c $depth -g $genome_size
fi

drprg predict ${snakemake_params[opts]} ${snakemake_params[tech_opts]} ${snakemake_params[filters]} \
    -i "$subreads" -o ${snakemake_params[outdir]} -x "${snakemake_input[index]}" -t ${snakemake[threads]}

rm -rf "$tmpout"
