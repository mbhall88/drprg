#!/usr/bin/env bash
set -xeo pipefail

DEFAULT_THREADS=1
IS_ILLUMINA=false

function usage {
    echo "usage: map_to_decontam.sh -r <STR> -i <FILE> -R <FILE> -o <FILE> -d <FILE> [OPTIONS]"
    echo "   "
    echo "  -r                       : Run accession [REQUIRED]"
    echo "  -i                       : Run info TSV file [REQUIRED]"
    echo "  -R                       : Input reads (fastq) file [REQUIRED]"
    echo "  -o                       : Output BAM file [REQUIRED]"
    echo "  -d                       : Database to map to [REQUIRED]"
    echo "  -I                       : Input is Illumina"
    echo "  -t                       : Number of threads [default: $DEFAULT_THREADS]"
    echo "  -h | --help              : This message"
}

function parse_args {
    # positional args
    args=()

    # named args
    while [ "$1" != "" ]; do
        case "$1" in
            -r)
                run_acc="$2"
                shift
                ;;
            -i)
                run_info="$2"
                shift
                ;;
            -R)
                reads="$2"
                shift
                ;;
            -o)
                output="$2"
                shift
                ;;
            -d)
                db="$2"
                shift
                ;;
            -t)
                threads="$2"
                shift
                ;;
            -I)
                IS_ILLUMINA=true
                ;;
            -h | --help)
                usage
                exit
                ;;             # quit and show usage
            *) args+=("$1") ;; # if no match, add it to the positional args
        esac
        shift # move to next kv pair
    done

    # restore positional args
    #    set -- "${args[@]}"

    # set positionals to vars
    #  positional_1="${args[0]}"
    #  positional_2="${args[1]}"

    # validate required args
    if [[ -z "${run_acc}" || -z "${run_info}" || -z "${reads}" || -z "${output}" || -z "$db" ]]; then
        echo "Invalid arguments"
        usage
        exit
    fi

    # set defaults
    if [[ -z "$threads" ]]; then
        threads="$DEFAULT_THREADS"
    fi
}

function run {
    parse_args "$@"

    files_str=$(grep "$run_acc" "$run_info" | cut -f2)
    IFS=';' read -r -a files <<< "$files_str"
    n_files="${#files[@]}"

    if [ "$IS_ILLUMINA" = true ]; then
        if [ "$n_files" -eq 2 ]; then
            args+=("-p")
        fi

        bwa mem ${args[*]} -t "$threads" "$db" "$reads" |
            samtools sort -n -@ "$threads" |
            samtools fixmate -m -@ "$threads" - - |
            samtools sort -@ "$threads" |
            samtools markdup -r -S -O bam - "$output"

    else
        minimap2 ${args[*]} -t "$threads" "$db" "$reads" |
            samtools sort -@ "$threads" -o "$output"
    fi

    samtools index -b -@ "$threads" "$output"
}

run "$@"
