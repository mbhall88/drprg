#!/usr/bin/env bash
set -eux

if ! command -v snakemake &>/dev/null; then
  echo "snakemake could not be found"
  exit 1
fi

JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')
LOG_DIR=logs/

if [[ ! -d "$LOG_DIR" ]]; then
  echo "Error: Log directory $LOG_DIR does not exist"
  exit 1
fi

MEMORY=16000
THREADS=2
PROFILE="lsf"
BINDS="/tmp,$HOME"
BINDS+=",/hps/scratch,/hps/nobackup/iqbal,/nfs/research/zi,$FASTSW_DIR --scratch /hps/scratch"

ARGS="--contain -B $BINDS"

SNAKEMAKE_TMP="tmp/snakemake"
mkdir -p "$SNAKEMAKE_TMP"

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
  -M "$MEMORY" \
  -n "$THREADS" \
  -o "$LOG_DIR"/"$JOB_NAME".o \
  -e "$LOG_DIR"/"$JOB_NAME".e \
  -J "$JOB_NAME" \
  snakemake --profile "$PROFILE" \
  --default-resources tmpdir=$SNAKEMAKE_TMP \
  --local-cores "$THREADS" \
  "$@" --singularity-args "$ARGS"

exit 0
