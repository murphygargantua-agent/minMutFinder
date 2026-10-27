#!/usr/bin/env bash
# run_multi.sh — Multi-sample orchestrator for the minMutFinder fork.
#
# Runs the existing single-sample minMutFinder pipeline once per sample (in
# parallel, concurrency-limited) and merges per-sample results into master
# tables. This reuses the tested single-sample pipeline verbatim; no internal
# Nextflow changes are required.
#
# Usage:
#   ./run_multi.sh \
#     --sample_dir DIR \
#     [--sample_ids s1,s2,s3] \
#     --ref_seq PATH \
#     [--annotate mutations.tsv] \
#     [--syn_muts yes|no] \
#     [--output_dir OUT] \
#     [--work_dir WORK] \
#     [--concurrency N]
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NF="${PIPELINE_DIR}/minMutFinder.nf"
CONFIG="${PIPELINE_DIR}/nextflow_resources.config"

# Allow overriding config for stubs / alternate pipelines
CONFIG_OVERRIDDEN=false

die() { echo "ERROR: $*" >&2; exit 1; }

[[ -f "$NF" ]] || die "pipeline not found: $NF"
command -v nextflow >/dev/null 2>&1 || die "nextflow not on PATH; install Nextflow (23.10.1+) or add it to PATH"

sample_dir=
sample_ids=
ref_seq=
annotate=
syn_muts=no
output_dir="${PIPELINE_DIR}/../multi_run_output"
work_dir="${PIPELINE_DIR}/../multi_run_work"
concurrency=2
dry_run=false
pipeline=""

while [[ $# -gt 0 ]]; do
  usage() {
  cat <<'EOF'
Usage: ./run_multi.sh --sample_dir DIR --ref_seq PATH [options]

Run minMutFinder on every sample found under --sample_dir (one subdirectory per
sample) and merge the per-sample results.

Options:
  --sample_dir     DIR with one subdirectory per sample, each containing a paired
                   FASTQ set ({sample}_R1.fq.gz + {sample}_R2.fq.gz, or r1.fq/r2.fq)
  --sample_ids     Optional comma-separated list of sample subdirs to run
                   (defaults to all subdirectories)
  --ref_seq        Absolute path to the CDS reference FASTA
  --annotate       Optional absolute path to the mutations annotation TSV
  --syn_muts       "yes" or "no" (default no)
  --output_dir     Where per-sample output dirs go (default: ../multi_run_output)
  --work_dir       Where per-sample Nextflow work dirs go (default: ../multi_run_work)
  --concurrency    Max parallel nextflow runs (default 2)
  --pipeline       Path to an alternate Nextflow pipeline (default: minMutFinder.nf)
  --config         Path to an alternate Nextflow config (default: nextflow_resources.config)
EOF
}


case "$1" in
    --sample_dir)  sample_dir="$2"; shift 2 ;;
    --sample_ids)  sample_ids="$2"; shift 2 ;;
    --ref_seq)     ref_seq="$2"; shift 2 ;;
    --annotate)    annotate="$2"; shift 2 ;;
    --syn_muts)    syn_muts="$2"; shift 2 ;;
    --output_dir)  output_dir="$2"; shift 2 ;;
    --work_dir)    work_dir="$2"; shift 2 ;;
    --concurrency) concurrency="$2"; shift 2 ;;
    --dry-run)     dry_run=true; shift 1 ;;
    --config)      CONFIG="$2"; CONFIG_OVERRIDDEN=true; shift 2 ;;
    --pipeline)    pipeline="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

if [[ -n "$pipeline" ]]; then NF="$pipeline"; fi

[[ -d "$sample_dir" ]] || die "sample_dir does not exist: $sample_dir"

# Resolve to absolute paths (pipeline requires absolute --out_path)
sample_dir="$(cd "$sample_dir" && pwd)"
output_dir="$(mkdir -p "$output_dir" && cd "$output_dir" && pwd)"
work_dir="$(mkdir -p "$work_dir" && cd "$work_dir" && pwd)"
ref_seq="$(cd "$(dirname "$ref_seq")" && pwd)/$(basename "$ref_seq")"
annotate_abs=""
if [[ -n "$annotate" ]]; then
  annotate_abs="$(cd "$(dirname "$annotate")" && pwd)/$(basename "$annotate")"
fi

# Build sample list
declare -a samples=()
if [[ -n "$sample_ids" ]]; then
  IFS=',' read -ra samples <<< "$sample_ids"
else
  for d in "$sample_dir"/*/; do
    [[ -d "$d" ]] || continue
    samples+=("$(basename "$d")")
  done
fi
[[ ${#samples[@]} -gt 0 ]] || die "no samples found in $sample_dir"

echo "=== run_multi.sh: ${#samples[@]} sample(s), concurrency=$concurrency ==="
echo "  sample_dir=$sample_dir  output_dir=$output_dir  work_dir=$work_dir"
if [[ "$dry_run" == true ]]; then
  echo "[DRY RUN] Would run the following per-sample pipelines:"
  for s in "${samples[@]}"; do
    echo "  nextflow -w '$work_dir/$s' run '$NF' -c '$CONFIG' \\"
    echo "    --ref_seq '$ref_seq' --r1 <pair> --r2 <pair> --out_path '$output_dir/$s' --syn_muts '$syn_muts' ${annotate:+--annotate '$annotate_abs'}"
  done
  exit 0
fi

# Find a paired FASTQ (R1/R2) inside a sample directory.
# Supports: {stem}_R1.fq(.gz) / {stem}_R2.fq(.gz), r1.fq/r2.fq, R1.fq/R2.fq.
find_pair() {
  local sdir="$1" r1="" r2="" stem=""
  # _R1/_R2 convention (stem = prefix before _R1); allow multi-extension like .fq.gz
  for f in "$sdir"/*; do
    local base="${f##*/}"
    if [[ "$base" =~ ^(.+)_R1(\.[a-zA-Z0-9]+)*$ ]]; then
      stem="${BASH_REMATCH[1]}"
      r1="$f"
    fi
  done
  if [[ -n "$r1" ]]; then
    local esc="${stem//./\\.}"
    for f in "$sdir"/*; do
      local base="${f##*/}"
      if [[ "$base" =~ ^${esc}_R2(\.[a-zA-Z0-9]+)*$ ]]; then r2="$f"; break; fi
    done
  fi
  if [[ -z "$r2" ]]; then
    # r1/r2 or R1/R2 convention (exact name prefix, case-insensitive)
    for f in "$sdir"/*; do
      local base="${f##*/}"
      case "$base" in
        r1*|R1*) r1="$f" ;;
        r2*|R2*) r2="$f" ;;
      esac
    done
  fi
  printf '%s %s' "$r1" "$r2"
}

run_sample() {
  local sample="$1"
  local sdir="$sample_dir/$sample"
  local pair
  pair=$(find_pair "$sdir")
  local r1 r2
  r1=$(echo "$pair" | awk '{print $1}')
  r2=$(echo "$pair" | awk '{print $2}')
  if [[ -z "$r1" || -z "$r2" || ! -f "$r1" || ! -f "$r2" ]]; then
    echo "ERROR: no valid paired FASTQ in $sdir (found: '$r1' '$r2')" >&2
    return 1
  fi
  local out_path="$output_dir/$sample"
  local wdir="$work_dir/$sample"
  mkdir -p "$wdir"
  echo "=== [$(date +%T)] Running sample: $sample ==="
  echo "  R1=$r1  R2=$r2  out=$out_path  work=$wdir"
  nextflow run "$NF" \
    -w "$wdir" \
    -c "$CONFIG" \
    --ref_seq "$ref_seq" \
    --r1 "$r1" \
    --r2 "$r2" \
    --out_path "$out_path" \
    --syn_muts "$syn_muts" \
    ${annotate:+--annotate "$annotate_abs"} \
    2>&1 | tee "$wdir/${sample}.log"
  local rc=${PIPESTATUS[0]}
  echo "=== [$(date +%T)] Sample $sample finished (rc=$rc) ==="
  return $rc
}
export -f run_sample find_pair
export sample_dir output_dir work_dir NF CONFIG ref_seq syn_muts annotate_abs concurrency

# Run samples in parallel, limited to --concurrency.
printf '%s\0' "${samples[@]}" | xargs -0 -P "$concurrency" -I {} bash -c 'run_sample "$@"' _ {}
xargs_rc=${PIPESTATUS[0]}

echo "=== Merging per-sample results into $output_dir ==="

HEADER="SampleID;Gene;Mutation_type;Aa_change;Type_of_aa_change;Nt_mutation;Mutation_frequency"
echo "$HEADER" > "$output_dir/all_mutations.csv"
for s in "${samples[@]}"; do
  f="$output_dir/$s/mutations/${s}_mutations.csv"
  if [[ -f "$f" ]]; then
    # Each per-sample file already has a header; skip it when concatenating.
    tail -n +2 "$f" >> "$output_dir/all_mutations.csv"
  fi
done

echo "sample;test;score" > "$output_dir/all_qc_metrics.csv"
for s in "${samples[@]}"; do
  f="$output_dir/$s/qc/qc_metrics.csv"
  if [[ -f "$f" ]]; then
    tail -n +2 "$f" >> "$output_dir/all_qc_metrics.csv"
  fi
done

echo ""
echo "=== Done. Per-sample outputs under $output_dir/<sample>/ ==="
echo "  merged: $output_dir/all_mutations.csv"
echo "  merged: $output_dir/all_qc_metrics.csv"
echo "  per-sample plots remain in $output_dir/<sample>/plots/"
exit $xargs_rc
