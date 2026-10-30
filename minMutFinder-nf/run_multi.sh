#!/usr/bin/env bash
# run_multi.sh — Multi-sample orchestrator for the minMutFinder fork.
#
# Runs the existing single-sample minMutFinder pipeline once per sample (in
# parallel, concurrency-limited) and merges per-sample results into master
# tables. Each nextflow run is launched in the background with `-bg` so the
# Nextflow TUI is detached from the terminal; the script polls per-sample log
# files for a compact status readout.
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
#     [--concurrency N] \
#     [--r1_pattern REGEX] [--r2_pattern REGEX] \
#     [--dry-run] [--list] [--recursive]
#
# Options:
#   --poll_interval    Seconds between status polls (default 5)
#   --resume           Resume previous Nextflow runs from their work directory

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NF="${PIPELINE_DIR}/minMutFinder.nf"
CONFIG="${PIPELINE_DIR}/nextflow_resources.config"

CONFIG_OVERRIDDEN=false

die() { echo "ERROR: $*" >&2; exit 1; }

[[ -f "$NF" ]] || die "pipeline not found: $NF"
# Use NEXTFLOW_BIN if set (allows tests/override with a fake nextflow).
NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"
command -v "$NEXTFLOW_BIN" >/dev/null 2>&1 || die "nextflow not on PATH; install Nextflow (26.04.0+) or add it to PATH"

# Check Nextflow version (we require >= 26.04.0).
_nf_version=$("${NEXTFLOW_BIN}" version 2>/dev/null | grep -oE 'version [0-9]+\.[0-9]+\.[0-9]+' | head -1 | awk '{print $2}') || true
if [[ -n "$_nf_version" ]]; then
  _nf_major=$(echo "$_nf_version" | cut -d. -f1)
  _nf_minor=$(echo "$_nf_version" | cut -d. -f2)
  if [[ "$_nf_major" -lt 26 ]] || { [[ "$_nf_major" -eq 26 ]] && [[ "$_nf_minor" -lt 4 ]]; }; then
    echo "WARNING: Nextflow $_nf_version is below the minimum supported version (26.04.0). Some features may not work." >&2
  fi
fi
unset _nf_version _nf_major _nf_minor

sample_dir=
sample_ids=
ref_seq=
annotate=
syn_muts=no
output_dir="${PIPELINE_DIR}/../multi_run_output"
work_dir="${PIPELINE_DIR}/../multi_run_work"
# Track whether defaults were used so we can append a run ID for namespacing.
default_output=true
default_work=true
# Optional: skip network filesystem checks (e.g. for testing on network drives).
no_lock_check=false
# Generate a unique run ID so each run's merged outputs are namespaced
# and don't overwrite previous runs.
run_id="$(date +%Y%m%d_%H%M%S)_$$"
concurrency=2
dry_run=false
list_only=false
recursive=false
pipeline=""
r1_pattern=""
r2_pattern=""
poll_interval=5
resume=false

while [[ $# -gt 0 ]]; do
  usage() {
  cat <<'EOF'
Usage: ./run_multi.sh --sample_dir DIR --ref_seq PATH [options]

Run minMutFinder on every sample found under --sample_dir (one subdirectory per
sample) and merge the per-sample results. Each nextflow run is backgrounded with
`-bg` so the Nextflow TUI is detached; a compact status line is polled until all
jobs finish.

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
  --poll_interval  Seconds between status polls (default 5)
  --pipeline       Path to an alternate Nextflow pipeline (default: minMutFinder.nf)
  --config         Path to an alternate Nextflow config (default: nextflow_resources.config)
  --r1_pattern     ERE pattern to match the R1 (forward) FASTQ within each sample dir
  --r2_pattern     ERE pattern to match the R2 (reverse) FASTQ within each sample dir
  --dry-run        Print the per-sample nextflow commands without executing them
  --list           Detect and print the paired FASTQ for each sample, then exit
  --recursive      Search for sample directories recursively (not just one level deep)
  --no_lock_check  Skip the network filesystem lock check
EOF
}

  case "$1" in
    --sample_dir)  sample_dir="$2"; shift 2 ;;
    --sample_ids)  sample_ids="$2"; shift 2 ;;
    --ref_seq)     ref_seq="$2"; shift 2 ;;
    --annotate)    annotate="$2"; shift 2 ;;
    --syn_muts)    syn_muts="$2"; shift 2 ;;
    --output_dir)  output_dir="$2"; default_output=false; shift 2 ;;
    --work_dir)    work_dir="$2"; default_work=false; shift 2 ;;
    --concurrency) concurrency="$2"; shift 2 ;;
    --poll_interval) poll_interval="$2"; shift 2 ;;
    --resume) resume=true; shift 1 ;;
    --dry-run)     dry_run=true; shift 1 ;;
    --config)      CONFIG="$2"; CONFIG_OVERRIDDEN=true; shift 2 ;;
    --pipeline)    pipeline="$2"; shift 2 ;;
    --r1_pattern)  r1_pattern="$2"; shift 2 ;;
    --r2_pattern)  r2_pattern="$2"; shift 2 ;;
    --list)        list_only=true; shift 1 ;;
    --recursive)   recursive=true; shift 1 ;;
    --no_lock_check) no_lock_check=true; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

if [[ -n "$pipeline" ]]; then NF="$pipeline"; fi

# Find a paired FASTQ (R1/R2) inside a sample directory.
find_pair() {
  local sdir="$1" r1="" r2="" stem=""
  if [[ -n "$r1_pattern" && -n "$r2_pattern" ]]; then
    for f in "$sdir"/*; do
      local base="${f##*/}"
      if [[ "$base" =~ $r1_pattern ]]; then r1="$f"; fi
      if [[ "$base" =~ $r2_pattern ]]; then r2="$f"; fi
    done
  else
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
      for f in "$sdir"/*; do
        local base="${f##*/}"
        case "$base" in
          r1*|R1*) r1="$f" ;;
          r2*|R2*) r2="$f" ;;
        esac
      done
    fi
  fi
  printf '%s %s' "$r1" "$r2"
}

[[ -d "$sample_dir" ]] || die "sample_dir does not exist: $sample_dir"

sample_dir="$(cd "$sample_dir" && pwd)"

# Append run ID to default output/work dirs so each run's merged outputs are
# namespaced and don't overwrite previous runs.
if [[ "$default_output" == "true" ]]; then
  output_dir="${output_dir}_${run_id}"
fi
if [[ "$default_work" == "true" ]]; then
  work_dir="${work_dir}_${run_id}"
fi

output_dir="$(mkdir -p "$output_dir" && cd "$output_dir" && pwd)"
work_dir="$(mkdir -p "$work_dir" && cd "$work_dir" && pwd)"
ref_seq="$(cd "$(dirname "$ref_seq")" && pwd)/$(basename "$ref_seq")"
annotate_abs=""
if [[ -n "$annotate" ]]; then
  annotate_abs="$(cd "$(dirname "$annotate")" && pwd)/$(basename "$annotate")"
fi

# Optional: skip network filesystem lock checks (e.g. for testing on network drives).
# Nextflow creates .nextflow in cwd and uses file locks; network drives often fail.
if [[ "$no_lock_check" != true ]]; then
  # Check if the work directory is on a mount that supports file locks.
  if [[ "$default_work" == "true" ]]; then
    _work_type=$(df -T "$work_dir" 2>/dev/null | tail -1 | awk '{print $2}')
    if [[ "$_work_type" =~ ^nfs|smb|cifs|gluster|webdav ]]; then
      echo "WARNING: work_dir is on a network filesystem ($_work_type). Nextflow needs file locks; use a local --work_dir /path/to/local." >&2
    fi
  fi
  unset _work_type

  # Warn if the current working directory is on a network filesystem.
  _cwd=$(pwd)
  _cwd_type=$(df -T "$_cwd" 2>/dev/null | tail -1 | awk '{print $2}')
  if [[ "$_cwd_type" =~ ^nfs|smb|cifs|gluster|webdav ]]; then
    echo "WARNING: current working directory is on a network filesystem ($_cwd_type)." >&2
    echo "  Nextflow creates .nextflow in cwd; network drives don't support file locks." >&2
    echo "  cd to a local directory and re-run, or use --no_lock_check to override." >&2
  fi
  unset _cwd _cwd_type
fi

declare -a samples=()
if [[ -n "$sample_ids" ]]; then
  IFS=',' read -ra samples <<< "$sample_ids"
else
  while IFS= read -r d; do
    [[ -z "$d" ]] && continue
    rel="${d#$sample_dir/}"
    rel="${rel%/}"
    pair=$(find_pair "$d")
    r1=$(echo "$pair" | awk '{print $1}')
    r2=$(echo "$pair" | awk '{print $2}')
    if [[ -n "$r1" && -n "$r2" && -f "$r1" && -f "$r2" ]]; then
      samples+=("$rel")
    fi
  done < <(
    if [[ "$recursive" == true ]]; then
      find "$sample_dir" -mindepth 1 -type d -print0 | sort -z | tr '\0' '\n'
    else
      for d in "$sample_dir"/*/; do echo "$d"; done | sort
    fi
  )
fi
[[ ${#samples[@]} -gt 0 ]] || die "no samples found in $sample_dir"

if [[ "$list_only" == true ]]; then
  echo "=== Detected FASTQ pairs per sample ==="
  for s in "${samples[@]}"; do
    pair=$(find_pair "$sample_dir/$s")
    r1=$(echo "$pair" | awk '{print $1}')
    r2=$(echo "$pair" | awk '{print $2}')
    echo "  $(basename "$s"): R1='$r1'  R2='$r2'"
  done
  exit 0
fi

echo "=== run_multi.sh: ${#samples[@]} sample(s), concurrency=$concurrency, poll_interval=$poll_interval ==="
echo "  sample_dir=$sample_dir  output_dir=$output_dir  work_dir=$work_dir"
if [[ "$dry_run" == true ]]; then
  echo "[DRY RUN] Would run the following per-sample pipelines:"
  for s in "${samples[@]}"; do
    bname=$(basename "$s")
    resume_opt=""
    if [[ "$resume" == true ]]; then resume_opt=" --resume"; fi
    echo "  ${NEXTFLOW_BIN} run -bg$resume_opt '$NF' -c '$CONFIG' --log '$work_dir/$bname/$bname.log' \\"
    echo "    --ref_seq '$ref_seq' --r1 <pair> --r2 <pair> --out_path '$output_dir/$bname' --syn_muts '$syn_muts' ${annotate:+--annotate '$annotate_abs'}"
  done
  exit 0
fi

# ---------------------------------------------------------------------------
# Backgrounded per-sample runner using `nextflow run -bg`.
# ---------------------------------------------------------------------------
run_sample_bg() {
  local sample_path="$1"
  local sample_name="${sample_path##*/}"
  local sdir="$sample_dir/$sample_path"
  local pair
  pair=$(find_pair "$sdir")
  local r1 r2
  r1=$(echo "$pair" | awk '{print $1}')
  r2=$(echo "$pair" | awk '{print $2}')
  if [[ -z "$r1" || -z "$r2" || ! -f "$r1" || ! -f "$r2" ]]; then
    echo "ERROR: no valid paired FASTQ in $sdir (found: '$r1' '$r2')" >&2
    return 1
  fi
  local out_path="$output_dir/$sample_name"
  local wdir="$work_dir/$sample_name"
  mkdir -p "$wdir"

  local log="$wdir/${sample_name}.log"
  local status_file="$wdir/${sample_name}.status"

  echo "=== [$(date +%T)] Launching sample: $sample_name ==="
  echo "  R1=$r1  R2=$r2  out=$out_path  work=$wdir"

  # Run in background, detached from terminal; log to a per-sample file.
  # If --resume is set and the work dir exists from a prior run, resume it.
  local resume_flag=""
  if [[ "$resume" == true && -d "$wdir" ]]; then
    resume_flag="--resume"
  fi
  ${NEXTFLOW_BIN} run -bg $resume_flag "$NF" \
    -w "$wdir" \
    -c "$CONFIG" \
    --ref_seq "$ref_seq" \
    --r1 "$r1" \
    --r2 "$r2" \
    --out_path "$out_path" \
    --syn_muts "$syn_muts" \
    ${annotate:+--annotate "$annotate_abs"} \
    --log "$log" \
    --name "minmutfinder_$sample_name" \
    > "$wdir/${sample_name}.stderr" 2>&1

  # Capture the launched PID for later polling.
  _pid=$(ps -eo pid,args | grep -i "minmutfinder_$sample_name" | grep -v grep | head -1 | awk '{print $1}')
  echo "  launched pid=$_pid"
  
  # If resuming, clear any previous status marker so we poll the resumed run.
  if [[ "$resume" == true && -f "$status_file" ]]; then
    rm -f "$status_file"
  fi
  echo "  status=running"
}
export -f run_sample_bg find_pair
export sample_dir output_dir work_dir NF CONFIG ref_seq syn_muts annotate_abs poll_interval r1_pattern r2_pattern resume

# Launch samples in parallel (concurrency-limited).
for i in "${!samples[@]}"; do
  if [[ $((i % concurrency)) -eq 0 ]]; then
    # Wait for a free slot before launching the next one.
    wait
  fi
  run_sample_bg "${samples[$i]}" &
done
wait

echo ""
echo "=== Polling samples for completion ==="
all_rc=0
declare -a statuses=()
for s in "${samples[@]}"; do
  bname=$(basename "$s")
  status_file="$work_dir/$bname/${bname}.status"
  log="$work_dir/$bname/${bname}.log"
  statuses+=("running")
done

# Poll until all are done.
all_done=true
while "$all_done"; do
  all_done=false
  for s in "${samples[@]}"; do
    bname=$(basename "$s")
    log="$work_dir/$bname/${bname}.log"
    status_file="$work_dir/$bname/${bname}.status"

    if [[ -f "$status_file" ]]; then
      continue  # already marked done
    fi

    # Check for completion markers in the log.
    if [[ -f "$log" ]] && grep -qi "Execution completed" "$log" 2>/dev/null; then
      echo "  $bname: completed"
      : > "$status_file"
    elif [[ -f "$log" ]] && grep -qi "Execution [fe]ailed\|Execution error\|Error" "$log" 2>/dev/null; then
      echo "  $bname: FAILED"
      echo "  --- $bname log tail (last 30 lines) ---"
      tail -n 30 "$log"
      echo "  --- end $bname log ---"
      echo "FAILED" > "$status_file"
      all_rc=1
    elif ps -eo args | grep -i "minmutfinder_$bname" | grep -v grep >/dev/null 2>&1; then
      : # still running
    else
      # Process exited but no completion marker in log (e.g. mock nextflow exits
      # immediately, or Nextflow failed before writing to the log).
      if [[ -d "$output_dir/$bname" ]]; then
        echo "  $bname: completed (pre-check)"
        : > "$status_file"
      else
        echo "  $bname: FAILED (no completion marker)"
        echo "  (Nextflow process exited without 'Execution completed' or 'Execution failed' in the log)"
        # Show the Nextflow stderr if available — it often contains the real error.
        stderr_file="$work_dir/$bname/${bname}.stderr"
        if [[ -f "$stderr_file" ]]; then
          if [[ -s "$stderr_file" ]]; then
            echo "  --- $bname stderr (last 30 lines) ---"
            tail -n 30 "$stderr_file"
            echo "  --- end $bname stderr ---"
          else
            echo "  $bname stderr file is empty."
          fi
        else
          echo "  $bname stderr file not found."
        fi
        # Show the Nextflow run log tail if available.
        if [[ -f "$log" ]]; then
          if [[ -s "$log" ]]; then
            echo "  --- $bname log tail (last 30 lines) ---"
            tail -n 30 "$log"
            echo "  --- end $bname log ---"
          else
            echo "  $bname log file is empty."
          fi
        else
          echo "  $bname log file not found (Nextflow may have failed to start)."
        fi
        echo "HINT: If Nextflow reports missing tools (minimap2, bcftools, lofreq, trimmomatic),"
        echo "  install them from bioconda. See README for details."
        echo "FAILED" > "$status_file"
        all_rc=1
      fi
    fi
  done
  if [[ "$all_done" == true ]]; then
    break
  fi
  sleep "$poll_interval"
done

echo ""
echo "=== Merging per-sample results into $output_dir ==="

HEADER="SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency"
echo "$HEADER" > "$output_dir/all_mutations.csv"
for s in "${samples[@]}"; do
  bname=$(basename "$s")
  f="$output_dir/$bname/mutations/${bname}_mutations.csv"
  if [[ -f "$f" ]]; then
    tail -n +2 "$f" >> "$output_dir/all_mutations.csv"
  fi
done

echo "sample,test,score" > "$output_dir/all_qc_metrics.csv"
for s in "${samples[@]}"; do
  bname=$(basename "$s")
  f="$output_dir/$bname/qc/qc_metrics.csv"
  if [[ -f "$f" ]]; then
    tail -n +2 "$f" >> "$output_dir/all_qc_metrics.csv"
  fi
done

echo ""
echo "=== Done. Per-sample outputs under $output_dir/<sample>/ ==="
echo "  merged: $output_dir/all_mutations.csv"
echo "  merged: $output_dir/all_qc_metrics.csv"
echo "  per-sample logs: $work_dir/*.log"
echo "  per-sample plots remain in $output_dir/<sample>/plots/"
exit $all_rc