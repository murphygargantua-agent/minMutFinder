#!/usr/bin/env bash
# tests/test_run_multi.sh — Unit tests for run_multi.sh
#
# Tests find_pair() matching conventions, --list/--dry-run flags, and the
# merge output structure. Does not require the real bioinformatics tools;
# nextflow is mocked.

set -euo pipefail
cd "$(dirname "$0")/.."

export PATH="$(pwd)/tests/mock:$PATH"

SCRIPT="$(pwd)/minMutFinder-nf/run_multi.sh"
TMP="${TMPDIR:-/tmp}/mmtest_$$"
mkdir -p "$TMP"

PASS=0
FAIL=0
assert() {
  local name="$1" expected="$2" actual="$3"
  if [[ "$expected" == "$actual" ]]; then
    echo "PASS: $name"
    PASS=$((PASS+1))
  else
    echo "FAIL: $name"
    echo "  expected: $expected"
    echo "  actual:   $actual"
    FAIL=$((FAIL+1))
  fi
}

# Extract the R1 or R2 filename from a sample line like:
#   "  S1: R1='/tmp/.../S1_R1.fq'  R2='/tmp/.../S1_R2.fq'"
get_r1() { echo "$1" | grep -oE '[A-Za-z0-9_.+-]+_R1\.fq' | tail -1; }
get_r2() { echo "$1" | grep -oE '[A-Za-z0-9_.+-]+_R2\.fq' | tail -1; }

# ---- Setup: create sample dirs with various FASTQ naming conventions ----
mkdir -p "$TMP/convention_r1r2/S1" "$TMP/convention_r1r2/S2"
echo "R1" > "$TMP/convention_r1r2/S1/S1_R1.fq"
echo "R2" > "$TMP/convention_r1r2/S1/S1_R2.fq"
echo "R1" > "$TMP/convention_r1r2/S2/S2_R1.fq"
echo "R2" > "$TMP/convention_r1r2/S2/S2_R2.fq"

mkdir -p "$TMP/convention_r1only/S1" "$TMP/convention_r1only/S2"
echo "R1" > "$TMP/convention_r1only/S1/R1.fq"
echo "R2" > "$TMP/convention_r1only/S1/R2.fq"
echo "R1" > "$TMP/convention_r1only/S2/R1.fq"
echo "R2" > "$TMP/convention_r1only/S2/R2.fq"

mkdir -p "$TMP/convention_forrev/S1" "$TMP/convention_forrev/S2"
echo "R1" > "$TMP/convention_forrev/S1/S1_for.fq"
echo "R2" > "$TMP/convention_forrev/S1/S1_rev.fq"
echo "R1" > "$TMP/convention_forrev/S2/S2_for.fq"
echo "R2" > "$TMP/convention_forrev/S2/S2_rev.fq"

# ---- Test 1: default _R1/_R2 convention via --list ----
out=$("$SCRIPT" --sample_dir "$TMP/convention_r1r2" --ref_seq ref.fasta --list 2>&1)
line1=$(echo "$out" | grep '  S1:' | head -1)
line2=$(echo "$out" | grep '  S2:' | head -1)
assert "default _R1/_R2: S1 R1" "S1_R1.fq" "$(get_r1 "$line1")"
assert "default _R1/_R2: S1 R2" "S1_R2.fq" "$(get_r2 "$line1")"
assert "default _R1/_R2: S2 R1" "S2_R1.fq" "$(get_r1 "$line2")"
assert "default _R1/_R2: S2 R2" "S2_R2.fq" "$(get_r2 "$line2")"

# ---- Test 2: R1/R2 convention (no stem) via --list ----
out=$("$SCRIPT" --sample_dir "$TMP/convention_r1only" --ref_seq ref.fasta --list 2>&1)
line1=$(echo "$out" | grep '  S1:' | head -1)
assert "R1/R2 convention: S1 R1" "R1.fq" "$(echo "$line1" | grep -oE 'R1\.fq')"
assert "R1/R2 convention: S1 R2" "R2.fq" "$(echo "$line1" | grep -oE 'R2\.fq')"

# ---- Test 3: custom --r1_pattern/--r2_pattern via --list ----
out=$("$SCRIPT" --sample_dir "$TMP/convention_forrev" --ref_seq ref.fasta \
  --r1_pattern '^.+_for\.fq$' --r2_pattern '^.+_rev\.fq$' --list 2>&1)
line1=$(echo "$out" | grep '  S1:' | head -1)
line2=$(echo "$out" | grep '  S2:' | head -1)
assert "custom for/rev: S1 R1" "S1_for.fq" "$(echo "$line1" | grep -oE 'S1_for\.fq')"
assert "custom for/rev: S1 R2" "S1_rev.fq" "$(echo "$line1" | grep -oE 'S1_rev\.fq')"
assert "custom for/rev: S2 R1" "S2_for.fq" "$(echo "$line2" | grep -oE 'S2_for\.fq')"
assert "custom for/rev: S2 R2" "S2_rev.fq" "$(echo "$line2" | grep -oE 'S2_rev\.fq')"

# ---- Test 4: --sample_ids selects specific samples ----
out=$("$SCRIPT" --sample_ids S1 --sample_dir "$TMP/convention_r1r2" --ref_seq ref.fasta --list 2>&1)
assert "--sample_ids: only S1" "1" "$(echo "$out" | grep -c '  S1:')"
assert "--sample_ids: no S2" "0" "$(echo "$out" | grep -c '  S2:')"

# ---- Test 5: --dry-run prints commands (no execution) ----
out=$("$SCRIPT" --sample_dir "$TMP/convention_r1r2" --ref_seq ref.fasta --dry-run 2>&1)
assert "--dry-run: prints nextflow" "2" "$(echo "$out" | grep -c 'nextflow')"
assert "--dry-run: does not run" "0" "$(echo "$out" | grep -c 'Sample .* finished')"

# ---- Test 6: merge produces correct CSVs (mock nextflow) ----
rm -rf "$TMP/output"
mkdir -p "$TMP/output"
"$SCRIPT" --sample_dir "$TMP/convention_r1r2" --ref_seq ref.fasta \
  --output_dir "$TMP/output" \
  --pipeline "$(pwd)/minMutFinder-nf/minMutFinder_stub.nf" --concurrency 2 2>&1 >/dev/null
assert "merge: all_mutations.csv exists" "1" "$(test -f "$TMP/output/all_mutations.csv" && echo 1 || echo 0)"
assert "merge: all_qc_metrics.csv exists" "1" "$(test -f "$TMP/output/all_qc_metrics.csv" && echo 1 || echo 0)"
mutations_header=$(head -1 "$TMP/output/all_mutations.csv")
assert "merge: mutations header" "SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency" "$mutations_header"
qc_header=$(head -1 "$TMP/output/all_qc_metrics.csv")
assert "merge: qc header" "sample,test,score" "$qc_header"
n_rows=$(wc -l < "$TMP/output/all_mutations.csv")
assert "merge: header + 2 samples rows" "5" "$n_rows"
n_headers=$(grep -c "^SampleID,Gene" "$TMP/output/all_mutations.csv")
assert "merge: no duplicate headers" "1" "$n_headers"

# ---- Test 7: recursive discovery (nested sample dirs) ----
mkdir -p "$TMP/recursive/batch1/S1" "$TMP/recursive/batch2/S2"
echo "R1" > "$TMP/recursive/batch1/S1/S1_R1.fq"
echo "R2" > "$TMP/recursive/batch1/S1/S1_R2.fq"
echo "R1" > "$TMP/recursive/batch2/S2/S2_R1.fq"
echo "R2" > "$TMP/recursive/batch2/S2/S2_R2.fq"
out=$("$SCRIPT" --sample_dir "$TMP/recursive" --ref_seq ref.fasta --recursive --list 2>&1)
line1=$(echo "$out" | grep '  S1:' | head -1)
line2=$(echo "$out" | grep '  S2:' | head -1)
assert "recursive: S1 R1" "S1_R1.fq" "$(echo "$line1" | grep -oE 'S1_R1\.fq')"
assert "recursive: S1 R2" "S1_R2.fq" "$(echo "$line1" | grep -oE 'S1_R2\.fq')"
assert "recursive: S2 R1" "S2_R1.fq" "$(echo "$line2" | grep -oE 'S2_R1\.fq')"
assert "recursive: S2 R2" "S2_R2.fq" "$(echo "$line2" | grep -oE 'S2_R2\.fq')"

# ---- Test 8: recursive merge (mock nextflow) ----
rm -rf "$TMP/routput"
mkdir -p "$TMP/routput"
"$SCRIPT" --sample_dir "$TMP/recursive" --ref_seq ref.fasta \
  --recursive --output_dir "$TMP/routput" \
  --pipeline "$(pwd)/minMutFinder-nf/minMutFinder_stub.nf" --concurrency 2 2>&1 >/dev/null
assert "recursive merge: all_mutations.csv exists" "1" "$(test -f "$TMP/routput/all_mutations.csv" && echo 1 || echo 0)"
assert "recursive merge: header + 2 samples rows" "5" "$(wc -l < "$TMP/routput/all_mutations.csv")"

# ---- Test 9: run ID prevents clobbering on repeated default runs ----
# Clean up any pre-existing default output dirs.
rm -rf "$(dirname "$SCRIPT")/../multi_run_output"* "$(dirname "$SCRIPT")/../multi_run_work"*
mkdir -p "$TMP/noclobber/S1" "$TMP/noclobber/S2"
echo "R1" > "$TMP/noclobber/S1/S1_R1.fq"
echo "R2" > "$TMP/noclobber/S1/S1_R2.fq"
echo "R1" > "$TMP/noclobber/S2/S2_R1.fq"
echo "R2" > "$TMP/noclobber/S2/S2_R2.fq"
# First run.
"$SCRIPT" --sample_dir "$TMP/noclobber" --ref_seq ref.fasta \
  --pipeline "$(pwd)/minMutFinder-nf/minMutFinder_stub.nf" --concurrency 2 --poll_interval 1 2>&1 >/dev/null
first_output=$(ls -d "$(dirname "$SCRIPT")/../multi_run_output"* 2>/dev/null | sort | head -1)
# Second run must create a DIFFERENT output dir (run ID suffix differs).
"$SCRIPT" --sample_dir "$TMP/noclobber" --ref_seq ref.fasta \
  --pipeline "$(pwd)/minMutFinder-nf/minMutFinder_stub.nf" --concurrency 2 --poll_interval 1 2>&1 >/dev/null
second_output=$(ls -d "$(dirname "$SCRIPT")/../multi_run_output"* 2>/dev/null | sort | tail -1)
assert "run ID: two separate output dirs created" "1" "$(if [[ "$first_output" != "$second_output" ]]; then echo 1; else echo 0; fi)"
# Both dirs must exist and contain comma-separated merged CSVs.
assert "run ID: first run output preserved" "1" "$(test -f "$first_output/all_mutations.csv" && echo 1 || echo 0)"
assert "run ID: second run output preserved" "1" "$(test -f "$second_output/all_mutations.csv" && echo 1 || echo 0)"
# Verify both are comma-separated (header + data rows use commas).
first_header=$(head -1 "$first_output/all_mutations.csv")
second_header=$(head -1 "$second_output/all_mutations.csv")
assert "run ID: first run header comma-separated" "SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency" "$first_header"
assert "run ID: second run header comma-separated" "SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency" "$second_header"

# ---- Test 10: error diagnostics show stderr and log tail on failure ----
# Create a sample dir with a fake nextflow that exits with an error (no completion marker).
rm -rf "$TMP/diag"
mkdir -p "$TMP/diag/S1/mutations" "$TMP/diag/S1/qc" "$TMP/diag/S1/plots"
echo "R1" > "$TMP/diag/S1/S1_R1.fq"
echo "R2" > "$TMP/diag/S1/S1_R2.fq"
# Custom fake nextflow that simulates an early failure (missing tool).
cat > "$TMP/fake_nextflow" <<'EOF'
#!/usr/bin/env bash
# Fake nextflow that simulates an early failure (missing tool).
echo "ERROR ~ [id] command not found: minimap2" >&2
echo "Command executed: minimap2 ..."
echo "Command exit status: 127"
exit 1
EOF
chmod +x "$TMP/fake_nextflow"
# Point the script's nextflow command to the fake one via NEXTFLOW_BIN.
diag_out=$(NEXTFLOW_BIN="$TMP/fake_nextflow" "$SCRIPT" --sample_dir "$TMP/diag" --ref_seq ref.fasta \
  --pipeline "$(pwd)/minMutFinder-nf/minMutFinder_stub.nf" --concurrency 1 --poll_interval 1 2>&1) || true
assert "diagnostics: stderr shown on failure" "3" "$(echo "$diag_out" | grep -c 'minimap2')"
assert "diagnostics: log status shown on failure" "1" "$(echo "$diag_out" | grep -c 'log file')"

# ---- Test 11: by_protein.nf uses shell vars, not Nextflow vars ----
# The byProteinAnalysis process loops over prot_names.txt and calls
# by_protein_analysis.py for each CDS. Using ${prot} would be interpreted
# as a Nextflow variable reference (which doesn't exist), causing a
# compilation error. It must use $prot (shell variable) instead.
bp_nf="$(cd "$(dirname "$SCRIPT")" && pwd)/modules/by_protein.nf"
prot_count=$(grep -F '${prot}' "$bp_nf" 2>/dev/null | wc -l || true)
prot_count=${prot_count:-0}
assert "by_protein.nf: no Nextflow var refs for shell loop vars" "0" "$prot_count"

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
if [[ $FAIL -gt 0 ]]; then exit 1; fi
