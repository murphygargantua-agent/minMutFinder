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
assert "merge: mutations header" "SampleID;Gene;Mutation_type;Aa_change;Type_of_aa_change;Nt_mutation;Mutation_frequency" "$mutations_header"
qc_header=$(head -1 "$TMP/output/all_qc_metrics.csv")
assert "merge: qc header" "sample;test;score" "$qc_header"
n_rows=$(wc -l < "$TMP/output/all_mutations.csv")
assert "merge: header + 2 samples rows" "5" "$n_rows"

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
if [[ $FAIL -gt 0 ]]; then exit 1; fi
