#!/usr/bin/env python3
"""Convert a minMutFinder per-sample QC metrics CSV into MultiQC custom-content format.

minMutFinder writes one row per QC test (long format):

    sample,test,score
    S1,reads_total,1000
    S1,reads_filt_qc,800

MultiQC custom content expects one row per sample (wide format) preceded by a YAML
comment header that configures the report section:

    # id: minmutfinder_qc
    # plot_type: table
    # section_name: minMutFinder QC metrics
    # description: Read counts from the minMutFinder QC step
    sample,reads_total,reads_filt_qc,reads_filt_qc_paired
    S1,1000,800,750

Each sample's converted file must be named ``<sample>_mqc.csv`` so that MultiQC's
custom content file search pattern (``.+_mqc\\.(csv|...)``) finds it. Files with
other names are ignored by the custom content module.

Usage:
    python3 qc_metrics_mqc.py --in qc_metrics.csv --out S1_mqc.csv --sample S1
"""

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in", dest="in_csv", required=True, help="Per-sample QC metrics CSV (long format)")
    parser.add_argument("--out", dest="out_csv", required=True, help="Output MultiQC custom-content CSV")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    # Read the long-format CSV
    rows = []
    with open(args.in_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row.get("sample"):
                continue
            rows.append(row)

    if not rows:
        raise SystemExit(f"No data found in {args.in_csv}")

    # Collect all test names (columns) from the data, preserving order
    cols = []
    for row in rows:
        key = row.get("test")
        if key and key not in cols:
            cols.append(key)

    # Write the wide-format CSV with the MultiQC custom-content config header.
    # Each input row (one per test) is pivoted into one column, so a sample that
    # had 3 tests is written as a single row of 3 values.
    with open(args.out_csv, "w", newline="") as f:
        f.write("# id: minmutfinder_qc\n")
        f.write("# plot_type: table\n")
        f.write("# section_name: minMutFinder QC metrics\n")
        f.write("# description: Read counts from the minMutFinder QC step\n")
        header = ["sample"] + cols
        f.write(",".join(header) + "\n")
        # Collect the value for each test, then write the sample row
        metrics = {r.get("test", ""): r.get("score", "") for r in rows if r.get("test")}
        vals = [args.sample] + [metrics.get(c, "") for c in cols]
        f.write(",".join(vals) + "\n")


if __name__ == "__main__":
    main()
