
![Static Badge](https://img.shields.io/badge/Version-Pre--Release-blue) ![License](https://img.shields.io/badge/License-GPL_V3-green)

# **minMutFinder**

## 📜 Table of Contents
- [🎯 Overview](#-overview)
- [🚀 Quickstart](#-quickstart)
- [🔍 Features](#-features)
- [🛠 Prerequisites](#-prerequisites)
- [📥 Installation](#-installation)
- [🚀 How to Run](#-how-to-run-minmutfinder)
- [⚙️ Arguments](#-arguments)
  - [📝 Notes](#-notes)
- [📊 MultiQC Report](#-multiqc-report)
- [🔏 License](#-license)
- [🖊️ Citing minMutFinder](#-citing-minmutfinder)
- [🔮 Future Work & Limitations](#-future-work-and-limitations)
- [✉️ Get in Touch](#-get-in-touch)
- [📚 References](#-references)

---

## 🎯 Overview

**minMutFinder** is a bioinformatics tool designed to help you identify minority mutations in population variants with precision and accuracy. Unlike other tools, **minMutFinder** considers the possibility of multiple nucleotide mutations within the same codon and provides comprehensive metrics for your sequences.

---

## 🚀 Quickstart

Get a multi-sample run going in three steps:

### Manifest mode (explicit r1/r2)

```bash
# 1. Clone the fork
git clone https://github.com/murphygargantua-agent/minMutFinder.git
cd minMutFinder

# 2. Prepare a manifest CSV
cat > samples.csv << 'EOF'
sample_name,r1,r2
SampleA,/path/to/SampleA_R1.fq.gz,/path/to/SampleA_R2.fq.gz
SampleB,/path/to/SampleB_R1.fq.gz,/path/to/SampleB_R2.fq.gz
EOF

# 3. Run (uses the real pipeline by default; pass --pipeline to use the stub)
nextflow run minMutFinder-nf/minMutFinder.nf \
  --manifest samples.csv \
  --ref_seq minMutFinder-nf/resources/cmv_cds.fasta \
  --annotate mutations.tsv \
  --syn_muts yes \
  --out_path out/manifest-test
```

### Scan mode (auto-pair by name)

```bash
# 1. Clone the fork
git clone https://github.com/murphygargantua-agent/minMutFinder.git
cd minMutFinder

# 2. Prepare a sample directory (one subdir per sample, each with paired FASTQs)
mkdir -p samples/SampleA samples/SampleB
cp /path/to/SampleA_R1.fq.gz samples/SampleA/
cp /path/to/SampleA_R2.fq.gz samples/SampleA/
cp /path/to/SampleB_R1.fq.gz samples/SampleB/
cp /path/to/SampleB_R2.fq.gz samples/SampleB/

# 3. Run — the pipeline pairs *_R1.* / *_R2.* automatically, no manifest needed
nextflow run minMutFinder-nf/minMutFinder.nf \
  --inputReads samples \
  --ref_seq minMutFinder-nf/resources/cmv_cds.fasta \
  --out_path out/scan-test
```

Outputs:
- `out/<sample>/` — per-sample pipeline outputs + plots
- `out/all_mutations.csv` — merged mutations table
- `out/all_qc_metrics.csv` — merged QC table
- `out/multiqc_report/` — consolidated MultiQC HTML report

To smoke-test the orchestration without installing bioinformatics tools, swap in the stub pipeline:

```bash
nextflow run minMutFinder-nf/minMutFinder.nf --manifest samples.csv --ref_seq minMutFinder-nf/resources/cmv_cds.fasta --pipeline minMutFinder-nf/minMutFinder_stub.nf
```

> **Note**: `run_multi.sh` (the legacy bash orchestrator) is still supported for
> single-sample runs launched in parallel, but the manifest/scan modes above run all
> samples within one Nextflow invocation and are the recommended path.

## 🔍 Features

- 🧬 **Advanced Mutation Detection**: Identifies minority mutations while accounting for multiple nucleotide changes within a single codon.
- 📊 **Comprehensive Analysis**: Provides detailed metrics and plots for a thorough understanding of your sequences.
- 🔧 **Customizable**: Supports various versions and annotated mutations for enhanced flexibility.

---

## 🛠 Prerequisites

Ensure the following programs are installed:

### Required Software
| Software      | Version | Installation |
| ------------- | ------- | ------------ |
| **Nextflow**  | 23.10.1 or higher | [Install](https://www.nextflow.io/docs/latest/install.html) |
| **Python**    | 3.9    |  |
| **libgcc-ng** | 12 or higher | [conda-forge](https://conda-forge.org/) |
| **Trimmomatic** | 0.39   | [Bioconda](https://bioconda.github.io/) |
| **Minimap2**  | 2.26    | [Bioconda](https://bioconda.github.io/) |
| **Lofreq**    | 2.1.5   | [Bioconda](https://bioconda.github.io/) |
| **Bcftools**  | 1.17    | [Bioconda](https://bioconda.github.io/) |
| **Samtools**  | 1.17    | [Bioconda](https://bioconda.github.io/) |

### Required Python Packages:
```bash
os, pandas, sys, csv, gzip, shutil, matplotlib, seaborn, Bio, re, plotly, numpy
```

---

## 📥 Installation

### Step 1: Clone the Repository
```bash
git clone https://github.com/murphygargantua-agent/minMutFinder.git
cd minMutFinder
```

### Step 2: Install Dependencies
```bash
# Install conda-forge dependencies
conda install -c conda-forge libgcc-ng>=12

# Install bioconda dependencies
conda install -c bioconda minimap2=2.26 htslib=1.17 bcftools=1.17 samtools=1.18 trimmomatic=0.39 lofreq=2.1.5=py39hb7ef6d5_10

# Install Python dependencies
pip install pandas>=2.1.4 Bio>=1.5.9 plotly>=5.18.0 numpy>=1.26.2 pysam>=0.21.0 matplotlib>=3.8.2 seaborn>=0.13.0
```

---

## 🚀 How to Run minMutFinder

```bash
nextflow run minMutFinder.nf --ref_seq <reference.fasta> --out_path <output_name> --r1 <forward_reads.fastq.gz> --r2 <reverse_reads.fastq.gz> --annotate <mutations.tsv> --syn_muts <"yes"/"no">
```

All path arguments (`--ref_seq`, `--r1`, `--r2`, `--out_path`, `--annotate`) accept
both relative and absolute paths. Relative paths are resolved against the current
working directory and converted to absolute internally, so the pipeline can be run
from any directory.

For example, from the project root:
```bash
nextflow run minMutFinder-nf/minMutFinder.nf   --ref_seq minMutFinder-nf/resources/cmv_cds.fasta   --r1 ../cmv_analysis/r1.fq   --r2 ../cmv_analysis/r2.fq   --out_path output/CMV2   --syn_muts yes
```

### Multi-sample

Run multiple samples at once. Two mutually exclusive input modes:

**1. Manifest mode** — a CSV with `sample_name`, `r1` and `r2`. Nextflow runs each
sample in parallel within a single pipeline run and merges the per-sample results
into master tables:

```bash
nextflow run minMutFinder-nf/minMutFinder.nf \
  --manifest samples.csv \
  --ref_seq minMutFinder-nf/resources/cmv_cds.fasta \
  --annotate mutations.tsv \
  --syn_muts yes \
  --out_path out/manifest-test
```

The manifest CSV (`samples.csv`) must have a header row and the columns:
`sample_name`, `r1`, `r2`. Paths may be relative (resolved against the CWD) or
absolute.

```bash
// samples.csv
sample_name,r1,r2
SampleA,/path/to/SampleA_R1.fq.gz,/path/to/SampleA_R2.fq.gz
SampleB,/path/to/SampleB_R1.fq.gz,/path/to/SampleB_R2.fq.gz
```

**2. Scan mode** — point at a directory of paired FASTQ and let the pipeline auto-pair
forward/reverse files by name (like `run_multi.sh`'s `find_pair`, or tpal's
`Channel.fromFilePairs`). No manifest needed:

```bash
nextflow run minMutFinder-nf/minMutFinder.nf \
  --inputReads samples \
  --ref_seq minMutFinder-nf/resources/cmv_cds.fasta \
  --out_path out/scan-test
```

Expect one subdirectory per sample, each containing `{sample}_R1.fq.gz` and
`{sample}_R2.fq.gz`. The naming convention matches `run_multi.sh`'s `find_pair`:

- `*_R1.*` + `*_R2.*` (any extension: `.fastq`, `.fq`, `.fq.gz`)
- files may sit directly in `--inputReads` or one level down inside per-sample dirs

Outputs are merged into `all_mutations.csv` and `all_qc_metrics.csv` under `--out_path`:

A consolidated MultiQC QC report is generated under `<out_path>/multiqc_report/`:
it merges the per-sample `qc/qc_metrics.csv` tables, `fastp` JSON logs and
`bcftools stats` reports into one HTML report.

### MultiQC report

Every run generates a consolidated MultiQC quality-control report:

- **Manifest mode** (one Nextflow run, all samples): the report merges per-sample
  `qc/qc_metrics.csv` tables, `fastp` JSON logs and `bcftools stats` reports into a
  single `multiqc_report.html` under `<out_path>/multiqc_report/`.
- **Single-sample mode** (`--r1/--r2`, or the per-sample runs launched by
  `run_multi.sh`): one report per sample, under `<out_dir>/<sample>/multiqc_report/`.

The per-sample QC metrics CSVs are converted to MultiQC custom-content format with
`bin/qc_metrics_mqc.py` (wide table: `sample,reads_total,reads_filt_qc,...`), and the
fastp/samtools/bcftools logs are produced inside `runMultiQC`
(`modules/multiqc.nf`). The report configuration is in
`minMutFinder-nf/multiqc_config.yaml`.

### Synthetic test data

`minMutFinder-nf/resources/synthetic/synth_R1.fq` + `synth_R2.fq` are small
deterministic paired-end reads (2,000 read pairs, 150 bp) generated with `wgsim`
from a single CDS contig of the CMV AD169 reference. They exercise the full
end-to-end chain — trimmomatic → minimap2 → lofreq → bcftools → byProteinAnalysis
→ MultiQC — on real data, and are intended for smoke-testing the pipeline without
requiring a large dataset:

```bash
nextflow run minMutFinder-nf/minMutFinder.nf \
  --manifest <(printf 'sample_name,r1,r2\nsynthetic,minMutFinder-nf/resources/synthetic/synth_R1.fq,minMutFinder-nf/resources/synthetic/synth_R2.fq') \
  --ref_seq minMutFinder-nf/resources/cmv_cds.fasta \
  --out_path out/synthetic-test
```

The synthetic reads carry known SNVs at RL1 positions 683, 700, 806 and 882
(e.g. `C683G`, `C700A`, `C806G`, `G882C`), which the pipeline detects at ~50-100%
allele frequency, so they also validate variant calling against ground truth.

### Smoke test: scan mode

Copy the synthetic FASTQs into a directory and run with `--inputReads` — the
pipeline auto-pairs them without a manifest:

```bash
mkdir -p out/synthetic-scan && cp minMutFinder-nf/resources/synthetic/*.fq out/synthetic-scan/
nextflow run minMutFinder-nf/minMutFinder.nf \
  --inputReads out/synthetic-scan \
  --ref_seq minMutFinder-nf/resources/cmv_cds.fasta \
  --out_path out/synthetic-scan-test
```

Expect `synth_R1.fq` + `synth_R2.fq` to be paired as sample `synth`.

### Custom FASTQ naming

If your files use a different convention (e.g. `_for`/`_rev`, `L001_R1`/`L001_R2`),
pass ERE patterns with `--r1_pattern` and `--r2_pattern`:

```bash
./run_multi.sh \
  --sample_dir samples \
  --ref_seq ref.fasta \
  --annotate mutations.tsv \
  --r1_pattern '^.+_for\.fq\.gz$' \
  --r2_pattern '^.+_rev\.fq\.gz$' \
  --concurrency 4
```

Before running, use `--list` to verify which files get detected:

```bash
./run_multi.sh --sample_dir samples --ref_seq ref.fasta \
  --r1_pattern '^.+_for\.fq\.gz$' --r2_pattern '^.+_rev\.fq\.gz$' --list
```

The default `_R1/_R2`, `r1/r2`, and `R1/R2` conventions continue to work without
any extra arguments. If your sample directories are nested (e.g. `batch1/S1`,
`batch2/S2`), add `--recursive` to search all subdirectories:

```bash
./run_multi.sh --sample_dir samples --ref_seq ref.fasta --annotate mutations.tsv \
  --recursive --concurrency 4
```

A test suite in `tests/` validates the pairing/merge logic (including recursive
discovery) without needing the real bioinformatics tools.

### Validation status

This multi-CDS pipeline and multi-sample orchestration were validated end-to-end:

- **Multi-CDS processing** (`byProteinAnalysis`): verified with the CMV2 sample
  (`r1.fq` + `r2.fq`, 2,696 read pairs) against the AD169 CMV GenBank CDS reference
  (`resources/cmv_cds.fasta`, 8 drug-target CDSs: UL54, UL56, UL89, UL91, UL92,
  UL93, UL94, UL97). The pipeline processes ALL 8 CDSs, producing per-CDS
  mutation CSVs, per-CDS QC metrics, and a merged `CMV2_mutations.csv`. Known CMV
  drug-resistance mutations were detected (e.g., UL54 T885A/A1108T/T1122A at AF=1.0;
  UL97 T75A/Q126L).

- **Wrapper orchestration** (`run_multi.sh`): verified end-to-end with the stub
  pipeline — per-sample runs execute in parallel (limited by `--concurrency`),
  per-sample output directories are created with the expected structure, and the
  merge step correctly produces `all_mutations.csv` and `all_qc_metrics.csv`.

- **Single-sample pipeline** (`minMutFinder.nf`): the explicit-arg refactor makes
  the DAG compile and instantiate in Nextflow 26.04.6 (the previous
  auto-flattening-dependent DAG no longer works in this Nextflow version).

- **Bug fixes**: the `csv.reader(delimiter="\n")` incompatibility with newer
  Python was replaced with manual tab-splitting; the `samfile_list` loop was fixed
  to use the pre-split list; the `byProteinAnalysis` process now iterates over all
  CDSs in `prot_names.txt` via a bash loop (restoring the original upstream per-CDS
  design).
---

## ⚙️ Arguments

- `--ref_seq`: Path and filename of the reference genome FASTA file (1)(2). Accepts relative or absolute paths.
- `--out_path`: Output name for the virus column. Accepts relative or absolute paths.
- `--manifest`: Path to a manifest CSV with `sample_name`, `r1` and `r2` columns (multi-sample mode).
- `--inputReads`: Path to a directory of paired FASTQ; the pipeline auto-pairs forward/reverse files by name (alternative multi-sample mode, mutually exclusive with `--manifest`).
- `--r1`: Path and filename of the forward FASTQ compressed file. Accepts relative or absolute paths.
- `--r2`: Path and filename of the reverse FASTQ compressed file. Accepts relative or absolute paths.
- `--annotate`: Path and filename of the TSV file containing the annotated mutations (3). Accepts relative or absolute paths.
- `--syn_muts`: "yes" or "no", depending on whether to include synonymous mutations in the output plot (default is "no")

### 📝 Notes

  1. The reference genome must contain the coding sequences (CDS) of the proteins. If there are multiple proteins, they should be separated in the FASTA file. A pre-built AD169 CMV CDS reference is included at `minMutFinder-nf/resources/cmv_cds.fasta`.
  2. FASTA headers must use underscores (`_`) between words, or pipes (`|`) if
     the reference is a multi-CDS GenBank FASTA. For example:
     `>NC_006273_2_UL96` or `>UL97|NC_006273.2`.
  3. The annotated mutation file should be tab-separated and contain a column named `mutation` for annotated mutations.
  4. If FASTA headers share patterns (e.g. H3N2_PA and H3N2_PA_X share H3N2_PA) they must be differentiated (e.g. H3N2_PA_prot, H3N2_PA_X).

---

## 🔏 License

This project, **minMutFinder** (fork), is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this software under the terms of this license. For more details, refer to the [LICENSE](./LICENSE) file.

---

## 🖊️ Citing minMutFinder

A research paper on **minMutFinder** is currently in progress. In the meantime, please cite this fork's GitHub repository. The original minMutFinder was developed by Ignasi Prats-Méndez at the Respiratory Viruses Unit, Vall d'Hebron Hospital Universitari & VHIR; this fork adds multi-CDS support and multi-sample orchestration.

---

## 🔮 Future Work and Limitations

### Limitations:
- As of now, **minMutFinder** has only been tested on viral sequencing data.
- Only available for Illumina® sequencing data.

### Current Thresholds:
- **Allele Frequency (AF) ≥ 5%**
- **Read depth per nucleotide position ≥ 20**
- **Number of threads = 64**

### Future Improvements:
- At the moment we are working on uploading **minMutFinder** to nf-core, for easier distribution and use.
- Support for user-provided VCF files along with SAM/BAM files, skipping the initial quality control, mapping, and variant calling steps.
- Make it possible for the user to choose the number of threads, AF, Read depth per nucleotide positon.
- Make it possible for the user to choose kind of output desired (complete or just final result files).
- Support for ONT® sequencing data. 

---

## ✉️ Get in Touch

If you encounter any issues, have feature requests, or need assistance, feel free to reach out:

- **Open an issue** directly in this repository by clicking [here](https://github.com/murphygargantua-agent/minMutFinder/issues).
- **Email us** at [ignasi.prats@vhir.org](mailto:ignasi.prats@vhir.org).

We're always happy to help!

---

## 📚 References

- [Nextflow](https://www.nextflow.io/)
- [Python](https://www.python.org/)
- [Bioconda](https://bioconda.github.io/)
- [Minimap2](https://github.com/lh3/minimap2)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [Lofreq](https://csb5.github.io/lofreq/)
- [Bcftools](http://samtools.github.io/bcftools/bcftools.html)
- [Samtools](http://www.htslib.org/)
