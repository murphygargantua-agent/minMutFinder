
![Static Badge](https://img.shields.io/badge/Version-Pre--Release-blue) ![License](https://img.shields.io/badge/License-GPL_V3-green)

# **minMutFinder**

## 📜 Table of Contents
- [🎯 Overview](#-overview)
- [🔍 Features](#-features)
- [🛠 Prerequisites](#-prerequisites)
- [📥 Installation](#-installation)
- [🚀 How to Run](#-how-to-run-minmutfinder)
- [⚙️ Arguments](#-arguments)
  - [📝 Notes](#-notes)
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

# 3. Run (uses the real pipeline by default; pass --pipeline to use the stub)
./minMutFinder-nf/run_multi.sh \
  --sample_dir samples \
  --ref_seq minMutFinder-nf/resources/cmv_cds.fasta \
  --annotate mutations.tsv \
  --syn_muts yes \
  --concurrency 2
```

Outputs:
- `multi_run_output/<sample>/` — per-sample pipeline outputs + plots
- `multi_run_output/all_mutations.csv` — merged mutations table
- `multi_run_output/all_qc_metrics.csv` — merged QC table

To smoke-test the orchestration without installing bioinformatics tools, swap in the stub pipeline:

```bash
./minMutFinder-nf/run_multi.sh \
  --sample_dir samples --ref_seq minMutFinder-nf/resources/cmv_cds.fasta \
  --pipeline minMutFinder-nf/minMutFinder_stub.nf --concurrency 2 --dry-run
```

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

Run multiple samples at once with `run_multi.sh`, which reuses the tested
single-sample pipeline in parallel (concurrency-limited) and merges per-sample
results into master tables:

```bash
./run_multi.sh \
  --sample_dir <dir_with_one_subdir_per_sample> \
  --ref_seq <reference.fasta> \
  --annotate <mutations.tsv> \
  --syn_muts yes \
  --concurrency 2
```

Each sample subdirectory should contain a paired FASTQ set
(`{sample}_R1.fq.gz` + `{sample}_R2.fq.gz`, or `r1.fq`/`r2.fq`).
Per-sample outputs land in `<output_dir>/<sample>/`; merged tables
`all_mutations.csv` and `all_qc_metrics.csv` are written to `<output_dir>`.
Per-sample plots remain under each `<output_dir>/<sample>/plots/`.

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

- `--ref_seq`: Path and filename of the reference genome FASTA file (1)(2)(4). Accepts relative or absolute paths.
- `--out_path`: Output name for the virus column. Accepts relative or absolute paths.
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
