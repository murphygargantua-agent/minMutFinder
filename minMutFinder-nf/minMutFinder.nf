#!/usr/bin/env nextflow

// This file is part of minMutFinder.
//
// minMutFinder is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// minMutFinder is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with minMutFinder. If not, see <https://www.gnu.org/licenses/>.
//
// Copyright (C) 2024 Respiratory Viruses Unit, Microbiology Department, Vall d’Hebron Hospital Universitari, Vall d’Hebron Institut de Recerca (VHIR), Vall d’Hebron Barcelona Hospital Campus, Passeig Vall d’Hebron 119-129, 08035 Barcelona, Spain

include { dirCreator; refCheck } from './modules/initialization'
include { fastqProcessing; mapping1; variant_calling; read_depth; mapping2; reads_qc_plot } from './modules/from_fq_to_vcf'
include { byProteinAnalysis; protNames } from './modules/by_protein'
include { waitForOutMuts; vizNoAnnot; annotateAndViz } from './modules/viz_and_annotate'

def versionMessage() 
{
	log.info"""
	 
	MINORITY MUTATION FINDER (minMutFinder) - Version: ${workflow.manifest.version} 
	""".stripIndent()
}

def helpMessage() 
{
	log.info"""

    minMutFinder v1.0.0 under GPL-3.0 license

Author: Ignasi Prats-Méndez
Supervisor: Alejandra González-Sánchez  
Institution: HUVH & VHIR  
Group: Servei de Microbiologia - Unitat de Virus Respiratoris  


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


## 🎯 Overview

**minMutFinder** is a bioinformatics tool designed to help you identify minority mutations in 
population variants with precision and accuracy. Unlike other tools, **minMutFinder** considers 
the possibility of multiple nucleotide mutations within the same codon and provides comprehensive 
metrics for your sequences.

---

## 🔍 Features

- 🧬 **Advanced Mutation Detection**: Identifies minority mutations while accounting for multiple 
nucleotide changes within a single codon.
- 📊 **Comprehensive Analysis**: Provides detailed metrics and plots for a thorough understanding 
of your sequences.
- 🔧 **Customizable**: Supports various versions and annotated mutations for enhanced flexibility.

---

## ⚙️ Arguments

- `--ref_seq`: Path and filename of the reference genome FASTA file (1)(2). Accepts relative or absolute paths.
- `--out_path`: Output name for the virus column. Accepts relative or absolute paths.
- `--r1`: Path and filename of the forward FASTQ compressed file. Accepts relative or absolute paths.
- `--r2`: Path and filename of the reverse FASTQ compressed file. Accepts relative or absolute paths.
- `--annotate`: Path and filename of the TSV file containing the annotated mutations (3). Accepts relative or absolute paths.
- `--syn_muts`: "yes" or "no", depending on whether to include synonymous mutations in the output 
plot (default is "no") (4)

### 📝 Notes

  1. The reference genome must contain the coding sequences (CDS) of the proteins. If there are multiple proteins, they should be separated in the FASTA file. A pre-built AD169 CMV CDS reference is included at `minMutFinder-nf/resources/cmv_cds.fasta`.
  2. FASTA headers must use underscores (`_`) between words, or pipes (`|`) if
     the reference is a multi-CDS GenBank FASTA. For example:
     `>NC_006273_2_UL96` or `>UL97|NC_006273.2`.
  3. The annotated mutation file should be tab-separated and contain a column named `mutation` for annotated mutations.
  4. Use "yes" to include synonymous mutations in the output plot, or "no" to exclude them (default).

### Multi-sample

Run multiple samples at once with `run_multi.sh`, which reuses the tested
single-sample pipeline in parallel (concurrency-limited) and merges per-sample
results into master tables:

```bash
./run_multi.sh \
  --sample_dir <dir_with_one_subdir_per_sample> \
  --ref_seq <ref.fa> \
  --annotate <mutations.tsv> \
  --syn_muts yes \
  --concurrency 2
```

`run_multi.sh` accepts `--pipeline <nf_file>`, `--config <config_file>`, and
`--output_dir <dir>` in addition to the standard arguments.

Each sample subdirectory should contain a paired FASTQ set
(`{sample}_R1.fq.gz` + `{sample}_R2.fq.gz`, or `r1.fq`/`r2.fq`).
Outputs are merged into `all_mutations.csv` and `all_qc_metrics.csv` under
`--output_dir`.

---

## 🔏 License

This project, **minMutFinder**, is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this software under the terms of this license. For more details, refer to the [LICENSE](./LICENSE) file.

---

## 🖊️ Citing minMutFinder

A research paper on **minMutFinder** is currently in progress. In the meantime, please cite this GitHub repository as follows:

- **minMutFinder** (fork): https://github.com/murphygargantua-agent/minMutFinder
  (original by Prats-Méndez I., Vall d'Hebron Hospital Universitari & VHIR).

---

## 🔮 Future Work and Limitations

### Current Thresholds:
- **Allele Frequency (AF) ≥ 5%**
- **Read depth per nucleotide position ≥ 20**

### Future Improvements:
- Support for user-provided VCF files along with SAM/BAM files, skipping the initial quality control, mapping, and variant calling steps.
- Currently, **minMutFinder** has been tested with viral datasets. Future developments will expand its capabilities to handle bacterial datasets and larger organisms. There are also plans to explore its potential in oncologic research for cancer datasets.

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

"""
}

workflow {
    if (params.version || params.v) {
        versionMessage()
        exit 0
    }

    if (params.help || params.h || params.isEmpty()) {
        helpMessage()
        exit 0
    }

    //Creates working dir
    def workingpath = params.out_path
    def workingdir = file(workingpath)
    if( !workingdir.exists() ) {
        if( !workingdir.mkdirs() ) {
            exit 1, "Cannot create working directory: $workingpath"
        }
    }

    // Header log info
    log.info """---------------------------------------------
MINORITY MUTATION FINDER (minMutFinder)
---------------------------------------------

Beginning of analysis:

"""

    def summary = [:]

    summary['Starting time'] = new java.util.Date()
    //Environment
    summary['Environment'] = ""
    summary['Pipeline Name'] = 'minMutFinder'
    summary['Pipeline Version'] = workflow.manifest.version

    println "${params}"

    // Resolve relative input paths to absolute canonical paths so that
    // scripts running in isolated work directories can always find the files.
    def out_path = new File(params.out_path).getCanonicalPath()
    def ref_seq = new File(params.ref_seq).getCanonicalPath()
    def annotate = params.annotate ? new File(params.annotate).getCanonicalPath() : null
    def r1 = params.r1 && !params.r1.equals(false) ? new File(params.r1).getCanonicalPath() : null
    def r2 = params.r2 && !params.r2.equals(false) ? new File(params.r2).getCanonicalPath() : null

    // dirCreator returns 8 paths: out_path, fastq, assembly, references, variant_calling, mutations, plots, qc
    def out_paths = dirCreator(out_path)
    def (p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc) = out_paths
    refCheck(ref_seq)

    // fastqProcessing outputs: ref_seq.fasta, ref_seq.fasta.fai, qc_R1, qc_R2, qc_R1_unpaired, qc_R2_unpaired
    def (rfasta, rfasta_fai, qR1, qR2, qR1u, qR2u) = fastqProcessing(
                r1, r2, ref_seq,
                p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc)

    reads_qc_plot(
        r1, r2,
        p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc,
        ref_seq, rfasta_fai, qR1, qR2, qR1u, qR2u)

    def (map_sam, map_bam) = mapping1(
        p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc,
        ref_seq, rfasta_fai, qR1, qR2, qR1u, qR2u)

    def (v_bam, v_bam_bai, v_vcf, v_sam, v_samfile,
         v_vcf_gz, v_vcf_gz_csi, v_af05_vcf, v_af05_vcf_gz, v_af05_vcf_gz_csi, v_af05_fasta) = variant_calling(
        p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc, map_sam, map_bam)

    def (cons_sam, cons_bam) = mapping2(
        p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc,
        v_bam, v_bam_bai, v_vcf, v_sam, v_samfile, v_vcf_gz, v_vcf_gz_csi,
        v_af05_vcf, v_af05_vcf_gz, v_af05_vcf_gz_csi, v_af05_fasta)

    def depth_tsv = read_depth(
        p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc, map_sam, map_bam)

    def prot_names = protNames(
        p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc,
        v_bam, v_bam_bai, v_vcf, v_sam, v_samfile, v_vcf_gz, v_vcf_gz_csi,
        v_af05_vcf, v_af05_vcf_gz, v_af05_vcf_gz_csi, v_af05_fasta)

    def (mut_csv, per_qc, per_depth) = byProteinAnalysis(
        p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc,
        ref_seq, prot_names, depth_tsv)

    def final_muts = waitForOutMuts(
        p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc,
        v_bam, v_bam_bai, v_vcf, v_sam, v_samfile, v_vcf_gz, v_vcf_gz_csi,
        v_af05_vcf, v_af05_vcf_gz, v_af05_vcf_gz_csi, v_af05_fasta,
        mut_csv, per_qc, per_depth)

    if ( annotate ){
        annotateAndViz(p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc,
                       ref_seq, annotate, final_muts, params.syn_muts)
    } else {
        vizNoAnnot(p_out, p_fq, p_asm, p_ref, p_vcf, p_mut, p_plt, p_qc,
                   ref_seq, final_muts, params.syn_muts)
    }
}
