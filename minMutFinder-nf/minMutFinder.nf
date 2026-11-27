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
include { ref_faidx; trimmomatic_qc; mapping1; lofreq_variant_calling; bcftools_consensus; read_depth; mapping2; reads_qc_plot } from './modules/from_fq_to_vcf'
include { byProteinAnalysis; protNames } from './modules/by_protein'
include { waitForOutMuts; vizNoAnnot; annotateAndViz } from './modules/viz_and_annotate'
include { runMultiQC } from './modules/multiqc'

def versionMessage() 
{
	log.info"""
	
	MINORITY MUTATION FINDER (minMutFinder) - Version: ${workflow.manifest.version} 
	""".stripIndent()
}

// Check that all required Python helper scripts are present. The original minMutFinder
// ships these in `bin/`; if they are missing the pipeline would fail with a cryptic
// "No such file or directory" error from inside a process, which gives no clue about
// what is actually wrong. This makes the failure explicit and tells the user where to
// get the files.
def checkScripts() {
    def scripts = [
        'bin/ref_check.py',
        'bin/read_metrics.py',
        'bin/by_protein_analysis.py',
        'bin/viz_no_annotate.py',
        'bin/annotate_and_viz.py'
    ]
    def missing = scripts.findAll { !file("$projectDir/$it").exists() }
    return missing
}

def validateInputs() {
    def missingScripts = checkScripts()
    if (!missingScripts.isEmpty()) {
        exit 1, """
        minMutFinder is missing required Python helper scripts: ${missingScripts.join(', ')}

        These files are not shipped with this repository. They come from the original
        minMutFinder (Prats-Méndez I., VHIR) and must be placed in the `bin/` directory
        of this project:

            $projectDir/bin/ref_check.py
            $projectDir/bin/read_metrics.py
            $projectDir/bin/by_protein_analysis.py
            $projectDir/bin/viz_no_annotate.py
            $projectDir/bin/annotate_and_viz.py

        Until they are present, the pipeline cannot run.
        """.stripIndent()
    }
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

Run multiple samples at once. Two mutually exclusive input modes:

**1. Manifest mode** — a CSV with `sample_name`, `r1` and `r2`:

```bash
nextflow run minMutFinder.nf \
  --manifest samples.csv \
  --ref_seq <ref.fa> \
  --annotate <mutations.tsv> \
  --syn_muts yes \
  --out_path <output_dir>
```

The manifest CSV (`samples.csv`) must have a header row and the columns:
`sample_name`, `r1`, `r2`. Paths may be relative (resolved against the CWD) or absolute.

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
nextflow run minMutFinder.nf \
  --inputReads /path/to/samples/ \
  --ref_seq <ref.fa> \
  --out_path <output_dir>
```

Expect one subdirectory per sample, each containing `{sample}_R1.fq.gz` and
`{sample}_R2.fq.gz`. The naming convention is the same as `run_multi.sh`'s
`find_pair` step:

- `*_R1.*` + `*_R2.*` (any extension: `.fastq`, `.fq`, `.fq.gz`)
- `r1*`/`r2*` or `R1*`/`R2*` as a fallback
- files may sit directly in `--inputReads` or one level down inside per-sample dirs

Outputs are merged into `all_mutations.csv` and `all_qc_metrics.csv` under `--out_path`:

```bash
Per-sample outputs land in `<out_path>/<sample_name>/`; merged tables
`all_mutations.csv` and `all_qc_metrics.csv` are written to `<out_path>`.

A consolidated MultiQC QC report is generated under `<out_path>/multiqc_report/`:
it merges the per-sample `qc/qc_metrics.csv` tables, `fastp` JSON logs and
`bcftools stats` reports into one HTML report.

To smoke-test the orchestration without installing bioinformatics tools, swap in the stub pipeline:

```bash
nextflow run minMutFinder.nf --manifest samples.csv --ref_seq minMutFinder-nf/resources/cmv_cds.fasta --pipeline minMutFinder-nf/minMutFinder_stub.nf
```

**Note**: The legacy single-sample mode (`--r1`/`--r2`) is still supported; each sample still
runs through the same per-sample pipeline, just driven by Nextflow instead of `run_multi.sh`.


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

process mergeMutations {
    input:
    path muts_csv
    path out_root

    output:
    path "${out_root}/all_mutations.csv"

    script:
    """
    echo "SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency" > ${out_root}/all_mutations.csv
    for f in ${muts_csv}; do
        if [[ -f "\$f" ]]; then
            tail -n +2 "\$f" >> ${out_root}/all_mutations.csv
        fi
    done
    """
}

process mergeQc {
    input:
    path qc_csv
    path out_root

    output:
    path "${out_root}/all_qc_metrics.csv"

    script:
    """
    echo "sample,test,score" > ${out_root}/all_qc_metrics.csv
    for f in ${qc_csv}; do
        if [[ -f "\$f" ]]; then
            tail -n +2 "\$f" >> ${out_root}/all_qc_metrics.csv
        fi
    done
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

    // Fail fast with a helpful message when the required Python helper scripts are
    // missing, instead of dying later with an opaque "No such file or directory"
    // from inside a process script.
    validateInputs()

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

    // ---------------------------------------------------------------------------
    // Sample input: two mutually exclusive modes.
    //   1. Manifest mode  --manifest samples.csv  → CSV with explicit r1/r2 paths
    //   2. Scan mode    --inputReads <dir>        → auto-pair *_R1.* / *_R2.* FASTQ
    // Falls back to the legacy single-sample --r1/--r2 when neither is given.
    // ---------------------------------------------------------------------------
    if (params.manifest) {
        // Manifest mode: read the CSV; each row becomes a [sample_name, r1, r2] tuple of
        // canonical paths so scripts in isolated work dirs can always find the files.
        def manifestChannel = Channel.fromPath(params.manifest).splitCsv(header: true, quote: '"')
            .map { row -> [row.sample_name, new File(row.r1).getCanonicalPath(), new File(row.r2).getCanonicalPath()] }
        manifestChannel.set { ch_samples }
    } else if (params.inputReads) {
        // Scan mode: auto-pair forward/reverse FASTQ files by name.
        // The default pairing convention mirrors run_multi.sh's find_pair():
        //   *_R1.* + *_R2.* (any extension: .fastq, .fq, .fq.gz)
        // Custom ERE patterns can be provided with --r1_pattern / --r2_pattern,
        // e.g. '^.+_R1_001\.fastq$' for files with an underscore suffix.
        def inputReads = new File(params.inputReads).getCanonicalPath()
        def r1_pattern = params.r1_pattern ? params.r1_pattern : '.*_R1\\..*'
        def r2_pattern = params.r2_pattern ? params.r2_pattern : '.*_R2\\..*'

        // Brace expansion {*,*/*} finds files directly in --inputReads or one level
        // down inside per-sample dirs (e.g. samples/S1/synth_R1.fq).
        // Sample id: filename prefix before the R1/R2 marker (flat layout), or the
        // parent directory name (one level down inside a sample dir).
        // Note: the suffix must be stripped regardless of what follows it (e.g.
        // SampleA_R1_001.fastq -> SampleA, not SampleA_R1_001.fastq), otherwise
        // the R1/R2 join never matches and the channel ends up empty.
        def r1_files = Channel.fromPath("$inputReads/{*,*/*}")
            .filter { file -> file.getName().matches(r1_pattern) }
            .map { file -> [
                file.getParent().toString() == inputReads ?
                    file.getName().replaceAll(/(_R1|_R2|_R1_|_R2_).*/, '') : // flat: prefix as sample id
                    file.getParent().getName(),                            // nested: dir name
                file.toString()
              ]
            }
        def r2_files = Channel.fromPath("$inputReads/{*,*/*}")
            .filter { file -> file.getName().matches(r2_pattern) }
            .map { file -> [
                file.getParent().toString() == inputReads ?
                    file.getName().replaceAll(/(_R1|_R2|_R1_|_R2_).*/, '') : // flat: prefix as sample id
                    file.getParent().getName(),                            // nested: dir name
                file.toString()
              ]
            }

        // Join R1/R2 files by sample id.
        r1_files.join(r2_files).set { ch_samples }
    } else {
        // Legacy single-sample mode: out_dir is the user-supplied out_path so the
        // per-sample naming behaviour (out_path.baseName) is preserved.
        def r1 = params.r1 && !params.r1.equals(false) ? new File(params.r1).getCanonicalPath() : null
        def r2 = params.r2 && !params.r2.equals(false) ? new File(params.r2).getCanonicalPath() : null
        Channel.of([null, r1, r2]).set { ch_samples }
    }

    // Attach the per-sample output dir (out_path/<sample_name>) and the resolved reference
    // to every sample. Nextflow zips the fields of the channel element-wise so each
    // instance of the chain below gets one complete sample.
    ch_samples.map { name, r1, r2 ->
        // All input modes must use a per-sample dir under out_path: if all samples
        // share one directory (scan/legacy mode), their outputs collide on the same
        // file names and overwrite each other. Manifest mode already did this; extend
        // it to every mode so multi-sample runs produce one subdir per sample.
        def out_dir = new File(out_path, name).getCanonicalPath()
        // Pre-create the per-sample dir: Nextflow symlinks path inputs into the work dir,
        // and process scripts can only mkdir through a symlink whose target exists.
        new File(out_dir).mkdirs()
        [name, r1, r2, ref_seq, out_dir]
    }.set { ch }

    // Split the sample tuple into individual channels; Nextflow pairs them element-wise
    // (the `ref_seq` field is broadcast to every instance).
    def ch_names = ch.map { it[0] }
    def ch_r1 = ch.map { it[1] }
    def ch_r2 = ch.map { it[2] }
    def ch_ref_seq = ch.map { it[3] }
    // The per-sample output dirs are reused by the MultiQC step at the end of the
    // workflow, so assign them to a value channel that can be subscribed to more
    // than once by the processes below.
    def ch_out_dir = ch.map { it[4] }

    // Create the per-sample output directories.
    dirCreator(ch_out_dir)
    refCheck(ref_seq)

    // ---------------------------------------------------------------------------
    // Per-sample pipeline chain. Nextflow fans out each process once per sample
    // (one DAG instance per channel element), then joins the channels again.
    // ---------------------------------------------------------------------------
    // Reference indexing (samtools env) and trimmomatic QC (python/trimmomatic env) are
    // now separate processes: a failure in one env never breaks the other.
    def (rfasta, rfasta_fai) = ref_faidx(ch_ref_seq, ch_out_dir)

    def (qR1, qR2, qR1u, qR2u) = trimmomatic_qc(ch_r1, ch_r2, ch_out_dir)

    reads_qc_plot(ch_r1, ch_r2, ch_out_dir, ch_ref_seq, rfasta_fai, qR1, qR2, qR1u, qR2u)

    def (map_sam, map_bam) = mapping1(ch_r1, ch_r2, ch_ref_seq, ch_out_dir, rfasta_fai, qR1, qR2, qR1u, qR2u)

    // lofreq does variant calling (bam/sam/vcf); bcftools does consensus/AF filtering.
    def (v_bam, v_bam_bai, v_vcf, v_sam, v_samfile) = lofreq_variant_calling(
        ch_out_dir, map_sam, map_bam)

    def (v_vcf_gz, v_vcf_gz_csi, v_af05_vcf, v_af05_vcf_gz, v_af05_vcf_gz_csi, v_af05_fasta) = bcftools_consensus(
        ch_out_dir, v_vcf)

    def (cons_sam, cons_bam) = mapping2(
        ch_out_dir, v_bam, v_bam_bai, v_vcf, v_sam, v_samfile, v_vcf_gz, v_vcf_gz_csi,
        v_af05_vcf, v_af05_vcf_gz, v_af05_vcf_gz_csi, v_af05_fasta)

    def depth_tsv = read_depth(ch_out_dir, map_sam, map_bam)

    def prot_names = protNames(
        ch_out_dir, v_bam, v_bam_bai, v_vcf, v_sam, v_samfile, v_vcf_gz, v_vcf_gz_csi,
        v_af05_vcf, v_af05_vcf_gz, v_af05_vcf_gz_csi, v_af05_fasta)

    def (mut_csv, per_qc, per_depth) = byProteinAnalysis(
        ch_out_dir, ch_ref_seq, prot_names, depth_tsv)

    def final_muts = waitForOutMuts(
        ch_out_dir, v_bam, v_bam_bai, v_vcf, v_sam, v_samfile, v_vcf_gz, v_vcf_gz_csi,
        v_af05_vcf, v_af05_vcf_gz, v_af05_vcf_gz_csi, v_af05_fasta,
        mut_csv, per_qc, per_depth)

    if ( annotate ){
        annotateAndViz(ch_out_dir, ch_ref_seq, final_muts, annotate, params.syn_muts)
    } else {
        vizNoAnnot(ch_out_dir, ch_ref_seq, final_muts, params.syn_muts)
    }

    // ---------------------------------------------------------------------------
    // Merge per-sample results into master tables under out_path.
    // The `.collect()` call reduces each channel to a single list of paths so the
    // merge process is executed once with every sample (rather than once per
    // sample), which is what the bash orchestrator did before.
    // ---------------------------------------------------------------------------
    def all_muts = final_muts.collect()
    def all_qc = per_qc.collect()
    // Collect from `per_qc` (produced by byProteinAnalysis) so runMultiQC waits
    // for every per-sample process to finish before it starts. `per_qc` emits a
    // list per CDS; some CDSs produce no output, so filter nulls and map each QC
    // metrics file path back to its per-sample dir, then dedupe.
    def all_dirs = per_qc.filter { it != null }.map {
        def f = it instanceof List ? it[0] : it
        f == null ? null : f.getParent().getParent()
    }.filter { it != null }.unique().collect()

    if (all_muts) { mergeMutations(all_muts, out_path) }
    if (all_qc) { mergeQc(all_qc, out_path) }
    // Aggregate the per-sample QC logs and metrics into a single MultiQC report.
    if (all_dirs) { runMultiQC(all_dirs, out_path) }
}
