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

process runMultiQC {
    tag "${sample_dirs}"

    conda 'bioconda::multiqc=1.35 bioconda::fastp=0.24.0 bioconda::samtools=1.13 bioconda::bcftools=1.15'

    input:
    path sample_dirs
    path out_root

    output:
    path "${out_root}/multiqc_report"

    script:
    """
    mkdir -p ${out_root}/multiqc_report

    # Resolve every sample dir to an absolute path so the script keeps working
    # no matter which directory Nextflow sets as its working directory.
    # (bash arrays are not supported in Nextflow task scripts, so use a
    # space-separated string instead).
    abs_dirs=""
    for dir in ${sample_dirs}; do
        if [[ -d \$dir ]]; then
            abs_dirs="\$abs_dirs \$(cd \$dir && pwd -P)"
        else
            abs_dirs="\$abs_dirs \$dir"
        fi
    done

    # Generate QC logs that MultiQC's built-in modules consume, plus the
    # per-sample custom-content CSVs. Every command runs unconditionally and
    # errors are caught with || so a missing input for one sample does not
    # abort the whole report (bash globs can be unreliable inside Nextflow
    # task scripts, so no `if` guards here).
    for dir in \$abs_dirs; do
        sample=\$(basename "\$dir")

        # fastp on the trimmed QC reads -> JSON log for the fastp module
        if ls "\$dir/fastq/"*_qc_R1.fastq.gz >/dev/null 2>&1; then
            fastp -i "\$dir/fastq/"*_qc_R1.fastq.gz \
                -I "\$dir/fastq/"*_qc_R2.fastq.gz \
                -o /tmp/mqc_out1.fastq.gz -O /tmp/mqc_out2.fastq.gz \
                --json "\$dir/fastq/\${sample}_fastp.json" \
                --html "\$dir/fastq/\${sample}_fastp.html" \
                2>&1 || echo "WARN: fastp failed for \$sample"
        else
            echo "WARN: no QC reads found for \$sample"
        fi

        # samtools stats on the mapped BAM -> log for the samtools module
        if [[ -f "\$dir/assembly/\${sample}_reads-mapped_to_prot.bam" ]]; then
            samtools stats "\$dir/assembly/\${sample}_reads-mapped_to_prot.bam" \
                > "\$dir/assembly/\${sample}_samtools_stats.txt" \
                2>&1 || echo "WARN: samtools stats failed for \$sample"
        else
            echo "WARN: no mapped BAM found for \$sample"
        fi

        # bcftools stats on the variant VCF -> log for the bcftools module
        if [[ -f "\$dir/variant_calling/\${sample}_prot_variants.vcf" ]]; then
            bcftools stats "\$dir/variant_calling/\${sample}_prot_variants.vcf" \
                > "\$dir/variant_calling/\${sample}_bcftools_stats.txt" \
                2>&1 || echo "WARN: bcftools stats failed for \$sample"
        else
            echo "WARN: no VCF found for \$sample"
        fi

        # Convert per-sample QC metrics CSV into MultiQC custom-content format
        if [[ -f "\$dir/qc/qc_metrics.csv" ]]; then
            python3 ${projectDir}/bin/qc_metrics_mqc.py \
                --in "\$dir/qc/qc_metrics.csv" \
                --out "\$dir/qc/\${sample}_mqc.csv" \
                --sample "\$sample" \
                2>&1 || echo "WARN: qc_metrics_mqc failed for \$sample"
        else
            echo "WARN: no qc_metrics.csv found for \$sample"
        fi
    done

    # Run MultiQC, scanning all per-sample output dirs and merging the results
    # into a single report under ${out_root}/multiqc_report.
    multiqc \$abs_dirs --outdir ${out_root}/multiqc_report --force --config ${projectDir}/multiqc_config.yaml
    """
}
