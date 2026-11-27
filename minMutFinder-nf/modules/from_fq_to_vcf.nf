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
// Copyright (C) 2024 Respiratory Viruses Unit, Microbiology Department, Vall d'Hebron Hospital Universitari, Vall d'Hebron Institut de Recerca (VHIR), Vall d'Hebron Barcelona Hospital Campus, Passeig Vall d'Hebron 119-129, 08035 Barcelona, Spain
//
// Every bioinformatics tool lives in its own process with its own conda directive.
// Mixing tools in one env means a single broken dependency (e.g. an old Python
// breaking the trimmomatic wrapper) poisons every step that shares that env.

process ref_faidx {
    // samtools 1.25 -- latest stable (bioconda).
    conda 'bioconda::samtools=1.25'

    input:
    path ref_seq
    path out_dir

    output:
    file "${out_dir}/assembly/ref_seq.fasta"
    file "${out_dir}/assembly/ref_seq.fasta.fai"

    script:
    """
    DIR_SAMPLE=${out_dir}
    ASSEMBLY=${out_dir}/assembly
    mkdir -p "\${ASSEMBLY}"
    cp ${ref_seq} "\${ASSEMBLY}/ref_seq.fasta"
    samtools faidx "\${ASSEMBLY}/ref_seq.fasta"
    """
}

process trimmomatic_qc {
    // trimmomatic 0.41 (latest, released Aug 2025) -- the conda wrapper is a Python
    // script using f-strings, so python=3.11 must be pinned to avoid "SyntaxError"
    // before Java ever runs.
    conda 'bioconda::trimmomatic=0.41 conda-forge::python=3.11'

    input:
    path r1
    path r2
    path out_dir

    output:
    file "${out_dir}/fastq/${out_dir.baseName}_qc_R1.fastq.gz"
    file "${out_dir}/fastq/${out_dir.baseName}_qc_R2.fastq.gz"
    file "${out_dir}/fastq/${out_dir.baseName}_qc_R1_unpaired.fastq.gz"
    file "${out_dir}/fastq/${out_dir.baseName}_qc_R2_unpaired.fastq.gz"

    script:
    """
    DIR_SAMPLE=${out_dir}
    FASTQ=${out_dir}/fastq
    mkdir -p "\${FASTQ}"
    SAMPLE=\$(basename \${DIR_SAMPLE})
    echo "\${SAMPLE}"

    R1=${r1}
    R2=${r2}
    FQ_R1="\${FASTQ}"/"\${SAMPLE}"_qc_R1.fastq.gz; FQ_R2="\${FASTQ}"/"\${SAMPLE}"_qc_R2.fastq.gz
    printf "\nTRIMMOMATIC: \nSoft filtering fastq files from request ID: \${SAMPLE}\n\n"
    time (trimmomatic PE -threads 64 -phred33 "\$R1" "\$R2" "\${FQ_R1}" "\${FASTQ}"/"\${SAMPLE}"_qc_R1_unpaired.fastq.gz "\${FQ_R2}" "\${FASTQ}"/"\${SAMPLE}"_qc_R2_unpaired.fastq.gz LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30)
    """
}

process mapping1 {
    // minimap2 2.31 (latest) + samtools 1.25 -- used together for the initial map
    // and the sort step, but kept in one env because they're both C tools with
    // compatible dependencies.
    conda 'bioconda::minimap2=2.31 bioconda::samtools=1.25'

    input:
    path r1
    path r2
    path ref_seq
    path out_dir
    file ref_seq_fai
    file fqr1
    file fqr2
    file fqUPr1
    file fqUPr2

    output:
    file "${out_dir}/assembly/${out_dir.baseName}_reads-mapped_to_prot.sam"
    file "${out_dir}/assembly/${out_dir.baseName}_reads-mapped_to_prot.bam"

    script:
    """
    DIR_SAMPLE=${out_dir}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${out_dir}/assembly
    FQ_R1=${fqr1}
    FQ_R2=${fqr2}
    echo -e "\nMINIMAP2\n->Mapping against \${ASSEMBLY}/ref_seq.fasta\n"
    time (minimap2 -ax sr "\${ASSEMBLY}/ref_seq.fasta" \$FQ_R1 \$FQ_R2 > "\${ASSEMBLY}"/"\${SAMPLE}"_reads-mapped_to_prot.sam)
    samtools view -bS "\${ASSEMBLY}"/"\${SAMPLE}"_reads-mapped_to_prot.sam | samtools sort - -o "\${ASSEMBLY}"/"\${SAMPLE}"_reads-mapped_to_prot.bam
    """
}

process lofreq_variant_calling {
    // lofreq 2.1.5-16 (latest) + samtools 1.25 -- variant calling and the bam/sam
    // index/view steps. bcftools is a separate process so its env stays clean.
    conda 'bioconda::lofreq=2.1.5-16 bioconda::samtools=1.25'

    input:
    path out_dir
    file sam
    file bam

    output:
    file "${out_dir}/variant_calling/${out_dir.baseName}_indelqual_prot.bam"
    file "${out_dir}/variant_calling/${out_dir.baseName}_indelqual_prot.bam.bai"
    file "${out_dir}/variant_calling/${out_dir.baseName}_prot_variants.vcf"
    file "${out_dir}/variant_calling/${out_dir.baseName}_indelqual_prot.sam"
    file "${out_dir}/variant_calling/samfile.sam"

    script:
    """
    DIR_SAMPLE=${out_dir}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${out_dir}/assembly
    VARIANT_CALLING=${out_dir}/variant_calling
    BAM2=${bam}

    # Call variants
    # tomake lofreq work
    echo -e "\nLOFREQ\n->Variant Calling Step: Mutation detection\n"

    lofreq indelqual --dindel -f "\${ASSEMBLY}/ref_seq.fasta" -o "\${VARIANT_CALLING}"/"\${SAMPLE}"_indelqual_prot.bam "\$BAM2"
    INDELQUAL_BAM="\${VARIANT_CALLING}"/"\${SAMPLE}"_indelqual_prot.bam
    samtools index "\${INDELQUAL_BAM}"
    # tomake lofreq work
    time (lofreq call-parallel --pp-threads 64 --call-indels -f "\${ASSEMBLY}/ref_seq.fasta" -o "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants.vcf "\${INDELQUAL_BAM}")

    samtools index "\${INDELQUAL_BAM}"
    samtools view -h -o "\${VARIANT_CALLING}"/"\${SAMPLE}"_indelqual_prot.sam "\${INDELQUAL_BAM}"

    awk '(\$6 != "*") && (\$12 != "") && (\$3 != "*")' "\${VARIANT_CALLING}"/"\${SAMPLE}"_indelqual_prot.sam > "\${VARIANT_CALLING}"/samfile.sam
    """
}

process bcftools_consensus {
    // bcftools 1.25 (latest) -- consensus/AF filtering/bgzip/indexing. Kept apart
    // from lofreq so a bcftools env issue never breaks variant calling.
    conda 'bioconda::bcftools=1.25'

    input:
    path out_dir
    file prot_variants_vcf

    output:
    file "${out_dir}/variant_calling/${out_dir.baseName}_prot_variants.vcf.gz"
    file "${out_dir}/variant_calling/${out_dir.baseName}_prot_variants.vcf.gz.csi"
    file "${out_dir}/variant_calling/${out_dir.baseName}_prot_variants_AF-0.05.vcf"
    file "${out_dir}/variant_calling/${out_dir.baseName}_prot_variants_AF-0.05.vcf.gz"
    file "${out_dir}/variant_calling/${out_dir.baseName}_prot_variants_AF-0.05.vcf.gz.csi"
    file "${out_dir}/assembly/${out_dir.baseName}_prot_variants_AF-0.05.fasta"

    script:
    """
    DIR_SAMPLE=${out_dir}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${out_dir}/assembly
    VARIANT_CALLING=${out_dir}/variant_calling
    FREQ_CUTOFF=0.05

    echo -e "\nBCFTOOLS\n->Variant Calling Step: Consensus sequence obtantion\n"

    bgzip -c "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants.vcf > "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants.vcf.gz
    bcftools index "\${VARIANT_CALLING}"/"\${SAMPLE}"_prot_variants.vcf.gz
    FREQ_CUTOFF=0.05
    bcftools filter -i "AF>=\${FREQ_CUTOFF} & DP>=20" "\${VARIANT_CALLING}/\${SAMPLE}_prot_variants.vcf.gz" > "\${VARIANT_CALLING}/\${SAMPLE}_prot_variants_AF-\${FREQ_CUTOFF}.vcf"
    VCF_FILE_CUTOFF="\${VARIANT_CALLING}/\${SAMPLE}_prot_variants_AF-\${FREQ_CUTOFF}.vcf"
    bgzip -c "\${VARIANT_CALLING}/\${SAMPLE}_prot_variants_AF-\${FREQ_CUTOFF}.vcf" > "\${VARIANT_CALLING}/\${SAMPLE}_prot_variants_AF-\${FREQ_CUTOFF}.vcf.gz"
    bcftools index \${VARIANT_CALLING}/"\${SAMPLE}"_prot_variants_AF-"\${FREQ_CUTOFF}".vcf.gz
    cat "\${ASSEMBLY}/ref_seq.fasta" | bcftools consensus \${VARIANT_CALLING}/\${SAMPLE}_prot_variants_AF-\${FREQ_CUTOFF}.vcf.gz -o \${ASSEMBLY}/\${SAMPLE}_prot_variants_AF-\${FREQ_CUTOFF}.fasta
    """
}

process mapping2 {
    // minimap2 2.31 (latest) + samtools 1.25 -- mapping against the consensus.
    conda 'bioconda::minimap2=2.31 bioconda::samtools=1.25'

    input:
    path out_dir
    file indelqual_prot_bam
    file indelqual_prot_bam_indexed
    file prot_variants_vcf
    file indelqual_prot_sam
    file samfile
    file prot_variants_vcf_gz
    file prot_variants_vcf_gz_indexed
    file prot_variants_AF05_vcf
    file prot_variants_AF05_vcf_gz
    file prot_variants_AF05_vcf_gz_indexed
    file prot_variants_AF05_fasta

    output:
    file "${out_dir}/assembly/${out_dir.baseName}_reads_mapped_consensus.sam"
    file "${out_dir}/assembly/${out_dir.baseName}_reads_mapped_consensus.bam"

    script:
    """
    DIR_SAMPLE=${out_dir}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${out_dir}/assembly
    VARIANT_CALLING=${out_dir}/variant_calling
    FASTQ=${out_dir}/fastq
    FREQ_CUTOFF=0.05
    FQ_R1="\${FASTQ}"/"\${SAMPLE}"_qc_R1.fastq.gz; FQ_R2="\${FASTQ}"/"\${SAMPLE}"_qc_R2.fastq.gz
    time (minimap2 -ax sr "\${prot_variants_AF05_fasta}" \$FQ_R1 \$FQ_R2 > \${ASSEMBLY}/"\${SAMPLE}"_reads_mapped_consensus.sam)
    samtools view -bS \${ASSEMBLY}/"\${SAMPLE}"_reads_mapped_consensus.sam | samtools sort - -o \${ASSEMBLY}/"\${SAMPLE}"_reads_mapped_consensus.bam
    """


}


process read_depth {
    // samtools 1.25 -- depth calculation only.
    conda 'bioconda::samtools=1.25'

    input:
    path out_dir
    file sam
    file bam


    output:
    file "${out_dir}/assembly/${out_dir.baseName}_depth_consensus.tsv"


    script:
    """
    DIR_SAMPLE=${out_dir}
    SAMPLE=\$(basename \${DIR_SAMPLE})
    ASSEMBLY=${out_dir}/assembly
    VARIANT_CALLING=${out_dir}/variant_calling
    FREQ_CUTOFF=0.05
    BAM2=${bam}
    echo -e "ref\tpos\tdepth" > \${ASSEMBLY}/"\${SAMPLE}"_depth_consensus.tsv
    samtools depth -a \$BAM2 >> \${ASSEMBLY}/"\${SAMPLE}"_depth_consensus.tsv
    """
}

process reads_qc_plot {
    conda 'conda-forge::python=3.11 conda-forge::matplotlib=3.8 bioconda::pysam conda-forge::pandas conda-forge::seaborn conda-forge::biopython'

    input:
    file r1
    file r2
    path out_dir
    file ref_seq
    file ref_seq_fai
    file fqr1
    file fqr2
    file fqUPr1
    file fqUPr2

    output:
    file "${out_dir}/qc/qc_metrics.csv"
    file "${out_dir}/plots/qc_metrics.svg"


    script:
    """
    python3 $projectDir/bin/read_metrics.py --out-dir ${out_dir} --r1 ${r1} --r2 ${r2} --sample ${out_dir.baseName} --r1qc ${fqr1} --r1UPqc ${fqUPr1} --r2qc ${fqr2} --r2UPqc ${fqUPr2}
    """
}
