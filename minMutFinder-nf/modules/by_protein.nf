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

process prot_names {
    input:
    path out_dir
    file indelqual_prot_bam
    file indelqual_prot_bam_indexed
    file prot_variants_vcf
    file indelqual_prot_sam
    file samfile
    file prot_variants_vcf_gz
    file prot_variants_vcf_gz_indexed
    file vcf_file
    file prot_variants_AF05_vcf_gz
    file prot_variants_AF05_vcf_gz_indexed
    file prot_variants_AF05_fasta


    output:
    file "${out_dir}/assembly/*_prot_name"

    script:
    """
    DIR_SAMPLE=${out_dir}
    SAMPLE=${out_dir.baseName}
    ASSEMBLY=${out_dir}/assembly
    grep -v "#" ${vcf_file} | cut -d\$'\t' -f1 | uniq > ${out_dir}/assembly/prot_names.txt
    grep -v "#" ${vcf_file} | cut -d\$'\t' -f1 | uniq | awk -F "|" '{close(F); ID=\$0; gsub("^>", "", ID); F="'\${ASSEMBLY}/'"ID"_prot_name"} {print >> F}'
    """
}


process protein_analysis {
    conda 'conda-forge::python conda-forge::pandas conda-forge::numpy conda-forge::plotly conda-forge::biopython'

    input:
    path out_dir
    path ref_seq
    file protein
    file depth

    output:
    path("${out_dir}/mutations/${out_dir.baseName}_*_mutations.csv"), optional: true
    // file "${out_dir}/mutations/${out_dir.baseName}_*_mutations.csv"
    // file "${out_dir}/variant_calling/sam_to_loop_*.sam"
    // file "${out_dir}/variant_calling/samfile_*.tsv"
    file "${out_dir}/qc/*_qc_metrics.csv"
    file "${out_dir}/assembly/*_depth.tsv"

    script:
    """
    MUTATIONS=${out_dir}/mutations
    SAMPLE=${out_dir.baseName}
    while IFS= read -r prot; do
        echo -e ""\$prot" is the protein fasta name"
        python $projectDir/bin/by_protein_analysis.py --out-dir ${out_dir} --ref-seq ${ref_seq} --sample ${out_dir.baseName} --prot "\$prot"
        echo \$prot >> ${out_dir}/mutations/stopper
    done < ${out_dir}/assembly/prot_names.txt
    """
}
