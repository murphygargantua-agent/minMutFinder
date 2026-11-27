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

process wait_for_out_muts {
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
    file mutations
    file qc_metrics
    file depth

    output:
    path("${out_dir}/mutations/${out_dir.baseName}_mutations.csv"), optional: true

    script:
    """
#    if [[ \$(wc -l ${out_dir}/assembly/prot_names.txt | cut -d\$' ' -f1) == \$(wc -l ${out_dir}/mutations/stopper | cut -d\$' ' -f1) ]]; then
    if [[ \$(awk 'END {print NR}' ${out_dir}/assembly/prot_names.txt) == \$(awk 'END {print NR}' ${out_dir}/mutations/stopper) ]]; then
        rm ${out_dir}/mutations/stopper
        echo "SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency" > ${out_dir}/mutations/${out_dir.baseName}_mutations.csv
        for file in ${out_dir}/mutations/"${out_dir.baseName}"_*_mutations.csv; do
            if [[ -f \$file ]]; then
                tail -n +2 \$file >> ${out_dir}/mutations/${out_dir.baseName}_mutations.csv
                rm \$file
            fi
        done
#        sed -i 1i"SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency" ${out_dir}/mutations/${out_dir.baseName}_mutations.csv
#        sed -i '' '1i\\SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency' ${out_dir}/mutations/${out_dir.baseName}_mutations.csv
    fi
    """
}

process visualize_no_annotate {
    conda 'conda-forge::python conda-forge::pandas conda-forge::numpy conda-forge::plotly conda-forge::biopython'

    input:
    path out_dir
    path ref_seq
    file muts
    val syn_muts

    output:
    stdout

    script:
    """
    python3 $projectDir/bin/viz_no_annotate.py --out-dir ${out_dir} --ref-seq ${ref_seq} --prot-names ${out_dir}/assembly/prot_names.txt --syn-muts ${params.syn_muts}
    """
}

process annotateAndViz {
    conda 'conda-forge::python conda-forge::pandas conda-forge::numpy conda-forge::plotly conda-forge::biopython'

    input:
    path out_dir
    path ref_seq
    file annotate
    file muts
    val syn_muts

    output:
    stdout

    script:
    """
    python3 $projectDir/bin/annotate_and_viz.py --out-dir ${out_dir} --annotate ${annotate} --ref-seq ${ref_seq} --prot-names ${out_dir}/assembly/prot_names.txt --syn-muts ${params.syn_muts}
    """
}
