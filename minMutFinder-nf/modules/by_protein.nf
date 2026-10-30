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

process protNames {
    input:
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path
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
    file "${assembly_path}/*_prot_name"
    
    script:
    """
    DIR_SAMPLE=${out_path}
    SAMPLE=${out_path.baseName}
    ASSEMBLY=${assembly_path}
    grep -v "#" ${vcf_file} | cut -d\$'\\t' -f1 | uniq > ${assembly_path}/prot_names.txt
    grep -v "#" ${vcf_file} | cut -d\$'\\t' -f1 | uniq | awk -v ASSEMBLY="${assembly_path}" -F "|" '{close(F); ID=$0; gsub("^>", "", ID); F=ASSEMBLY"/"ID"_prot_name"; print >> F}'
    """
}


process byProteinAnalysis {
    input:    
    path out_path
    path fq_path
    path assembly_path
    path ref_path
    path vcf_path
    path muts_path
    path plots_path
    path qc_path
    path ref_seq
    file protein
    file depth

    output:
    path("${muts_path}/${out_path.baseName}_*_mutations.csv"), optional: true
    // file "${muts_path}/${out_path.baseName}_*_mutations.csv"
    // file "${vcf_path}/sam_to_loop_*.sam"
    // file "${vcf_path}/samfile_*.tsv"
    file "${qc_path}/*_qc_metrics.csv"
    file "${assembly_path}/*_depth.tsv"


    script:
    """
    MUTATIONS=${muts_path}
    awk -v pd="$projectDir" -v op="$out_path" -v rs="$ref_seq" -v sm="${out_path.baseName}" -v mp="$muts_path" '
    {
        cmd = "python " pd "/bin/by_protein_analysis.py --out-dir " op " --ref-seq " rs " --sample " sm " --prot \"" $0 "\"; echo \"" $0 "\" >> " mp "/stopper"
        system(cmd)
    }' ${assembly_path}/prot_names.txt
    """
}