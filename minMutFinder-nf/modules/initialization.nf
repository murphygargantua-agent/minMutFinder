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

process create_dirs {

    input:
    path out_dir

    output:
    path "${out_dir}"
    path "${out_dir}/fastq"
    path "${out_dir}/assembly"
    path "${out_dir}/assembly/references"
    path "${out_dir}/variant_calling"
    path "${out_dir}/mutations"
    path "${out_dir}/plots"
    path "${out_dir}/qc"

    script:
    """
    #!/bin/bash
    echo Creating directories..
    DIR_SAMPLE="${out_dir}"
    if [[ -d "\${DIR_SAMPLE}" ]]; then echo -e ""\${DIR_SAMPLE}" output folder already exists"; else mkdir -p "\${DIR_SAMPLE}"; fi
    QC_DIR="${out_dir}/qc"
    if [[ -d "\${QC_DIR}" ]]; then echo -e ""\${QC_DIR}" output folder already exists"; else mkdir "\${QC_DIR}"; fi
    FASTQ="${out_dir}/fastq"
    if [[ -d "\${FASTQ}" ]]; then echo -e ""\${FASTQ}" output folder already exists"; else mkdir "\${FASTQ}"; fi
    ASSEMBLY="${out_dir}/assembly"
    REF_PATH="${out_dir}/assembly/references"
    if [[ -d "\${ASSEMBLY}" ]]; then echo -e ""\${ASSEMBLY}" output folder already exists"; else mkdir "\${ASSEMBLY}"; fi
    if [[ -d "\${REF_PATH}" ]]; then echo -e ""\${REF_PATH}" output folder already exists"; else mkdir "\${REF_PATH}"; fi
    VARIANT_CALLING="${out_dir}/variant_calling"
    if [[ -d "\$VARIANT_CALLING" ]]; then echo -e ""\${VARIANT_CALLING}" output folder already exists"; else mkdir "\${VARIANT_CALLING}"; fi
    MUTATIONS="${out_dir}/mutations"
    if [[ -d "\${MUTATIONS}" ]]; then echo -e ""\${MUTATIONS}" output folder already exists"; else mkdir \${MUTATIONS}; fi
    touch ${out_dir}/mutations/stopper
    PLOTS="${out_dir}/plots"
    if [[ -d "\${PLOTS}" ]]; then echo -e ""\${PLOTS}" output folder already exists"; else mkdir "\${PLOTS}"; fi
    """
}


process check_ref {
    conda 'conda-forge::python'

    input:
    val ref_seq

    output:
    stdout

    script:
    """
    python3 $projectDir/bin/ref_check.py --ref-seq "${ref_seq}"
    """
}
