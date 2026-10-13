#!/usr/bin/env python3

# This file is part of minMutFinder.
#
# minMutFinder is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# minMutFinder is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with minMutFinder. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2024 Respiratory Viruses Unit, Microbiology Department, Vall d’Hebron Hospital Universitari, Vall d’Hebron Institut de Recerca (VHIR), Vall d’Hebron Barcelona Hospital Campus, Passeig Vall d’Hebron 119-129, 08035 Barcelona, Spain

import re
from cigar_remove import cigar_remove

def reads_assembler(samfile, all_fastas):
    print('Beginning reads assembly of ' + str(samfile))
    infile = open(samfile, 'r')
    while infile:
        line = infile.readline()
        if (line == ""):
            break
        # Create segments
        node, flag = line.split('\t')[0], line.split('\t')[1]
        if (int(flag) == -1):
            isit = 0
        else:
            isit = 1
            pos, pos_end = int(line.split('\t')[3]), int(line.split('\t')[3])
            cigar, sequence, p, idx, seq_init = line.split('\t')[5], line.split('\t')[9], re.compile("[a-zA-Z]"), 0, 0
            class_del, number_class, new_cig, startH, H_exists = 0, 0, '', 0, False

            if cigar == "*":
                continue

            cigar = cigar_remove(cigar, startH, new_cig, idx, number_class, H_exists, p)

            idx = 0
            segments = []  # (INI, END, CLASS, SEQ)
            for c in p.finditer(cigar):
                end, seq_end = int(pos) + int(cigar[idx:c.start()]) - 1, seq_init + int(cigar[idx:c.start()])
                if (number_class == 0) & (c.group() == "S"):
                    pos = int(pos) - int(cigar[idx:c.start()])

                if c.group() == "D":
                    seq, class_del, count_del = "X" * int(cigar[idx:c.start()]), 1, int(cigar[idx:c.start()])
                else:
                    if class_del == 1:
                        seq_init, seq_end = (seq_init - count_del), (seq_end - count_del)
                    seq, class_del = sequence[seq_init:seq_end], 0

                segments.append((int(pos), end, c.group(), seq))
                pos, idx, seq_init = (int(pos) + int(cigar[idx:c.start()])), (c.start() + 1), seq_end

            # 2. Process the data: remove S and H categories (clipped sequences) and sort by start
            segments = sorted((s for s in segments if s[2] not in ('S', 'H')), key=lambda s: s[0])
            startN = segments[0][0]  # Correct the start of the consensus sequence to the first position defined
            # (L'ajust d'insercions de l'original nomes actuava entre lectures diferents; aqui sempre es una sola lectura)

            # 3. Generate consensus sequence (la primera base per posicio guanya, com el drop_duplicates original)
            consensus, first, last = {}, None, None
            for ini, end, cls, seq in segments:
                ini, end = ini - startN, end - startN
                for j in range(ini, end + 1):
                    k = j - ini
                    if 0 <= k < len(seq):
                        if first is None:
                            first = j
                        last = j
                        consensus.setdefault(j, seq[k])

            # Add 'N' in the positions that are not defined.
            for q in range(first, last):
                consensus.setdefault(q, 'N')

            consensus_seq = ''.join(consensus[q] for q in sorted(consensus))

        spaces = "-" * (pos_end - 1)
        ending = "-" * (3 - (len(spaces + consensus_seq) % 3))
        fasta_seq = '>' + node + '\n' + spaces + consensus_seq + ending + '\n'
        if (isit == 1):
            all_fastas += fasta_seq
    return all_fastas