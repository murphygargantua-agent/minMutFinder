#!/usr/bin/env python3

# Optimized reads_assembler: uses numpy arrays instead of pandas for speed.

import re
import numpy as np

def cigar_remove(cigar, startH, new_cig, idx, number_class, H_exists, p):
    """Simplified cigar parser - returns cleaned cigar string."""
    # Parse cigar into (number, letter) pairs
    pairs = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    return pairs

def reads_assembler(samfile, all_fastas=''):
    """Read SAM file and produce multi-FASTA of per-read consensus sequences.
    
    Each read is converted to a CIGAR-aware consensus aligned to the reference.
    """
    try:
        with open(samfile, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        return all_fastas
    
    all_fastas_list = []
    
    for line in lines:
        line = line.rstrip('\n')
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) < 11:
            continue
        
        node = fields[0]
        flag = int(fields[1])
        pos = int(fields[3])
        cigar = fields[5]
        sequence = fields[9]
        
        if cigar == '*' or not sequence:
            continue
        
        isit = 1 if flag >= 0 else 0
        pos_end = pos + len(sequence)
        
        # Parse CIGAR into (number, type) pairs
        cig_pairs = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        if not cig_pairs:
            continue
        
        # Build per-position consensus using a dict
        ref_pos = pos
        seq_idx = 0
        consensus = {}  # ref_pos -> base
        
        for num, ctype in cig_pairs:
            n = int(num)
            if ctype in ('M', '=', 'X'):
                for i in range(n):
                    if seq_idx + i < len(sequence):
                        consensus[ref_pos + i] = sequence[seq_idx + i]
                ref_pos += n
                seq_idx += n
            elif ctype == 'D':
                ref_pos += n
            elif ctype == 'I':
                seq_idx += n
            elif ctype == 'S':
                seq_idx += n
            elif ctype == 'H':
                pass  # clip, skip
            elif ctype == 'N':
                ref_pos += n
        
        if not consensus:
            continue
        
        min_pos = min(consensus.keys())
        max_pos = max(consensus.keys())
        
        # Build consensus string
        cons_list = [consensus.get(p, 'N') for p in range(min_pos, max_pos + 1)]
        consensus_seq = ''.join(cons_list)
        
        # Pad to multiple of 3
        spaces = '-' * (min_pos - 1)
        ending = '-' * ((3 - (len(spaces + consensus_seq) % 3)) % 3)
        fasta_seq = '>' + node + '\n' + spaces + consensus_seq + ending + '\n'
        all_fastas_list.append(fasta_seq)
    
    all_fastas += ''.join(all_fastas_list)
    return all_fastas
