#!/usr/bin/env python3
"""CDS mapping utilities for minMutFinder CMV analysis.

Maps genome (NC_006273.2) coordinates to CDS coordinates and translates
variants to amino acid changes using the correct coding sequences.
"""

import json
import bisect
import os

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]

def load_cds_coords(path):
    """Load CDS coordinates from JSON file.
    
    Returns list of dicts with keys: gene, start, end, cdna_start, cdna_end, strand, length.
    All positions are 1-based inclusive (genome coordinates).
    """
    with open(path) as f:
        coords = json.load(f)
    # Sort by start position for binary search
    coords.sort(key=lambda x: x["start"])
    return coords

def find_cds(genome_pos, cds_coords):
    """Find the CDS that contains the given genome position (1-based).
    
    genome_pos: 1-based genome coordinate.
    Returns CDS dict or None.
    """
    starts = [c["start"] for c in cds_coords]
    idx = bisect.bisect_right(starts, genome_pos) - 1
    if idx >= 0:
        cds = cds_coords[idx]
        if cds["start"] <= genome_pos <= cds["end"]:
            return cds
    return None

def genome_to_cds_pos(genome_pos, cds):
    """Convert genome position (1-based) to CDS position (1-based).
    
    For forward strand: cds_pos = genome_pos - start + 1
    For reverse strand: cds_pos = end - genome_pos + 1
    """
    if cds["strand"] == 1:
        return genome_pos - cds["start"] + 1
    else:
        return cds["end"] - genome_pos + 1

def get_cds_seq(cds, cds_fasta):
    """Get the CDS sequence (forward strand, 5'->3' as it would be translated)."""
    gene = cds["gene"]
    if gene in cds_fasta:
        seq = cds_fasta[gene]
        if cds["strand"] == -1:
            seq = revcomp(seq)
        return seq
    return None

def load_cds_fasta(path):
    """Load CDS FASTA file.
    
    Returns dict: gene_name -> sequence (forward strand, 5'->3').
    """
    sequences = {}
    current_gene = None
    current_seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current_gene is not None:
                    sequences[current_gene] = "".join(current_seq)
                current_gene = line[1:].strip().split("|")[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_gene is not None:
        sequences[current_gene] = "".join(current_seq)
    return sequences

def build_position_lookup(cds_coords, cds_fasta):
    """Build a lookup: genome_pos (1-based) -> (gene, cds_pos_1based, strand, ref_cds_seq, cds_full_seq).
    
    For reverse strand, cds_pos is the position from the 3' end of the CDS (i.e., 
    the position in the forward-strand CDS sequence).
    """
    lookup = {}
    for cds in cds_coords:
        gene = cds["gene"]
        strand = cds["strand"]
        seq = get_cds_seq(cds, cds_fasta)
        if seq is None:
            continue
        length = cds["end"] - cds["start"] + 1  # 1-based inclusive
        for i in range(length):
            genome_pos = cds["start"] + i  # 1-based genome position
            cds_pos = i + 1  # 1-based CDS position
            lookup[genome_pos] = {
                "gene": gene,
                "cds_pos": cds_pos,
                "strand": strand,
                "seq": seq,
            }
    return lookup

def get_ref_alt_at_genome_pos(genome_pos, ref, alt, cds, cds_fasta):
    """Get the reference and alternate alleles at a genome position,
    already oriented on the CDS forward strand.
    
    For forward strand CDS: REF and ALT are as-is.
    For reverse strand CDS: REF and ALT need to be reverse-complemented.
    """
    if cds["strand"] == 1:
        return ref, alt
    else:
        return revcomp(ref), revcomp(alt)

def aa_pos_from_cds_pos(cds_pos):
    """Calculate amino acid position from CDS position (1-based)."""
    return (cds_pos - 1) // 3 + 1
