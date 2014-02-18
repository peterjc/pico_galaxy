#!/usr/bin/env python
"""Back-translate a protein alignment to nucleotides

This tool is a short Python script (using Biopython library functions) to
load a protein alignment, and matching nucleotide FASTA file of unaligned
sequences, in order to produce a codon aware nucleotide alignment - which
can be viewed as a back translation.

The development repository for this tool is here:

* https://github.com/peterjc/pico_galaxy/tree/master/tools/align_back_trans  

This tool is available with a Galaxy wrapper from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/align_back_trans

See accompanying text file for licence details (MIT licence).

This is version 0.0.2 of the script.
"""

import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.2"
    sys.exit(0)

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

def sequence_back_translate(aligned_protein_record, unaligned_nucleotide_record, gap):
    #TODO - Separate arguments for protein gap and nucleotide gap?
    if not gap or len(gap) != 1:
        raise ValueError("Please supply a single gap character")

    alpha = unaligned_nucleotide_record.seq.alphabet
    if hasattr(alpha, "gap_char"):
        gap_codon = alpha.gap_char * 3
        assert len(gap_codon) == 3
    else:
        from Bio.Alphabet import Gapped
        alpha = Gapped(alpha, gap)
        gap_codon = gap*3

    if len(aligned_protein_record.seq.ungap(gap))*3 != len(unaligned_nucleotide_record.seq):
        stop_err("Inconsistent lengths for %s, ungapped protein %i, "
                 "tripled %i vs ungapped nucleotide %i" %
                 (len(aligned_protein_record.seq.ungap(gap)),
                  len(aligned_protein_record.seq.ungap(gap))*3,
                  len(unaligned_nucleotide_record.seq)))

    seq = []
    nuc = str(unaligned_nucleotide_record.seq)
    for amino_acid in aligned_protein_record.seq:
        if amino_acid == gap:
            seq.append(gap_codon)
        else:
            seq.append(nuc[:3])
            nuc = nuc[3:]
    assert not nuc, "Nucleotide sequence for %r longer than protein %s" \
        % (unaligned_nucleotide_record.id, aligned_protein_record.id)

    aligned_nuc = unaligned_nucleotide_record[:] #copy for most annotation
    aligned_nuc.letter_annotation = {} #clear this
    aligned_nuc.seq = Seq("".join(seq), alpha) #replace this
    assert len(aligned_protein_record.seq) * 3 == len(aligned_nuc)
    return aligned_nuc

def alignment_back_translate(protein_alignment, nucleotide_records, key_function=None, gap=None):
    """Thread nucleotide sequences onto a protein alignment."""
    #TODO - Separate arguments for protein and nucleotide gap characters?
    if key_function is None:
        key_function = lambda x: x
    if gap is None:
        gap = "-"

    aligned = []
    for protein in protein_alignment:
        try:
            nucleotide = nucleotide_records[key_function(protein.id)]
        except KeyError:
            raise ValueError("Could not find nucleotide sequence for protein %r" \
                             % protein.id)
        aligned.append(sequence_back_translate(protein, nucleotide, gap))
    return MultipleSeqAlignment(aligned)


if len(sys.argv) == 4:
    align_format, prot_align_file, nuc_fasta_file = sys.argv[1:]
    nuc_align_file = sys.stdout
elif len(sys.argv) == 5:
    align_format, prot_align_file, nuc_fasta_file, nuc_align_file = sys.argv[1:]
else:
    stop_err("""This is a Python script for 'back-translating' a protein alignment,

It requires three or four arguments:
- alignment format (e.g. fasta, clustal),
- aligned protein file (in specified format),
- unaligned nucleotide file (in fasta format).
- aligned nucleotiode output file (in same format), optional.

The nucleotide alignment is printed to stdout if no output filename is given.

Example usage:

$ python align_back_trans.py fasta demo_prot_align.fasta demo_nucs.fasta demo_nuc_align.fasta

Warning: If the output file already exists, it will be overwritten.

This script is available with sample data and a Galaxy wrapper here:
https://github.com/peterjc/pico_galaxy/tree/master/tools/align_back_trans
http://toolshed.g2.bx.psu.edu/view/peterjc/align_back_trans
""")

prot_align = AlignIO.read(prot_align_file, align_format, alphabet=generic_protein)
nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
nuc_align = alignment_back_translate(prot_align, nuc_dict, gap="-")
AlignIO.write(nuc_align, nuc_align_file, align_format)
