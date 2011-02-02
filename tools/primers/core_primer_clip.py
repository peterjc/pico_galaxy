"""Shared code for primer/barcode trimming.

This script is copyright 2011 by Peter Cock, SCRI, UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script. Currently it uses Python's regular
expression engine for finding the primers, which for my needs is fast enough.
"""
import sys
import re
from galaxy_utils.sequence.fasta import fastaReader

def stop_err(msg, err=1):
    sys.stderr.write(msg)
    sys.exit(err)

#This includes all the IUPAC codes:
_dna_comp_table = '\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TVGHEFCDIJMLKNOPQYSAUBWXRZ[\\]^_`tvghefcdijmlknopqysaubwxrz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff'
def reverse_complement_dna(seq):
    return seq.translate(_dna_comp_table)[::-1]

ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "ACM",
    "R": "AGR",
    "W": "ATW",
    "S": "CGS",
    "Y": "CTY",
    "K": "GTK",
    "V": "ACGMRSV",
    "H": "ACTMWYH",
    "D": "AGTRWKD",
    "B": "CGTSYKB",
    "X": ".", #faster than [GATCMRWSYKVVHDBXN] or even [GATC]
    "N": ".",
    }

ambiguous_dna_re = {}
for letter, values in ambiguous_dna_values.iteritems():
    if len(values) == 1:
        ambiguous_dna_re[letter] = values
    else:
        ambiguous_dna_re[letter] = "[%s]" % values


def make_reg_ex(seq):
    return "".join(ambiguous_dna_re[letter] for letter in seq)

def make_reg_ex_mm(seq, mm):
    if mm > 2:
        raise NotImplementedError("At most 2 mismatches allowed!")
    seq = seq.upper()
    yield make_reg_ex(seq)
    if mm >= 1:
        for i,letter in enumerate(seq):
            #We'll use a set to remove any duplicate patterns
            #if letter not in "NX":
            pattern = seq[:i] + "N" + seq[i+1:]
            assert len(pattern) == len(seq), "Len %s is %i, len %s is %i" \
                   % (pattern, len(pattern), seq, len(seq))
            yield make_reg_ex(pattern)
    if mm >=2:
        for i,letter in enumerate(seq):
            #We'll use a set to remove any duplicate patterns
            #if letter not in "NX":
            for k,letter in enumerate(seq[i+1:]):
                #We'll use a set to remove any duplicate patterns
                #if letter not in "NX":
                pattern = seq[:i] + "N" + seq[i+1:i+1+k] + "N" + seq[i+k+2:]
                assert len(pattern) == len(seq), "Len %s is %i, len %s is %i" \
                       % (pattern, len(pattern), seq, len(seq))
                yield make_reg_ex(pattern)

def load_primers_as_re(primer_fasta, mm, rc=False):
    #Read primer file and record all specified sequences
    primers = set()
    in_handle = open(primer_fasta, "rU")
    reader = fastaReader(in_handle)
    count = 0
    for record in reader:
        if rc:
            seq = reverse_complement(record.sequence)
        else:
            seq = record.sequence
        #primers.add(re.compile(make_reg_ex(seq)))
        count += 1
        for pattern in make_reg_ex_mm(seq, mm):
            primers.add(pattern)
    in_handle.close()
    return count, re.compile("|".join(sorted(set(primers)))) #make one monster re!
