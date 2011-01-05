#!/usr/bin/env python
"""Looks for the given primer sequences and clips matching reads.


Takes seven command line options, input FASTQ read filename, input primer FASTA
filename, type of primers (forward, reverse or reverse-complement), number of
mismatches (keep this small, for example 1 or 2), minimum length to keep a
read (after primer trimming), should primer-less reads be kept, and finally the
output FASTQ filename.

This script is copyright 2011 by Peter Cock, SCRI, UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script.
"""
import sys
from galaxy_utils.sequence.fasta import fastaReader
from galaxy_utils.sequence.fastq import fastqReader, fastqWriter

def stop_err(msg, err=1):
    sys.stderr.write(msg)
    sys.exit(err)

#Parse Command Line
try:
    input_fastq, primer_fasta, primer_type, mm, min_len, keep_negatives, output_fastq = sys.argv[1:]
except ValueError:
    stop_err("Expected seven arguments, got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))

try:
    mm = int(mm)
except ValueError:
    stop_err("Expected non-negative integer number of mismatches (e.g. 0 or 1), not %r" % mm)
if mm < 0:
    stop_err("Expected non-negtive integer number of mismatches (e.g. 0 or 1), not %r" % mm)
if mm != 0:
    raise NotImplementedError

try:
    min_len = int(min_len)
except ValueError:
    stop_err("Expected non-negative integer min_len (e.g. 0 or 1), not %r" % min_len)
if min_len < 0:
    stop_err("Expected non-negtive integer min_len (e.g. 0 or 1), not %r" % min_len)


if keep_negatives.lower() in ["true", "yes", "on"]:
    keep_negatives = True
elif keep_negatives.lower() in ["false", "no", "off"]:
    keep_negatives = False
else:
    stop_err("Expected boolean for keep_negatives (e.g. true or false), not %r" % keep_negatives)


if primer_type.lower() == "forward":
    forward = True
    rc = False
elif primer_type.lower() == "reverse":
    forward = False
    rc = False
elif primer_type.lower() == "reverse-complement":
    forward = False
    rc = True
else:
    stop_err("Expected foward, reverse or reverse-complement not %r" % primer_type)

#This includes all the IUPAC codes:
_dna_comp_table = '\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TVGHEFCDIJMLKNOPQYSAUBWXRZ[\\]^_`tvghefcdijmlknopqysaubwxrz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff'
def reverse_complement_dna(seq):
    return seq.translate(_dna_comp_table)[::-1]

#Read primer file and record all specified identifiers
primers = set()
in_handle = open(primer_fasta, "rU")
reader = fastaReader(in_handle)
count = 0
for record in reader:
    if rc:
        seq = reverse_complement(record.sequence)
    else:
        seq = record.sequence
    primers.add(seq)
    count += 1
    #TODO - generate mismatches
    #TODO - generate variants for IUPAC ambiguity codes
in_handle.close()
print "%i primer sequences giving %i unique matches" % (count, len(primers))


in_handle = open(input_fastq, "rU")
out_handle = open(output_fastq, "w")
reader = fastqReader(in_handle)
writer = fastqWriter(out_handle)
short = 0
clipped = 0
negs = 0
if forward:
    for record in reader:
        seq = record.sequence
        match = False
        for primer in primers:
            i = seq.find(primer)
            if i != -1:
                match=True
                cut = i + len(primer)
                record.sequence = seq[cut:]
                record.quality = record.quality[cut:]
                if len(record.sequence) >= min_len:
                    clipped += 1
                    writer.write(record)
                else:
                    short += 1
                break
        if not match and keep_negatives:
            negs += 1
            writer.write(record)
else:
    for record in reader:
        seq = record.sequence
        match = False
        for primer in primers:
            i = seq.find(primer)
            if i != -1:
                match=True
                record.sequence = seq[:i]
                record.quality = record.quality[:i]
                if len(record.sequence) >= min_len:
                    clipped += 1
                    writer.write(record)
                else:
                    short += 1
                break
        if not match and keep_negatives:
            negs += 1
            writer.write(record)
in_handle.close()
out_handle.close()

print "Kept %i clipped reads" % clipped
print "Discarded %i short reads" % short
if keep_negatives:
    print "Kept %i non-matching reads" % negs
