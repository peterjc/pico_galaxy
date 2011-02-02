#!/usr/bin/env python
"""Looks for the given primer sequences and clips matching FASTQ reads.

Takes seven command line options, input FASTQ read filename, input primer FASTA
filename, type of primers (forward, reverse or reverse-complement), number of
mismatches (currently only 0, 1 and 2 are supported), minimum length to keep a
read (after primer trimming), should primer-less reads be kept, and finally the
output FASTQ filename.

This can also be used for stripping off (and optionally filtering on) barcodes.

This script is copyright 2011 by Peter Cock, SCRI, UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script.
"""
import sys
from core_primer_clip import stop_err, load_primers_as_re
from galaxy_utils.sequence.fastq import fastqReader, fastqWriter

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
if mm not in [0,1,2]:
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


#Read primer file and record all specified sequences
count, primer = load_primers_as_re(primer_fasta, mm, rc)
print "%i primer sequences" % count

in_handle = open(input_fastq, "rU")
out_handle = open(output_fastq, "w")
reader = fastqReader(in_handle)
writer = fastqWriter(out_handle)
short = 0
clipped = 0
negs = 0
if forward:
    for record in reader:
        seq = record.sequence.upper()
        result = primer.search(seq)
        if result:
            #Forward primer, take everything after it
            cut = result.end()
            record.sequence = seq[cut:]
            record.quality = record.quality[cut:]
            if len(record.sequence) >= min_len:
                clipped += 1
                writer.write(record)
            else:
                short += 1
        elif keep_negatives:
            negs += 1
            writer.write(record)
else:
    for record in reader:
        seq = record.sequence.upper()
        result = primer.search(seq)
        if result:
            #Reverse primer, take everything before it
            cut = result.start()
            record.sequence = seq[cut:]
            record.quality = record.quality[cut:]
            if len(record.sequence) >= min_len:
                clipped += 1
                writer.write(record)
            else:
                short += 1
        elif keep_negatives:
            negs += 1
            writer.write(record)
in_handle.close()
out_handle.close()

print "Kept %i clipped reads" % clipped
print "Discarded %i short reads" % short
if keep_negatives:
    print "Kept %i non-matching reads" % negs
