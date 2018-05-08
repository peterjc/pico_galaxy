#!/usr/bin/env python
"""Compute length of FASTA, QUAL, FASTQ or SSF sequences.

Takes three command line options: input sequence filename, input type
(e.g. FASTA or SFF) and the output filename (tabular).

This tool is a short Python script which requires Biopython 1.54 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2018 by Peter Cock, The James Hutton Institute UK.
All rights reserved. See accompanying text file for licence details (MIT
license).
"""

from __future__ import print_function

import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Missing required Python library Biopython.")


# Parse Command Line
try:
    in_file, seq_format, out_file = sys.argv[1:]
except ValueError:
    sys.exit("Expected three arguments (input file, format, output file), "
             "got %i:\n%s" % (len(sys.argv) - 1, " ".join(sys.argv)))


if seq_format.startswith("fastq"):
    # We don't care about the quality score encoding, just
    # need to translate Galaxy format name into something
    # Biopython will accept:
    format = "fastq"
elif seq_format.lower() == "csfasta":
    # I have not tested with colour space FASTA
    format = "fasta"
elif seq_format.lower == "sff":
    # The masked/trimmed numbers are more interesting
    format = "sff-trim"
elif seq_format.lower() in ["fasta", "qual"]:
    format = seq_format.lower()
else:
    # TODO: Does Galaxy understand GenBank, EMBL, etc yet?
    sys.exit("Unexpected format argument: %r" % seq_format)


count = 0
total = 0
with open(out_file, "w") as out_handle:
    out_handle.write("#Identifier\tLength\n")
    for record in SeqIO.parse(in_file, format):
        count += 1
        length = len(record)
        total += length
        out_handle.write("%s\t%i\n" % (record.id, length))
print("%i sequences, total length %i" % (count, total))
