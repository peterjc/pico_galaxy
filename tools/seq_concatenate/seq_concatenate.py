#!/usr/bin/env python
"""Concatenate sequences in GenBank format.

This tool is a short Python script which requires Biopython 1.64 or later.
If you use this tool in scientific work leading to a publication, please
cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2014 by Peter Cock, The James Hutton Institute
(formerly the Scottish Crop Research Institute, SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT license).

Use -v or --version to get the version, -h or --help for help.
"""
import sys
from optparse import OptionParser

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# Parse Command Line
usage = """Example usage:

$ python seq_concatenate.py -o my_output.gbk -n 20 file1.gbk file2.gbk

At least one input sequence file is required (in GenBank format).
"""
# TODO - FASTA / GenBank / EMBL input? Needs the datatype defined...
# TODO - Option for number of Ns to use as linker
# TODO - Option for explicit linker sequence
# TODO - Support protein input?
parser = OptionParser(usage=usage)
parser.add_option('-o', '--output', dest='output',
                  default=None, help='Output filename (tabular)',
                  metavar="FILE")
parser.add_option("-v", "--version", dest="version",
                  default=False, action="store_true",
                  help="Show version and quit")
options, args = parser.parse_args()

if options.version:
    print("v0.0.1")
    sys.exit(0)

if not args:
    sys.exit("Require an input filename")
if not options.output:
    sys.exit("Require an output filename")


file_count = 0
seq_count = 0
counts = dict()

# TODO - itertools chain?


def get_records(filenames):
    for f in filenames:
        for r in SeqIO.parse(f, "genbank"):
            yield r


gap = 1000
spacer = SeqRecord(Seq("N" * gap, generic_dna))
spacer.features.append(SeqFeature(FeatureLocation(0, len(spacer)), type="gap"))

record = None
count = 0
for r in get_records(args):
    if record is None:
        record = r
    else:
        record += spacer
        record += r
    count += 1
record.id = "union"
record.description = "Concatenation of %i records" % count

x = SeqIO.write(record, options.output, "genbank")
assert x == 1, "Error writing record"
print("Done, merged %i records into a single record of length %i" % (count, len(record)))
del record
