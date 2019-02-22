#!/usr/bin/env python
"""Record sequence composition from FASTA, FASTQ or SFF files.

This tool is a short Python script which requires Biopython 1.62 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

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

# Parse Command Line
usage = """Example usage:

$ python seq_composition.py -o my_output.tsv -q input1.fastq -q input2.fastq

At least one input sequence file is required (using the -f, -q, or -s options).
If the expected alphabet is given, the sequence composition is verfied against
it.
"""
# TODO - Case senstivity?
# TODO - GenBank / EMBL input? Needs the datatype defined...
# TODO - Handle all the FASTQ datatype subclasses in the XML cheetah code?
parser = OptionParser(usage=usage)
parser.add_option(
    "-f",
    "--fasta",
    dest="fasta",
    action="append",
    default=[],
    help="Input sequence filename in FASTA format",
)
parser.add_option(
    "-q",
    "--fastq",
    "--fastqsanger",
    "--fastqillumina",
    "--fastqsolexa",
    dest="fastq",
    action="append",
    default=[],
    help="Input sequence filename in FASTQ format",
)
parser.add_option(
    "-s",
    "--sff",
    dest="sff",
    action="append",
    default=[],
    help="Input sequence filename in SFF format",
)
parser.add_option(
    "-o",
    "--output",
    dest="output",
    default=None,
    help="Output filename (tabular)",
    metavar="FILE",
)
parser.add_option(
    "-v",
    "--version",
    dest="version",
    default=False,
    action="store_true",
    help="Show version and quit",
)
options, args = parser.parse_args()

if options.version:
    print("v0.0.3")
    sys.exit(0)

if not (options.fasta or options.fastq or options.sff):
    sys.exit("Requires an input filename")
if not options.output:
    sys.exit("Requires an output filename")


file_count = 0
seq_count = 0
counts = dict()

for format, filenames in [
    ("fasta", options.fasta),
    ("fastq", options.fastq),
    ("sff-trim", options.sff),
]:
    for filename in filenames:
        file_count += 1
        for record in SeqIO.parse(filename, format):
            seq_count += 1
            for letter in record:
                try:
                    counts[letter] += 1
                except KeyError:
                    counts[letter] = 1

total = sum(counts.values())
sys.stderr.write(
    "Counted %i sequence letters from %i records from %i files\n"
    % (total, seq_count, file_count)
)

scale = 100.0 / total
with open(options.output, "w") as handle:
    handle.write("Letter\tCount\tPercentage\n")
    for letter, count in sorted(counts.items()):
        handle.write("%s\t%i\t%0.2f\n" % (letter, count, count * scale))
