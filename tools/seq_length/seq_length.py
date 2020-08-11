#!/usr/bin/env python
"""Compute length of FASTA, QUAL, FASTQ or SSF sequences.

Takes three command line options: input sequence filename, input type
(e.g. FASTA or SFF) and the output filename (tabular).

This tool is a short Python script which requires Biopython 1.54 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2018 by Peter Cock, The James Hutton Institute UK.
All rights reserved. See accompanying text file for licence details (MIT
license).
"""

from __future__ import print_function

import sys
from collections import defaultdict
from optparse import OptionParser

usage = r"""Use as follows to compute all the lengths in a sequence file:

$ python seq_length.py -i example.fasta -f fasta -o lengths.tsv
"""

parser = OptionParser(usage=usage)
parser.add_option(
    "-i",
    "--input",
    dest="input",
    default=None,
    help="Input sequence filename (FASTA, FASTQ, etc)",
    metavar="FILE",
)
parser.add_option(
    "-f",
    "--format",
    dest="format",
    default=None,
    help="Input sequence format (FASTA, QUAL, FASTQ, SFF)",
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
    "-s",
    "--stats",
    dest="stats",
    default=False,
    action="store_true",
    help="Compute statistics (median, N50 - will require much more RAM).",
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
    print("v0.0.4")
    sys.exit(0)
if not options.input:
    sys.exit("Require an input filename")
if not options.format:
    sys.exit("Require the input format")
if not options.output:
    sys.exit("Require an output filename")

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Missing required Python library Biopython.")

try:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError:
    sys.exit("Biopython tool old?, missing Bio.SeqIO.QualityIO.FastqGeneralIterator")

try:
    from Bio.SeqIO.FastaIO import SimpleFastaParser
except ImportError:
    sys.exit("Biopython tool old?, missing Bio.SeqIO.FastaIO.SimpleFastaParser")

in_file = options.input
out_file = options.output

if options.format.startswith("fastq"):
    # We don't care about the quality score encoding, just
    # need to translate Galaxy format name into something
    # Biopython will accept:
    format = "fastq"
elif options.format.lower() == "csfasta":
    # I have not tested with colour space FASTA
    format = "fasta"
elif options.format.lower() == "sff":
    # The masked/trimmed numbers are more interesting
    format = "sff-trim"
elif options.format.lower() in ["fasta", "qual"]:
    format = options.format.lower()
else:
    # TODO: Does Galaxy understand GenBank, EMBL, etc yet?
    sys.exit("Unexpected format argument: %r" % options.format)


def median_from_counts_dict(counts_dict, count=None):
    sorted_lengths = sorted(counts_dict)  # i.e. sort the keys
    if count is None:
        count = sum(counts_dict.values())
    index = count / 2
    if count % 2:
        # Odd, easy case - will be an exact value
        # within one of the tally entries
        for value in sorted_lengths:
            index -= counts_dict[value]
            if index < 0:
                return value
    else:
        # Even, hard case - may have to take mean
        special = None
        for value in sorted_lengths:
            if special is not None:
                # We were right at boundary
                return (special + value) / 2.0
            index -= counts_dict[value]
            if index == 0:
                # Special case, want mean of this value
                # (final entry of this tally) and the next
                # value (first entry of the next tally).
                special = value
            elif index < 0:
                # Typical case, the two middle values
                # are equal and fall into same dict entry
                return value
    return None


if sys.version_info[0] >= 3:
    from statistics import median

    for test in [
        {1: 4, 2: 3},
        {1: 4, 3: 6},
        {1: 4, 5: 4},
        {0: 5, 1: 1, 2: 1, 3: 5},
        {0: 5, 1: 1, 2: 1, 3: 1, 4: 5},
    ]:
        test_list = []
        for v, c in test.items():
            test_list.extend([v] * c)
        # print(test)
        # print(test_list)
        assert median_from_counts_dict(test) == median(test_list)


def n50_from_counts_dict(counts_dict):
    """Calculate N50.

    N50 is a statistical measure of average length of a set of sequences.
    It is used widely in genomics, especially in reference to contig or
    supercontig lengths within a draft assembly.

    Given a set of sequences of varying lengths, the N50 length is defined
    as the length N for which 50% of all bases in the sequences are in a
    sequence of length L < N. This can be found mathematically as follows:
    Take a list L of positive integers. Create another list L' , which is
    identical to L, except that every element n in L has been replaced with
    n copies of itself. Then the median of L' is the N50 of L. For example:
    If L = {2, 2, 2, 3, 3, 4, 8, 8}, then L' consists of six 2's, six 3's,
    four 4's, and sixteen 8's; the N50 of L is the median of L' , which is 6.

    https://web.archive.org/web/20160726124802/http://www.broadinstitute.org/crd/wiki/index.php/N50
    """
    # Continuing the above example, input L would be {2:3, 3:2, 4:1, 8:2}
    # and L' becomes {2:6, 3:6, 4:4, 8:16}} as tally tables.
    l_prime = {v: v * c for v, c in counts_dict.items()}
    return median_from_counts_dict(l_prime)


count = 0
total = 0
stats = bool(options.stats)
length_counts = defaultdict(int)  # used if stats requested
length_min = sys.maxsize  # used if stats not requested
length_max = 0

with open(out_file, "w") as out_handle:
    out_handle.write("#Identifier\tLength\n")
    if format == "fastq":
        with open(in_file) as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                count += 1
                length = len(seq)
                total += length
                identifier = title.split(None, 1)[0]
                out_handle.write("%s\t%i\n" % (identifier, length))
                if stats:
                    length_counts[length] += 1
                else:
                    length_min = min(length_min, length)
                    length_max = max(length_max, length)
    elif format == "fasta":
        with open(in_file) as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                count += 1
                length = len(seq)
                total += length
                identifier = title.split(None, 1)[0]
                out_handle.write("%s\t%i\n" % (identifier, length))
                if stats:
                    length_counts[length] += 1
                else:
                    length_min = min(length_min, length)
                    length_max = max(length_max, length)
    else:
        for record in SeqIO.parse(in_file, format):
            count += 1
            length = len(record)
            total += length
            out_handle.write("%s\t%i\n" % (record.id, length))
            if stats:
                length_counts[length] += 1
            else:
                length_min = min(length_min, length)
                length_max = max(length_max, length)
print(
    "%i sequences, total length %i, mean %0.1f" % (count, total, float(total) / count)
)
if not count:
    pass
elif not stats:
    print("Shortest %i, longest %i" % (length_min, length_max))
elif count and stats:
    print("Shortest %i, longest %i" % (min(length_counts), max(length_counts)))
    median = median_from_counts_dict(length_counts, count)
    n50 = n50_from_counts_dict(length_counts)
    print("Median length %0.1f, N50 %i" % (median, n50))
