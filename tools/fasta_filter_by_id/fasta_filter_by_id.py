#!/usr/bin/env python
"""Filter a FASTA file with IDs from a tabular file, e.g. from BLAST.

NOTE - This script is now OBSOLETE, having been replaced by a new verion
which handles FASTA, FASTQ and SFF all in one.

Takes five command line options, tabular filename, ID column numbers
(comma separated list using one based counting), input FASTA filename, and
two output FASTA filenames (for records with and without the given IDs).

If either output filename is just a minus sign, that file is not created.
This is intended to allow output for just the matched (or just the non-matched)
records.

Note in the default NCBI BLAST+ tabular output, the query sequence ID is
in column one, and the ID of the match from the database is in column two.
Here sensible values for the column numbers would therefore be "1" or "2".

This tool is copyright 2010-2023 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See accompanying text file for licence details (MIT license).
"""

from __future__ import print_function

import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.7")
    sys.exit(0)

from galaxy_utils.sequence.fasta import fastaReader, fastaWriter

# Parse Command Line
try:
    tabular_file, cols_arg, in_file, out_positive_file, out_negative_file = sys.argv[1:]
except ValueError:
    sys.exit(
        "Expected five arguments, got %i:\n%s" % (len(sys.argv) - 1, " ".join(sys.argv))
    )
try:
    columns = [int(arg) - 1 for arg in cols_arg.split(",")]
except ValueError:
    sys.exit("Expected list of columns (comma separated integers), got %s" % cols_arg)

# Read tabular file and record all specified identifiers
ids = set()
handle = open(tabular_file)
if len(columns) > 1:
    # General case of many columns
    for line in handle:
        if line.startswith("#"):
            # Ignore comments
            continue
        parts = line.rstrip("\n").split("\t")
        for col in columns:
            ids.add(parts[col])
    print("Using %i IDs from %i columns of tabular file" % (len(ids), len(columns)))
else:
    # Single column, special case speed up
    col = columns[0]
    for line in handle:
        if not line.startswith("#"):
            ids.add(line.rstrip("\n").split("\t")[col])
    print("Using %i IDs from tabular file" % (len(ids)))
handle.close()

# Write filtered FASTA file based on IDs from tabular file
reader = fastaReader(open(in_file))
if out_positive_file != "-" and out_negative_file != "-":
    print("Generating two FASTA files")
    positive_writer = fastaWriter(open(out_positive_file, "w"))
    negative_writer = fastaWriter(open(out_negative_file, "w"))
    for record in reader:
        # The [1:] is because the fastaReader leaves the > on the identifer.
        if record.identifier and record.identifier.split()[0][1:] in ids:
            positive_writer.write(record)
        else:
            negative_writer.write(record)
    positive_writer.close()
    negative_writer.close()
elif out_positive_file != "-":
    print("Generating matching FASTA file")
    positive_writer = fastaWriter(open(out_positive_file, "w"))
    for record in reader:
        # The [1:] is because the fastaReader leaves the > on the identifer.
        if record.identifier and record.identifier.split()[0][1:] in ids:
            positive_writer.write(record)
    positive_writer.close()
elif out_negative_file != "-":
    print("Generating non-matching FASTA file")
    negative_writer = fastaWriter(open(out_negative_file, "w"))
    for record in reader:
        # The [1:] is because the fastaReader leaves the > on the identifer.
        if not record.identifier or record.identifier.split()[0][1:] not in ids:
            negative_writer.write(record)
    negative_writer.close()
else:
    sys.exit("Neither output file requested")
reader.close()
