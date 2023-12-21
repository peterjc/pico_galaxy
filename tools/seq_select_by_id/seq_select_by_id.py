#!/usr/bin/env python
"""Select FASTA, QUAL, FASTQ or SSF sequences by IDs from a tabular file.

Takes five command line options, tabular filename, ID column number (using
one based counting), input filename, input type (e.g. FASTA or SFF) and the
output filename (same format as input sequence file).

When selecting from an SFF file, any Roche XML manifest in the input file is
preserved in both output files.

This tool is a short Python script which requires Biopython 1.54 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2011-2023 by Peter Cock, The James Hutton Institute UK.
All rights reserved. See accompanying text file for licence details (MIT
license).
"""

from __future__ import print_function

import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.15")
    sys.exit(0)

# Parse Command Line
try:
    tabular_file, col_arg, in_file, seq_format, out_file = sys.argv[1:]
except ValueError:
    sys.exit(
        "Expected five arguments, got %i:\n%s" % (len(sys.argv) - 1, " ".join(sys.argv))
    )
try:
    if col_arg.startswith("c"):
        column = int(col_arg[1:]) - 1
    else:
        column = int(col_arg) - 1
except ValueError:
    sys.exit("Expected column number, got %s" % col_arg)

if seq_format == "fastqcssanger":
    sys.exit("Colorspace FASTQ not supported.")
elif seq_format.lower() in ["sff", "fastq", "qual", "fasta"]:
    seq_format = seq_format.lower()
elif seq_format.lower().startswith("fastq"):
    # We don't care how the qualities are encoded
    seq_format = "fastq"
elif seq_format.lower().startswith("qual"):
    # We don't care what the scores are
    seq_format = "qual"
else:
    sys.exit("Unrecognised file format %r" % seq_format)


try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Biopython 1.54 or later is required")


def parse_ids(tabular_file, col):
    """Read tabular file and record all specified identifiers.

    Will print a single warning to stderr if any of the fields have
    non-trailing white space (only the first word will be used as
    the identifier).
    """
    handle = open(tabular_file)
    warn = False
    for line in handle:
        if line.strip() and not line.startswith("#"):
            field = line.rstrip("\n").split("\t")[col].strip()
            parts = field.split(None, 1)
            if len(parts) > 1 and not warn:
                warn = (
                    "WARNING: Some of your identifiers had white space in them, "
                    + "using first word only. e.g.:\n%s\n" % field
                )
            yield parts[0]
    handle.close()
    if warn:
        sys.stderr.write(warn)


# Index the sequence file.
# If very big, could use SeqIO.index_db() to avoid memory bottleneck...
records = SeqIO.index(in_file, seq_format)
print("Indexed %i sequences" % len(records))

if seq_format.lower() == "sff":
    # Special case to try to preserve the XML manifest
    try:
        from Bio.SeqIO.SffIO import SffWriter
    except ImportError:
        sys.exit("Requires Biopython 1.54 or later")

    try:
        from Bio.SeqIO.SffIO import ReadRocheXmlManifest
    except ImportError:
        # Prior to Biopython 1.56 this was a private function
        from Bio.SeqIO.SffIO import _sff_read_roche_index_xml as ReadRocheXmlManifest

    in_handle = open(in_file, "rb")  # must be binary mode!
    try:
        manifest = ReadRocheXmlManifest(in_handle)
    except ValueError:
        manifest = None
    in_handle.close()

    out_handle = open(out_file, "wb")
    writer = SffWriter(out_handle, xml=manifest)
    count = 0
    # This does have the overhead of parsing into SeqRecord objects,
    # but doing the header and index at the low level is too fidly.
    name = None  # We want the variable to leak from the iterator's scope...
    iterator = (records[name] for name in parse_ids(tabular_file, column))
    try:
        count = writer.write_file(iterator)
    except KeyError:
        out_handle.close()
        if name not in records:
            sys.exit("Identifier %r not found in sequence file" % name)
        else:
            raise
    out_handle.close()
else:
    # Avoid overhead of parsing into SeqRecord objects,
    # just re-use the original formatting from the input file.
    out_handle = open(out_file, "wb")
    count = 0
    for name in parse_ids(tabular_file, column):
        try:
            out_handle.write(records.get_raw(name))
        except KeyError:
            out_handle.close()
            sys.exit("Identifier %r not found in sequence file" % name)
        count += 1
    out_handle.close()

print("Selected %i sequences by ID" % count)
