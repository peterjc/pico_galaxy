#!/usr/bin/env python
"""Filter a FASTA, FASTQ or SSF file with IDs from a tabular file.

Takes six command line options, tabular filename, ID column numbers (comma
separated list using one based counting), input filename, input type (e.g.
FASTA or SFF) and up to two output filenames (for records with and without
the given IDs, same format as input sequence file).

When filtering an SFF file, any Roche XML manifest in the input file is
preserved in both output files.

Note in the default NCBI BLAST+ tabular output, the query sequence ID is
in column one, and the ID of the match from the database is in column two.
Here sensible values for the column numbers would therefore be "1" or "2".

This tool is a short Python script which requires Biopython 1.54 or later.
If you use this tool in scientific work leading to a publication, please
cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2010-2017 by Peter Cock, The James Hutton Institute
(formerly the Scottish Crop Research Institute, SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT license).

Use -v or --version to get the version, -h or --help for help.
"""

from __future__ import print_function

import os
import re
import sys

from optparse import OptionParser

# Parse Command Line
usage = """Use as follows:

$ python seq_filter_by_id.py [options] tab1 cols1 [, tab2 cols2, ...]

e.g. Positive matches using column one from tabular file:

$ seq_filter_by_id.py -i my_seqs.fastq -f fastq -p matches.fastq ids.tabular 1

Multiple tabular files and column numbers may be given, or replaced with
the -t or --text option.
"""
parser = OptionParser(usage=usage)
parser.add_option('-i', '--input', dest='input',
                  default=None, help='Input sequences filename',
                  metavar="FILE")
parser.add_option('-f', '--format', dest='format',
                  default=None,
                  help='Input sequence format (e.g. fasta, fastq, sff)')
parser.add_option('-t', '--text', dest='id_list',
                  default=None, help="Lists of white space separated IDs (instead of a tabular file)")
parser.add_option('-p', '--positive', dest='output_positive',
                  default=None,
                  help='Output filename for matches',
                  metavar="FILE")
parser.add_option('-n', '--negative', dest='output_negative',
                  default=None,
                  help='Output filename for non-matches',
                  metavar="FILE")
parser.add_option("-l", "--logic", dest="logic",
                  default="UNION",
                  help="How to combined multiple ID columns (UNION or INTERSECTION)")
parser.add_option("-s", "--suffix", dest="suffix",
                  action="store_true",
                  help="Ignore pair-read suffices for matching names")
parser.add_option("-v", "--version", dest="version",
                  default=False, action="store_true",
                  help="Show version and quit")

options, args = parser.parse_args()

if options.version:
    print("v0.2.7")
    sys.exit(0)

in_file = options.input
seq_format = options.format
out_positive_file = options.output_positive
out_negative_file = options.output_negative
logic = options.logic
drop_suffices = bool(options.suffix)

if in_file is None or not os.path.isfile(in_file):
    sys.exit("Missing input file: %r" % in_file)
if out_positive_file is None and out_negative_file is None:
    sys.exit("Neither output file requested")
if seq_format is None:
    sys.exit("Missing sequence format")
if logic not in ["UNION", "INTERSECTION"]:
    sys.exit("Logic agrument should be 'UNION' or 'INTERSECTION', not %r" % logic)
if options.id_list and args:
    sys.exit("Cannot accept IDs via both -t in the command line, and as tabular files")
elif not options.id_list and not args:
    sys.exit("Expected matched pairs of tabular files and columns (or -t given)")
if len(args) % 2:
    sys.exit("Expected matched pairs of tabular files and columns, not: %r" % args)


# Cope with three widely used suffix naming convensions,
# Illumina: /1 or /2
# Forward/revered: .f or .r
# Sanger, e.g. .p1k and .q1k
# See http://staden.sourceforge.net/manual/pregap4_unix_50.html
# re_f = re.compile(r"(/1|\.f|\.[sfp]\d\w*)$")
# re_r = re.compile(r"(/2|\.r|\.[rq]\d\w*)$")
re_suffix = re.compile(r"(/1|\.f|\.[sfp]\d\w*|/2|\.r|\.[rq]\d\w*)$")
assert re_suffix.search("demo.f")
assert re_suffix.search("demo.s1")
assert re_suffix.search("demo.f1k")
assert re_suffix.search("demo.p1")
assert re_suffix.search("demo.p1k")
assert re_suffix.search("demo.p1lk")
assert re_suffix.search("demo/2")
assert re_suffix.search("demo.r")
assert re_suffix.search("demo.q1")
assert re_suffix.search("demo.q1lk")

identifiers = []
for i in range(len(args) // 2):
    tabular_file = args[2 * i]
    cols_arg = args[2 * i + 1]
    if not os.path.isfile(tabular_file):
        sys.exit("Missing tabular identifier file %r" % tabular_file)
    try:
        columns = [int(arg) - 1 for arg in cols_arg.split(",")]
    except ValueError:
        sys.exit("Expected list of columns (comma separated integers), got %r" % cols_arg)
    if min(columns) < 0:
        sys.exit("Expect one-based column numbers (not zero-based counting), got %r" % cols_arg)
    identifiers.append((tabular_file, columns))

name_warn = False


def check_white_space(name):
    """Check identifier for white space, take first word only."""
    parts = name.split(None, 1)
    global name_warn
    if not name_warn and len(parts) > 1:
        name_warn = "WARNING: Some of your identifiers had white space in them, " + \
                    "using first word only. e.g.:\n%s\n" % name
    return parts[0]


if drop_suffices:
    def clean_name(name):
        """Remove suffix."""
        name = check_white_space(name)
        match = re_suffix.search(name)
        if match:
            # Use the fact this is a suffix, and regular expression will be
            # anchored to the end of the name:
            return name[:match.start()]
        else:
            # Nothing to do
            return name
    assert clean_name("foo/1") == "foo"
    assert clean_name("foo/2") == "foo"
    assert clean_name("bar.f") == "bar"
    assert clean_name("bar.r") == "bar"
    assert clean_name("baz.p1") == "baz"
    assert clean_name("baz.q2") == "baz"
else:
    # Just check the white space
    clean_name = check_white_space


mapped_chars = {
    '>': '__gt__',
    '<': '__lt__',
    "'": '__sq__',
    '"': '__dq__',
    '[': '__ob__',
    ']': '__cb__',
    '{': '__oc__',
    '}': '__cc__',
    '@': '__at__',
    '\n': '__cn__',
    '\r': '__cr__',
    '\t': '__tc__',
    '#': '__pd__',
}

# Read tabular file(s) and record all specified identifiers
ids = None  # Will be a set
if options.id_list:
    assert not identifiers
    ids = set()
    id_list = options.id_list
    # Galaxy turns \r into __cr__ (CR) etc
    for k in mapped_chars:
        id_list = id_list.replace(mapped_chars[k], k)
    for x in options.id_list.split():
        ids.add(clean_name(x.strip()))
    print("Have %i unique identifiers from list" % len(ids))
for tabular_file, columns in identifiers:
    file_ids = set()
    handle = open(tabular_file, "rU")
    if len(columns) > 1:
        # General case of many columns
        for line in handle:
            if line.startswith("#"):
                # Ignore comments
                continue
            parts = line.rstrip("\n").split("\t")
            for col in columns:
                name = clean_name(parts[col])
                if name:
                    file_ids.add(name)
    else:
        # Single column, special case speed up
        col = columns[0]
        for line in handle:
            if not line.strip():  # skip empty lines
                continue
            if not line.startswith("#"):
                name = clean_name(line.rstrip("\n").split("\t")[col])
                if name:
                    file_ids.add(name)
    print("Using %i IDs from column %s in tabular file" % (len(file_ids), ", ".join(str(col + 1) for col in columns)))
    if ids is None:
        ids = file_ids
    if logic == "UNION":
        ids.update(file_ids)
    else:
        ids.intersection_update(file_ids)
    handle.close()
if len(identifiers) > 1:
    if logic == "UNION":
        print("Have %i IDs combined from %i tabular files" % (len(ids), len(identifiers)))
    else:
        print("Have %i IDs in common from %i tabular files" % (len(ids), len(identifiers)))
if name_warn:
    sys.stderr.write(name_warn)


def fasta_filter(in_file, pos_file, neg_file, wanted):
    """FASTA filter producing 60 character line wrapped outout."""
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    pos_count = neg_count = 0
    # Galaxy now requires Python 2.5+ so can use with statements,
    with open(in_file) as in_handle:
        # Doing the if statement outside the loop for speed
        # (with the downside of three very similar loops).
        if pos_file is not None and neg_file is not None:
            print("Generating two FASTA files")
            with open(pos_file, "w") as pos_handle:
                with open(neg_file, "w") as neg_handle:
                    for identifier, record in SimpleFastaParser(in_handle):
                        if clean_name(identifier) in wanted:
                            pos_handle.write(record)
                            pos_count += 1
                        else:
                            neg_handle.write(record)
                            neg_count += 1
        elif pos_file is not None:
            print("Generating matching FASTA file")
            with open(pos_file, "w") as pos_handle:
                for identifier, record in SimpleFastaParser(in_handle):
                    if clean_name(identifier) in wanted:
                        pos_handle.write(record)
                        pos_count += 1
                    else:
                        neg_count += 1
        else:
            print("Generating non-matching FASTA file")
            assert neg_file is not None
            with open(neg_file, "w") as neg_handle:
                for identifier, record in SimpleFastaParser(in_handle):
                    if clean_name(identifier) in wanted:
                        pos_count += 1
                    else:
                        neg_handle.write(record)
                        neg_count += 1
    return pos_count, neg_count


def fastq_filter(in_file, pos_file, neg_file, wanted):
    """FASTQ filter."""
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    handle = open(in_file, "r")
    if pos_file is not None and neg_file is not None:
        print("Generating two FASTQ files")
        positive_handle = open(pos_file, "w")
        negative_handle = open(neg_file, "w")
        print(in_file)
        for title, seq, qual in FastqGeneralIterator(handle):
            print("%s --> %s" % (title, clean_name(title.split(None, 1)[0])))
            if clean_name(title.split(None, 1)[0]) in wanted:
                positive_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            else:
                negative_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        positive_handle.close()
        negative_handle.close()
    elif pos_file is not None:
        print("Generating matching FASTQ file")
        positive_handle = open(pos_file, "w")
        for title, seq, qual in FastqGeneralIterator(handle):
            if clean_name(title.split(None, 1)[0]) in wanted:
                positive_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        positive_handle.close()
    elif neg_file is not None:
        print("Generating non-matching FASTQ file")
        negative_handle = open(neg_file, "w")
        for title, seq, qual in FastqGeneralIterator(handle):
            if clean_name(title.split(None, 1)[0]) not in wanted:
                negative_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        negative_handle.close()
    handle.close()
    # This does not currently bother to record record counts (faster)


def sff_filter(in_file, pos_file, neg_file, wanted):
    """SFF filter."""
    try:
        from Bio.SeqIO.SffIO import SffIterator, SffWriter
    except ImportError:
        sys.exit("SFF filtering requires Biopython 1.54 or later")

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

    # This makes two passes though the SFF file with isn't so efficient,
    # but this makes the code simple.
    pos_count = neg_count = 0
    if pos_file is not None:
        out_handle = open(pos_file, "wb")
        writer = SffWriter(out_handle, xml=manifest)
        in_handle.seek(0)  # start again after getting manifest
        pos_count = writer.write_file(rec for rec in SffIterator(in_handle) if clean_name(rec.id) in wanted)
        out_handle.close()
    if neg_file is not None:
        out_handle = open(neg_file, "wb")
        writer = SffWriter(out_handle, xml=manifest)
        in_handle.seek(0)  # start again
        neg_count = writer.write_file(rec for rec in SffIterator(in_handle) if clean_name(rec.id) not in wanted)
        out_handle.close()
    # And we're done
    in_handle.close()
    # At the time of writing, Galaxy doesn't show SFF file read counts,
    # so it is useful to put them in stdout and thus shown in job info.
    return pos_count, neg_count


if seq_format.lower() == "sff":
    # Now write filtered SFF file based on IDs wanted
    pos_count, neg_count = sff_filter(in_file, out_positive_file, out_negative_file, ids)
    # At the time of writing, Galaxy doesn't show SFF file read counts,
    # so it is useful to put them in stdout and thus shown in job info.
elif seq_format.lower() == "fasta":
    # Write filtered FASTA file based on IDs from tabular file
    pos_count, neg_count = fasta_filter(in_file, out_positive_file, out_negative_file, ids)
    print("%i with and %i without specified IDs" % (pos_count, neg_count))
elif seq_format.lower().startswith("fastq"):
    # Write filtered FASTQ file based on IDs from tabular file
    fastq_filter(in_file, out_positive_file, out_negative_file, ids)
    # This does not currently track the counts
else:
    sys.exit("Unsupported file type %r" % seq_format)
