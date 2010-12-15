#!/usr/bin/env python
"""Filter an SSF file with IDs from a tabular file, e.g. from BLAST.

Takes five command line options, tabular filename, ID column numbers
(comma separated list using one based counting), input SFF filename, and
two output SFF filenames (for records with and without the given IDs).

If either output filename is just a minus sign, that file is not created.
This is intended to allow output for just the matched (or just the non-matched)
records.

Any Roche XML manifest in the input file is preserved in both output files.

Note in the default NCBI BLAST+ tabular output, the query sequence ID is
in column one, and the ID of the match from the database is in column two.
Here sensible values for the column numbers would therefore be "1" or "2".

This tool is a short Python script which requires Biopython 1.54 or later.
If you use this tool in scientific work leading to a publication, please cite
the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2010 by Peter Cock, SCRI, UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script.
"""
import sys

def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

try:
    from Bio.SeqIO.SffIO import SffIterator, SffWriter
except ImportError:
    stop_err("Requires Biopython 1.54 or later")

try:
    from Bio.SeqIO.SffIO import ReadRocheXmlManifest
except ImportError:
    #Prior to Biopython 1.56 this was a private function
    from Bio.SeqIO.SffIO import _sff_read_roche_index_xml as ReadRocheXmlManifest

#Parse Command Line
try:
    tabular_file, cols_arg, in_file, out_positive_file, out_negative_file = sys.argv[1:]
except ValueError:
    stop_err("Expected five arguments, got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))
try:
    columns = [int(arg)-1 for arg in cols_arg.split(",")]
except ValueError:
    stop_err("Expected list of columns (comma separated integers), got %s" % cols_arg)

#Read tabular file and record all specified identifiers
ids = set()
handle = open(tabular_file, "rU")
if len(columns)>1:
    #General case of many columns
    for line in handle:
        if line.startswith("#"):
            #Ignore comments
            continue
        parts = line.rstrip("\n").split("\t")
        for col in columns:
            ids.add(parts[col])
    print "Using %i IDs from %i columns of tabular file" % (len(ids), len(columns))
else:
    #Single column, special case speed up
    col = columns[0]
    for line in handle:
        if not line.startswith("#"):
            ids.add(line.rstrip("\n").split("\t")[col])
    print "Using %i IDs from tabular file" % (len(ids))
handle.close()

#Now write filtered SFF file based on IDs from BLAST file
in_handle = open(in_file, "rb") #must be binary mode!
try:
    manifest = ReadRocheXmlManifest(in_handle)
except ValueError:
    manifest = None

#This makes two passes though the SFF file with isn't so efficient,
#but this makes the code simple.

if out_positive_file != "-":
    out_handle = open(out_positive_file, "wb")
    writer = SffWriter(out_handle, xml=manifest)
    in_handle.seek(0) #start again after getting manifest
    pos_count = writer.write_file(rec for rec in SffIterator(in_handle) if rec.id in ids)
    out_handle.close()

if out_negative_file != "-":
    out_handle = open(out_negative_file, "wb")
    writer = SffWriter(out_handle, xml=manifest)
    in_handle.seek(0) #start again
    neg_count = writer.write_file(rec for rec in SffIterator(in_handle) if rec.id not in ids)
    out_handle.close()

#And we're done
in_handle.close()

if out_positive_file != "-" and out_negative_file != "-":
    print "%i with and %i without specified IDs" % (pos_count, neg_count)
elif out_positive_file != "-":
    print "%i with specified IDs" % pos_count
elif out_negative_file != "-":
    print "%i without specified IDs" % neg_count
else:
    stop_err("Neither output file requested")
