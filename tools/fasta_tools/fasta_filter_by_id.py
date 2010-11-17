#!/usr/bin/env python
"""Filter a FASTA file with IDs from a tabular file, e.g. from BLAST.

Takes five command line options, tabular BLAST filename, ID column numbers
(comma separated list using one based counting), input FASTA filename, and
two output FASTA filenames (for records with and without any BLAST hits).

Note in the default NCBI BLAST+ tabular output, the query sequence ID is
in column one, and the ID of the match from the database is in column two.
"""
import sys
from galaxy_utils.sequence.fasta import fastaReader, fastaWriter

#Parse Command Line
tabular_file, cols_arg, in_file, out_positive_file, out_negative_file = sys.argv[1:]
columns = [int(arg)-1 for arg in cols_arg.split(",")]

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

#Write filtered FASTA file based on IDs from BLAST file
reader = fastaReader(open(in_file, "rU"))
positive_writer = fastaWriter(open(out_positive_file, "w"))
negative_writer = fastaWriter(open(out_negative_file, "w"))
for record in reader:
    #The [1:] is because the fastaReader leaves the > on the identifer.
    if record.identifier and record.identifier.split()[0][1:] in ids:
        positive_writer.write(record)
    else:
        negative_writer.write(record)
positive_writer.close()
negative_writer.close()
reader.close()
