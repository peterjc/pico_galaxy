#!/usr/bin/env python
"""Rename FASTA, QUAL, FASTQ or SSF sequences with ID mapping from tabular file.

Takes six command line options, tabular filename, current (old)  ID column
number (using one based counting), new ID column number (also using one based
counting), input sequence filename, input type (e.g. FASTA or SFF) and the
output filename (same format as input sequence file).

When selecting from an SFF file, any Roche XML manifest in the input file is
preserved in both output files.

This tool is a short Python script which requires Biopython 1.54 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2011-2013 by Peter Cock, The James Hutton Institute UK.
All rights reserved. See accompanying text file for licence details (MIT
license).

This is version 0.0.4 of the script.
"""
import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.4"
    sys.exit(0)

def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

#Parse Command Line
try:
    tabular_file, old_col_arg, new_col_arg, in_file, seq_format, out_file = sys.argv[1:]
except ValueError:
    stop_err("Expected six arguments (tabular file, old col, new col, input file, format, output file), got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))

try:
    if old_col_arg.startswith("c"):
        old_column = int(old_col_arg[1:])-1
    else:
        old_column = int(old_col_arg)-1
except ValueError:
    stop_err("Expected column number, got %s" % old_col_arg)
try:
    if old_col_arg.startswith("c"):
        new_column = int(new_col_arg[1:])-1
    else:
        new_column = int(new_col_arg)-1
except ValueError:
    stop_err("Expected column number, got %s" % new_col_arg)
if old_column == new_column:
    stop_err("Old and new column arguments are the same!")

def parse_ids(tabular_file, old_col, new_col):
    """Read tabular file and record all specified ID mappings."""
    handle = open(tabular_file, "rU")
    for line in handle:
        if not line.startswith("#"):
            parts = line.rstrip("\n").split("\t")
            yield parts[old_col].strip(), parts[new_col].strip()
    handle.close()

#Load the rename mappings
rename = dict(parse_ids(tabular_file, old_column, new_column))
print "Loaded %i ID mappings" % len(rename)
              
#Rewrite the sequence file
if seq_format.lower()=="sff":
    #Use Biopython for this format
    renamed = 0
    def rename_seqrecords(records, mapping):
        global renamed #nasty, but practical!
        for record in records:
            try:
                record.id = mapping[record.id]
                renamed += 1
            except KeyError:
                pass
            yield record
                                                                
    try:
        from Bio.SeqIO.SffIO import SffIterator, SffWriter
    except ImportError:
        stop_err("Requires Biopython 1.54 or later")

    try:
        from Bio.SeqIO.SffIO import ReadRocheXmlManifest
    except ImportError:
        #Prior to Biopython 1.56 this was a private function
        from Bio.SeqIO.SffIO import _sff_read_roche_index_xml as ReadRocheXmlManifest

    in_handle = open(in_file, "rb") #must be binary mode!
    try:
        manifest = ReadRocheXmlManifest(in_handle)
    except ValueError:
        manifest = None
    out_handle = open(out_file, "wb")
    writer = SffWriter(out_handle, xml=manifest)
    in_handle.seek(0) #start again after getting manifest
    count = writer.write_file(rename_seqrecords(SffIterator(in_handle), rename))
    out_handle.close()
    in_handle.close()
else:
    #Use Galaxy for FASTA, QUAL or FASTQ
    if seq_format.lower() in ["fasta", "csfasta"] \
    or seq_format.lower().startswith("qual"):
        from galaxy_utils.sequence.fasta import fastaReader, fastaWriter
        reader = fastaReader(open(in_file, "rU"))
        writer = fastaWriter(open(out_file, "w"))
        marker = ">"
    elif seq_format.lower().startswith("fastq"):
        from galaxy_utils.sequence.fastq import fastqReader, fastqWriter
        reader = fastqReader(open(in_file, "rU"))
        writer = fastqWriter(open(out_file, "w"))
        marker = "@"
    else:
        stop_err("Unsupported file type %r" % seq_format)
    #Now do the renaming
    count = 0
    renamed = 0
    for record in reader:
        #The [1:] is because the fastaReader leaves the > on the identifier,
        #likewise the fastqReader leaves the @ on the identifier
        try:
            idn, descr = record.identifier[1:].split(None, 1)
        except ValueError:
            idn = record.identifier[1:]
            descr = None
        if idn in rename:
            if descr:
                record.identifier = "%s%s %s" % (marker, rename[idn], descr)
            else:
                record.identifier = "%s%s" % (marker, rename[idn])
            renamed += 1
        writer.write(record)
        count += 1
    writer.close()
    reader.close()

print "Renamed %i out of %i records" % (renamed, count)
