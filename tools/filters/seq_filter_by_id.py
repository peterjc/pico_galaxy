#!/usr/bin/env python
"""Filter a FASTA, FASTQ or SSF file with IDs from a tabular file.

Takes six command line options, tabular filename, ID column numbers (comma
separated list using one based counting), input filename, input type (e.g.
FASTA or SFF) and two output filenames (for records with and without the
given IDs, same format as input sequence file).

If either output filename is just a minus sign, that file is not created.
This is intended to allow output for just the matched (or just the non-matched)
records.

When filtering an SFF file, any Roche XML manifest in the input file is
preserved in both output files.

Note in the default NCBI BLAST+ tabular output, the query sequence ID is
in column one, and the ID of the match from the database is in column two.
Here sensible values for the column numbers would therefore be "1" or "2".

This tool is a short Python script which requires Biopython 1.54 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2010-2013 by Peter Cock, The James Hutton Institute
(formerly the Scottish Crop Research Institute, SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.6 of the script, use -v or --version to get the version.
"""
import sys

def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.6"
    sys.exit(0)

#Parse Command Line
try:
    tabular_file, cols_arg, in_file, seq_format, out_positive_file, out_negative_file = sys.argv[1:]
except ValueError:
    stop_err("Expected six arguments (tab file, columns, in seq, seq format, out pos, out neg), got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))
try:
    columns = [int(arg)-1 for arg in cols_arg.split(",")]
except ValueError:
    stop_err("Expected list of columns (comma separated integers), got %s" % cols_arg)
if min(columns) < 0:
    stop_err("Expect one-based column numbers (not zero-based counting), got %s" % cols_arg)

if out_positive_file == "-" and out_negative_file == "-":
    stop_err("Neither output file requested")


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


def crude_fasta_iterator(handle):
    """Yields tuples, record ID and the full record as a string."""
    while True:
        line = handle.readline()
        if line == "":
            return # Premature end of file, or just empty?
        if line[0] == ">":
            break

    no_id_warned = False
    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in Fasta files should start with '>' character")
        try:
            id = line[1:].split(None, 1)[0]
        except IndexError:
            if not no_id_warned:
                sys.stderr.write("WARNING - Malformed FASTA entry with no identifier\n")
                no_id_warned = True
            id = None
        lines = [line]
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line)
            line = handle.readline()
        yield id, "".join(lines)
        if not line:
            return # StopIteration


def fasta_filter(in_file, pos_file, neg_file, wanted):
    """FASTA filter producing 60 character line wrapped outout."""
    pos_count = neg_count = 0
    #Galaxy now requires Python 2.5+ so can use with statements,
    with open(in_file) as in_handle:
        #Doing the if statement outside the loop for speed
        #(with the downside of three very similar loops).
        if pos_file != "-" and neg_file != "-":
            print "Generating two FASTA files"
            with open(pos_file, "w") as pos_handle:
                with open(neg_file, "w") as neg_handle:
                    for identifier, record in crude_fasta_iterator(in_handle):
                        if identifier in wanted:
                            pos_handle.write(record)
                            pos_count += 1
                        else:
                            neg_handle.write(record)
                            neg_count += 1
        elif pos_file != "-":
            print "Generating matching FASTA file"
            with open(pos_file, "w") as pos_handle:
                for identifier, record in crude_fasta_iterator(in_handle):
                    if identifier in wanted:
                        pos_handle.write(record)
                        pos_count += 1
                    else:
                        neg_count += 1
        else:
            print "Generating non-matching FASTA file"
            assert neg_file != "-"
            with open(neg_file, "w") as neg_handle:
                for identifier, record in crude_fasta_iterator(in_handle):
                    if identifier in wanted:
                        pos_count += 1
                    else:
                        neg_handle.write(record)
                        neg_count += 1
    return pos_count, neg_count


if seq_format.lower()=="sff":
    #Now write filtered SFF file based on IDs from BLAST file
    try:
        from Bio.SeqIO.SffIO import SffIterator, SffWriter
    except ImportError:
        stop_err("SFF filtering requires Biopython 1.54 or later")

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
    #This makes two passes though the SFF file with isn't so efficient,
    #but this makes the code simple.
    pos_count = neg_count = 0
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
    #At the time of writing, Galaxy doesn't show SFF file read counts,
    #so it is useful to put them in stdout and thus shown in job info.
    print "%i with and %i without specified IDs" % (pos_count, neg_count)
elif seq_format.lower()=="fasta":
    #Write filtered FASTA file based on IDs from tabular file
    pos_count, neg_count = fasta_filter(in_file, out_positive_file, out_negative_file, ids)
    print "%i with and %i without specified IDs" % (pos_count, neg_count)
elif seq_format.lower().startswith("fastq"):
    #Write filtered FASTQ file based on IDs from tabular file
    from galaxy_utils.sequence.fastq import fastqReader, fastqWriter
    reader = fastqReader(open(in_file, "rU"))
    if out_positive_file != "-" and out_negative_file != "-":
        print "Generating two FASTQ files"
        positive_writer = fastqWriter(open(out_positive_file, "w"))
        negative_writer = fastqWriter(open(out_negative_file, "w"))
        for record in reader:
            #The [1:] is because the fastaReader leaves the > on the identifier.
            if record.identifier and record.identifier.split()[0][1:] in ids:
                positive_writer.write(record)
            else:
                negative_writer.write(record)
        positive_writer.close()
        negative_writer.close()
    elif out_positive_file != "-":
        print "Generating matching FASTQ file"
        positive_writer = fastqWriter(open(out_positive_file, "w"))
        for record in reader:
            #The [1:] is because the fastaReader leaves the > on the identifier.
            if record.identifier and record.identifier.split()[0][1:] in ids:
                positive_writer.write(record)
        positive_writer.close()
    elif out_negative_file != "-":
        print "Generating non-matching FASTQ file"
        negative_writer = fastqWriter(open(out_negative_file, "w"))
        for record in reader:
            #The [1:] is because the fastaReader leaves the > on the identifier.
            if not record.identifier or record.identifier.split()[0][1:] not in ids:
                negative_writer.write(record)
        negative_writer.close()
    reader.close()
else:
    stop_err("Unsupported file type %r" % seq_format)
