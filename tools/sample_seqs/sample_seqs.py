#!/usr/bin/env python
"""Sub-sample sequence from a FASTA, FASTQ or SFF file.

This tool is a short Python script which requires Biopython 1.62 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2010-2013 by Peter Cock, The James Hutton Institute
(formerly the Scottish Crop Research Institute, SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT license).

This is version 0.1.0 of the script, use -v or --version to get the version.
"""
import os
import sys

def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.1.0"
    sys.exit(0)

#Parse Command Line
if len(sys.argv) < 5:
    stop_err("Requires at least four arguments: seq_format, in_file, out_file, mode, ...")
seq_format, in_file, out_file, mode = sys.argv[1:5]
if mode == "everyNth":
    if len(sys.argv) != 6:
        stop_err("If using everyNth, just need argument N")
    try:
        N = int(sys.argv[5])
    except:
        stop_err("Bad N argument %r" % sys.argv[5])
    if N < 2:
        stop_err("Bad N argument %r" % sys.argv[5])
    if (N % 10) == 1:
        print("Sampling every %ist sequence" % N)
    elif (N % 10) == 2:
        print("Sampling every %ind sequence" % N)
    elif (N % 10) == 3:
        print("Sampling every %ird sequence" % N)
    else:
        print("Sampling every %ith sequence" % N)
else:
    stop_err("Unsupported mode %r" % mode)
if not os.path.isfile(in_file):
    stop_err("Missing input file %r" % in_file)


def pick_every_N(iterator, N):
    count = 0
    for record in iterator:
        count += 1
        if count % N == 1:
            yield record

def raw_fasta_iterator(handle):
    """Yields raw FASTA records as multi-line strings."""
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
        yield "".join(lines)
        if not line:
            return # StopIteration 

def fasta_filter_every_N(in_file, out_file, N):
    count = 0
    #Galaxy now requires Python 2.5+ so can use with statements,
    with open(in_file) as in_handle:
        with open(out_file, "w") as pos_handle:
            for record in pick_every_N(raw_fasta_iterator(in_handle), N):
                count += 1
                pos_handle.write(record)
    return count

try:
    from galaxy_utils.sequence.fastq import fastqReader, fastqWriter
    def fastq_filter_every_N(in_file, out_file, N):
        count = 0
        #from galaxy_utils.sequence.fastq import fastqReader, fastqWriter
        reader = fastqReader(open(in_file, "rU"))
        writer = fastqWriter(open(out_file, "w"))
        for record in pick_every_N(reader, N):
            count += 1
            writer.write(record)
        writer.close()
        reader.close()
        return count
except ImportError:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    def fastq_filter_every_N(in_file, out_file, N):
        count = 0
        with open(in_file) as in_handle:
            with open(out_file, "w") as pos_handle:
                for title, seq, qual in pick_every_N(FastqGeneralIterator(in_handle), N):
                    count += 1
                    pos_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        return count

def sff_filter_every_N(in_file, out_file, N):
    count = 0
    try:
        from Bio.SeqIO.SffIO import SffIterator, SffWriter
    except ImportError:
        stop_err("SFF filtering requires Biopython 1.54 or later")
    try:
        from Bio.SeqIO.SffIO import ReadRocheXmlManifest
    except ImportError:
        #Prior to Biopython 1.56 this was a private function
        from Bio.SeqIO.SffIO import _sff_read_roche_index_xml as ReadRocheXmlManifest
    with open(in_file, "rb") as in_handle:
        try:
            manifest = ReadRocheXmlManifest(in_handle)
        except ValueError:
            manifest = None
        in_handle.seek(0)
        with open(out_file, "wb") as out_handle:
            writer = SffWriter(out_handle, xml=manifest)
            in_handle.seek(0) #start again after getting manifest
            count = writer.write_file(pick_every_N(SffIterator(in_handle), N))
            #count = writer.write_file(SffIterator(in_handle))
    return count

if seq_format.lower()=="sff":
    count = sff_filter_every_N(in_file, out_file, N)
elif seq_format.lower()=="fasta":
    count = fasta_filter_every_N(in_file, out_file, N)
elif seq_format.lower().startswith("fastq"):
    count = fastq_filter_every_N(in_file, out_file, N)
else:
    stop_err("Unsupported file type %r" % seq_format)

print("Sampled %i records" % count)
