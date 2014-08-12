#!/usr/bin/env python
"""Sub-sample sequence from a FASTA, FASTQ or SFF file.

This tool is a short Python script which requires Biopython 1.62 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2014 by Peter Cock, The James Hutton Institute
(formerly the Scottish Crop Research Institute, SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT license).

Use -v or --version to get the version, -h or --help for help.
"""
import os
import sys
from optparse import OptionParser


def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

#Parse Command Line
usage = """Use as follows:

$ python sample_seqs.py [options]

e.g. Sample 20% of the reads:

$ python sample_seqs.py -i my_seq.fastq -f fastq -p 20.0 -o sample.fastq
"""
parser = OptionParser(usage=usage)
parser.add_option('-i', '--input', dest='input',
                  default=None, help='Input sequences filename',
                  metavar="FILE")
parser.add_option('-f', '--format', dest='format',
                  default=None,
                  help='Input sequence format (e.g. fasta, fastq, sff)')
parser.add_option('-o', '--output', dest='output',
                  default=None, help='Output sampled sequenced filename',
                  metavar="FILE")
parser.add_option('-p', '--percent', dest='percent',
                  default=None,
                  help='Take this percent of the reads')
parser.add_option('-n', '--everyn', dest='everyn',
                  default=None,
                  help='Take every N-th read')
parser.add_option("--interleaved", dest="interleaved",
                  default=False, action="store_true",
                  help="Input is interleaved reads, preserve the pairings")
parser.add_option("-v", "--version", dest="version",
                  default=False, action="store_true",
                  help="Show version and quit")
options, args = parser.parse_args()

if options.version:
    print("v0.1.2")
    sys.exit(0)

seq_format = options.format
in_file = options.input
out_file = options.output
interleaved = options.interleaved

if not in_file:
    stop_err("Require an input filename")
if in_file != "/dev/stdin" and not os.path.isfile(in_file):
    stop_err("Missing input file %r" % in_file)
if not out_file:
    stop_err("Require and output filename")


if options.percent and options.everyn:
    stop_err("Cannot combine -p and -n options")
elif options.everyn:
    try:
        N = int(options.everyn)
    except:
        stop_err("Bad N argument %r" % options.everyn)
    if N < 2:
        stop_err("Bad N argument %r" % options.everyn)
    if (N % 10) == 1:
        sys.stderr.write("Sampling every %ist sequence\n" % N)
    elif (N % 10) == 2:
        sys.stderr.write("Sampling every %ind sequence\n" % N)
    elif (N % 10) == 3:
        sys.stderr.write("Sampling every %ird sequence\n" % N)
    else:
        sys.stderr.write("Sampling every %ith sequence\n" % N)
    def sampler(iterator):
        global N
        count = 0
        for record in iterator:
            count += 1
            if count % N == 1:
                yield record
elif options.percent:
    try:
        percent = float(options.percent) / 100.0
    except:
        stop_err("Bad percent argument %r" % options.percent)
    if percent <= 0.0 or 1.0 <= percent:
        stop_err("Bad percent argument %r" % options.percent)
    sys.stderr.write("Sampling %0.3f%% of sequences\n" % (100.0 * percent))
    def sampler(iterator):
        global percent
        count = 0
        taken = 0
        for record in iterator:
            count += 1
            if percent * count > taken:
                taken += 1
                yield record
else:
    stop_err("Must use either -n or -p")


def pair(iterator):
    """Quick and dirty pair batched iterator."""
    while True:
        a = next(iterator)
        b = next(iterator)
        if not b:
            assert not a, "Odd number of records?"
            break
        yield (a, b)


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

def fasta_filter(in_file, out_file, iterator_filter, inter):
    count = 0
    #Galaxy now requires Python 2.5+ so can use with statements,
    with open(in_file) as in_handle:
        with open(out_file, "w") as pos_handle:
            if inter:
                for r1, r2 in iterator_filter(pair(raw_fasta_iterator(in_handle))):
                    count += 1
                    pos_handle.write(r1)
                    pos_handle.write(r2)
            else:
                for record in iterator_filter(raw_fasta_iterator(in_handle)):
                    count += 1
                    pos_handle.write(record)
    return count

try:
    from galaxy_utils.sequence.fastq import fastqReader, fastqWriter
    def fastq_filter(in_file, out_file, iterator_filter, inter):
        count = 0
        #from galaxy_utils.sequence.fastq import fastqReader, fastqWriter
        reader = fastqReader(open(in_file, "rU"))
        writer = fastqWriter(open(out_file, "w"))
        if inter:
            for r1, r2 in iterator_filter(pair(reader)):
                count += 1
                writer.write(r1)
                writer.write(r2)
        else:
            for record in iterator_filter(reader):
                count += 1
                writer.write(record)
        writer.close()
        reader.close()
        return count
except ImportError:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    def fastq_filter(in_file, out_file, iterator_filter, inter):
        count = 0
        with open(in_file) as in_handle:
            with open(out_file, "w") as pos_handle:
                if inter:
                    for r1, r2 in iterator_filter(pair(FastqGeneralIterator(in_handle))):
                        count += 1
                        pos_handle.write("@%s\n%s\n+\n%s\n" % r1)
                        pos_handle.write("@%s\n%s\n+\n%s\n" % r2)
                else:
                    for title, seq, qual in iterator_filter(FastqGeneralIterator(in_handle)):
                        count += 1
                        pos_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        return count

def sff_filter(in_file, out_file, iterator_filter, inter):
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
            if inter:
                from itertools import chain
                count = writer.write_file(chain.from_iterable(iterator_filter(pair(SffIterator(in_handle)))))
                assert count % 2 == 0, "Odd number of records? %i" % count
                count /= 2
            else:
                count = writer.write_file(iterator_filter(SffIterator(in_handle)))
                #count = writer.write_file(SffIterator(in_handle))
    return count

if seq_format.lower()=="sff":
    count = sff_filter(in_file, out_file, sampler, interleaved)
elif seq_format.lower()=="fasta":
    count = fasta_filter(in_file, out_file, sampler, interleaved)
elif seq_format.lower().startswith("fastq"):
    count = fastq_filter(in_file, out_file, sampler, interleaved)
else:
    stop_err("Unsupported file type %r" % seq_format)

if interleaved:
    sys.stderr.write("Selected %i pairs\n" % count)
else:
    sys.stderr.write("Selected %i records\n" % count)
