#!/usr/bin/env python
"""Sub-sample sequence from a FASTA, FASTQ or SFF file.

This tool is a short Python script which requires Biopython 1.62 or later
for sequence parsing. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2014-2015 by Peter Cock, The James Hutton Institute
(formerly the Scottish Crop Research Institute, SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT license).

Use -v or --version to get the version, -h or --help for help.
"""
import os
import sys
from optparse import OptionParser

# Parse Command Line
usage = """Use as follows:

$ python sample_seqs.py [options]

e.g. Sample 20% of the reads:

$ python sample_seqs.py -i my_seq.fastq -f fastq -p 20.0 -o sample.fastq

This samples uniformly though the file, rather than at random, and therefore
should be reproducible.

If you have interleaved paired reads, use the --interleaved switch. If
instead you have two matched files (one for each pair), run the two
twice with the same sampling options to make to matched smaller files.
"""
parser = OptionParser(usage=usage)
parser.add_option(
    "-i",
    "--input",
    dest="input",
    default=None,
    help="Input sequences filename",
    metavar="FILE",
)
parser.add_option(
    "-f",
    "--format",
    dest="format",
    default=None,
    help="Input sequence format (e.g. fasta, fastq, sff)",
)
parser.add_option(
    "-o",
    "--output",
    dest="output",
    default=None,
    help="Output sampled sequenced filename",
    metavar="FILE",
)
parser.add_option(
    "-p",
    "--percent",
    dest="percent",
    default=None,
    help="Take this percent of the reads",
)
parser.add_option(
    "-n", "--everyn", dest="everyn", default=None, help="Take every N-th read"
)
parser.add_option(
    "-c", "--count", dest="count", default=None, help="Take exactly N reads"
)
parser.add_option(
    "--interleaved",
    dest="interleaved",
    default=False,
    action="store_true",
    help="Input is interleaved reads, preserve the pairings",
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
    print("v0.2.4")
    sys.exit(0)

try:
    from Bio import SeqIO
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.SeqIO.SffIO import SffIterator, SffWriter
except ImportError:
    sys.exit("This script requires Biopython.")

in_file = options.input
out_file = options.output
interleaved = options.interleaved

if not in_file:
    sys.exit("Require an input filename")
if in_file != "/dev/stdin" and not os.path.isfile(in_file):
    sys.exit("Missing input file %r" % in_file)
if not out_file:
    sys.exit("Require an output filename")
if not options.format:
    sys.exit("Require the sequence format")
seq_format = options.format.lower()


def count_fasta(filename):
    count = 0
    with open(filename) as handle:
        for title, seq in SimpleFastaParser(handle):
            count += 1
    return count


def count_fastq(filename):
    count = 0
    with open(filename) as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            count += 1
    return count


def count_sff(filename):
    # If the SFF file has a built in index (which is normal),
    # this will be parsed and is the quicker than scanning
    # the whole file.
    return len(SeqIO.index(filename, "sff"))


def count_sequences(filename, format):
    if format == "sff":
        return count_sff(filename)
    elif format == "fasta":
        return count_fasta(filename)
    elif format.startswith("fastq"):
        return count_fastq(filename)
    else:
        sys.exit("Unsupported file type %r" % format)


if options.percent and options.everyn:
    sys.exit("Cannot combine -p and -n options")
elif options.everyn and options.count:
    sys.exit("Cannot combine -p and -c options")
elif options.percent and options.count:
    sys.exit("Cannot combine -n and -c options")
elif options.everyn:
    try:
        N = int(options.everyn)
    except ValueError:
        sys.exit("Bad -n argument %r" % options.everyn)
    if N < 2:
        sys.exit("Bad -n argument %r" % options.everyn)
    if (N % 10) == 1:
        sys.stderr.write("Sampling every %ist sequence\n" % N)
    elif (N % 10) == 2:
        sys.stderr.write("Sampling every %ind sequence\n" % N)
    elif (N % 10) == 3:
        sys.stderr.write("Sampling every %ird sequence\n" % N)
    else:
        sys.stderr.write("Sampling every %ith sequence\n" % N)

    def sampler(iterator):
        """Sample every Nth sequence."""
        global N
        count = 0
        for record in iterator:
            count += 1
            if count % N == 1:
                yield record


elif options.percent:
    try:
        percent = float(options.percent) / 100.0
    except ValueError:
        sys.exit("Bad -p percent argument %r" % options.percent)
    if not (0.0 <= percent <= 1.0):
        sys.exit("Bad -p percent argument %r" % options.percent)
    sys.stderr.write("Sampling %0.3f%% of sequences\n" % (100.0 * percent))

    def sampler(iterator):
        """Sample given percentage of sequences."""
        global percent
        count = 0
        taken = 0
        for record in iterator:
            count += 1
            if percent * count > taken:
                taken += 1
                yield record


elif options.count:
    try:
        N = int(options.count)
    except ValueError:
        sys.exit("Bad -c count argument %r" % options.count)
    if N < 1:
        sys.exit("Bad -c count argument %r" % options.count)
    total = count_sequences(in_file, seq_format)
    sys.stderr.write("Input file has %i sequences\n" % total)
    if interleaved:
        # Paired
        if total % 2:
            sys.exit(
                "Paired mode, but input file has an odd number of sequences: %i" % total
            )
        elif N > total // 2:
            sys.exit(
                "Requested %i sequence pairs, "
                "but file only has %i pairs (%i sequences)." % (N, total // 2, total)
            )
        total = total // 2
        if N == 1:
            sys.stderr.write("Sampling just first sequence pair!\n")
        elif N == total:
            sys.stderr.write("Taking all the sequence pairs\n")
        else:
            sys.stderr.write("Sampling %i sequence pairs\n" % N)
    else:
        # Not paired
        if total < N:
            sys.exit("Requested %i sequences, but file only has %i." % (N, total))
        if N == 1:
            sys.stderr.write("Sampling just first sequence!\n")
        elif N == total:
            sys.stderr.write("Taking all the sequences\n")
        else:
            sys.stderr.write("Sampling %i sequences\n" % N)
    if N == total:

        def sampler(iterator):
            """No-operation dummy filter, taking everything."""
            global N
            taken = 0
            for record in iterator:
                taken += 1
                yield record
            assert taken == N, "Picked %i, wanted %i" % (taken, N)

    else:

        def sampler(iterator):
            """Sample given number of sequences."""
            # Mimic the percentage sampler, with double check on final count
            global N, total
            # Do we need a floating point fudge factor epsilon?
            # i.e. What if percentage comes out slighty too low, and
            # we could end up missing last few desired sequences?
            percentage = float(N) / float(total)
            # print("DEBUG: Want %i out of %i sequences/pairs, as a percentage %0.2f"
            #      % (N, total, percentage * 100.0))
            count = 0
            taken = 0
            for record in iterator:
                count += 1
                # Do we need the extra upper bound?
                if percentage * count > taken and taken < N:
                    taken += 1
                    yield record
                elif total - count + 1 <= N - taken:
                    # remaining records (incuding this one) <= what we still need.
                    # This is a safey check for floating point edge cases where
                    # we need to take all remaining sequences to meet target
                    taken += 1
                    yield record
            assert taken == N, "Picked %i, wanted %i" % (taken, N)


else:
    sys.exit("Must use either -n, -p or -c")


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
    """Yield raw FASTA records as multi-line strings."""
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break

    no_id_warned = False
    while True:
        if line[0] != ">":
            raise ValueError("Records in Fasta files should start with '>' character")
        try:
            line[1:].split(None, 1)[0]
        except IndexError:
            if not no_id_warned:
                sys.stderr.write("WARNING - Malformed FASTA entry with no identifier\n")
                no_id_warned = True
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
            return  # StopIteration


def fasta_filter(in_file, out_file, iterator_filter, inter):
    count = 0
    # Galaxy now requires Python 2.5+ so can use with statements,
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
                for title, seq, qual in iterator_filter(
                    FastqGeneralIterator(in_handle)
                ):
                    count += 1
                    pos_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    return count


def sff_filter(in_file, out_file, iterator_filter, inter):
    count = 0
    try:
        from Bio.SeqIO.SffIO import ReadRocheXmlManifest
    except ImportError:
        # Prior to Biopython 1.56 this was a private function
        from Bio.SeqIO.SffIO import _sff_read_roche_index_xml as ReadRocheXmlManifest
    with open(in_file, "rb") as in_handle:
        try:
            manifest = ReadRocheXmlManifest(in_handle)
        except ValueError:
            manifest = None
        in_handle.seek(0)
        with open(out_file, "wb") as out_handle:
            writer = SffWriter(out_handle, xml=manifest)
            in_handle.seek(0)  # start again after getting manifest
            if inter:
                from itertools import chain

                count = writer.write_file(
                    chain.from_iterable(iterator_filter(pair(SffIterator(in_handle))))
                )
                assert count % 2 == 0, "Odd number of records? %i" % count
                count /= 2
            else:
                count = writer.write_file(iterator_filter(SffIterator(in_handle)))
    return count


if seq_format == "sff":
    count = sff_filter(in_file, out_file, sampler, interleaved)
elif seq_format == "fasta":
    count = fasta_filter(in_file, out_file, sampler, interleaved)
elif seq_format.startswith("fastq"):
    count = fastq_filter(in_file, out_file, sampler, interleaved)
else:
    sys.exit("Unsupported file type %r" % seq_format)

if interleaved:
    sys.stderr.write("Selected %i pairs\n" % count)
else:
    sys.stderr.write("Selected %i records\n" % count)
