#!/usr/bin/env python
"""Filter a FASTA, FASTQ or SSF file with SAM/BAM mapping information.

When filtering an SFF file, any Roche XML manifest in the input file is
preserved in both output files.

This tool is a short Python script which requires Biopython 1.54 or later
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
import re
import subprocess
from optparse import OptionParser

#Parse Command Line
usage = """Use as follows:

$ python seq_filter_by_mapping.py [options] mapping.sam/bam [more mappings]

e.g. Positive matches using column one from a single BAM file:

$ seq_filter_by_mapping.py -i my_seqs.fastq -f fastq -p matches.fastq mapping.bam

Multiple SAM/BAM mapping files may be given.
"""
parser = OptionParser(usage=usage)
parser.add_option('-i', '--input', dest='input',
                  default=None, help='Input sequences filename',
                  metavar="FILE")
parser.add_option('-f', '--format', dest='format',
                  default=None,
                  help='Input sequence format (e.g. fasta, fastq, sff)')
parser.add_option('-p', '--positive', dest='output_positive',
                  default=None,
                  help='Output filename for mapping reads',
                  metavar="FILE")
parser.add_option('-n', '--negative', dest='output_negative',
                  default=None,
                  help='Output filename for non-mapping reads',
                  metavar="FILE")
parser.add_option("-m", "--pair-mode", dest="pair_mode",
                  default="lax",
                  help="How to treat paired reads (lax or strict, default lax)")
parser.add_option("-v", "--version", dest="version",
                  default=False, action="store_true",
                  help="Show version and quit")

options, args = parser.parse_args()

if options.version:
    print "v0.0.3"
    sys.exit(0)

in_file = options.input
seq_format = options.format
out_positive_file = options.output_positive
out_negative_file = options.output_negative
pair_mode = options.pair_mode

if in_file is None or not os.path.isfile(in_file):
    sys.exit("Missing input file: %r" % in_file)
if out_positive_file is None and out_negative_file is None:
    sys.exit("Neither output file requested")
if seq_format is None:
    sys.exit("Missing sequence format")
if pair_mode not in ["lax", "strict"]:
    sys.exit("Pair mode argument should be 'lax' or 'strict', not %r" % pair_mode)
for mapping in args:
    if not os.path.isfile(mapping):
        sys.exit("Mapping file %r not found" % mapping)
if not args:
    sys.exit("At least one SAM/BAM mapping file is required")


#Cope with three widely used suffix naming convensions,
#Illumina: /1 or /2
#Forward/revered: .f or .r
#Sanger, e.g. .p1k and .q1k
#See http://staden.sourceforge.net/manual/pregap4_unix_50.html
#re_f = re.compile(r"(/1|\.f|\.[sfp]\d\w*)$")
#re_r = re.compile(r"(/2|\.r|\.[rq]\d\w*)$")
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

def clean_name(name):
    """Remove suffix."""
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

mapped_chars = { '>' :'__gt__',
                 '<' :'__lt__',
                 "'" :'__sq__',
                 '"' :'__dq__',
                 '[' :'__ob__',
                 ']' :'__cb__',
                 '{' :'__oc__',
                 '}' :'__cc__',
                 '@' : '__at__',
                 '\n' : '__cn__',
                 '\r' : '__cr__',
                 '\t' : '__tc__',
                 '#' : '__pd__'
                 }

def load_mapping_ids(filename, pair_mode, ids):
    """Parse SAM/BAM file, updating given set of ids.

    Parses BAM files via call out to samtools view command.
    """
    handle = open(filename, "rb")
    magic = handle.read(4)
    if magic == b"\x1f\x8b\x08\x04":
        # Presumably a BAM file then...
        handle.close()
        # Call samtools view, don't need header so no -h added:
        child = subprocess.Popen(["samtools", "view", filename],
                                 stdin=None,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        handle = child.stdout
    else:
        # Presumably a SAM file...
        child = None
        handle.seek(0)
    # Handle should now contain SAM records
    for line in handle:
        # Ignore header lines
        if line[0] != "@":
            qname, flag, rest = line.split("\t", 2)
            flag = int(flag)
            if pair_mode == "lax":
                # If either read or its partner is mapped, take it!
                # Being lazy, since we will look at both reads
                # can just check if (either) has 0x4 clear.
                if not (flag & 0x4):
                    ids.add(qname)
            elif pair_mode == "strict":
                # For paired reads, require BOTH be mapped.
                if (flag & 0x4):
                    # This is unmapped, ignore it
                    pass
                elif not (flag & 0x1):
                    # Unpaired (& mapped) - take it
                    ids.add(qname)
                elif not (flag & 0x8):
                    # Paired and partner also mapped, good
                    ids.add(qname)
    if child:
        # Check terminated normally.
        stdout, stderr = child.communicate()
        assert child.returncode is not None
        if child.returncode:
            msg = "Error %i from 'samtools view %s'\n%s" % (child.returncode,
                                                            filename, stderr)
            sys.exit(msg.strip(), child.returncode)
    else:
        handle.close()


# Read mapping file(s) and record all mapped identifiers
ids = set()
for filename in args:
    load_mapping_ids(filename, pair_mode, ids)
# TODO - If want to support naive paired mode, have to record
# more than just qname (need /1 or /2 indicator)
print("Loaded %i mapped IDs" % (len(ids)))
if len(ids) < 10:
    print("Looking for %s" % ", ".join(sorted(ids)))

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
        if pos_file is not None and neg_file is not None:
            print "Generating two FASTA files"
            with open(pos_file, "w") as pos_handle:
                with open(neg_file, "w") as neg_handle:
                    for identifier, record in crude_fasta_iterator(in_handle):
                        if clean_name(identifier) in wanted:
                            pos_handle.write(record)
                            pos_count += 1
                        else:
                            neg_handle.write(record)
                            neg_count += 1
        elif pos_file is not None:
            print "Generating matching FASTA file"
            with open(pos_file, "w") as pos_handle:
                for identifier, record in crude_fasta_iterator(in_handle):
                    if clean_name(identifier) in wanted:
                        pos_handle.write(record)
                        pos_count += 1
                    else:
                        neg_count += 1
        else:
            print "Generating non-matching FASTA file"
            assert neg_file is not None
            with open(neg_file, "w") as neg_handle:
                for identifier, record in crude_fasta_iterator(in_handle):
                    if clean_name(identifier) in wanted:
                        pos_count += 1
                    else:
                        neg_handle.write(record)
                        neg_count += 1
    return pos_count, neg_count


def fastq_filter(in_file, pos_file, neg_file, wanted):
    """FASTQ filter."""
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    pos_count = neg_count = 0
    handle = open(in_file, "r")
    if out_positive_file is not None and out_negative_file is not None:
        print "Generating two FASTQ files"
        positive_handle = open(out_positive_file, "w")
        negative_handle = open(out_negative_file, "w")
        print in_file
        for title, seq, qual in FastqGeneralIterator(handle):
            # print("%s --> %s" % (title, clean_name(title.split(None, 1)[0])))
            if clean_name(title.split(None, 1)[0]) in ids:
                positive_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                pos_count += 1
            else:
                negative_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                neg_count += 1
        positive_handle.close()
        negative_handle.close()
    elif out_positive_file is not None:
        print "Generating matching FASTQ file"
        positive_handle = open(out_positive_file, "w")
        for title, seq, qual in FastqGeneralIterator(handle):
            if clean_name(title.split(None, 1)[0]) in ids:
                positive_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                pos_count += 1
            else:
                neg_count += 1
        positive_handle.close()
    elif out_negative_file is not None:
        print "Generating non-matching FASTQ file"
        negative_handle = open(out_negative_file, "w")
        for title, seq, qual in FastqGeneralIterator(handle):
            if clean_name(title.split(None, 1)[0]) in ids:
                pos_count += 1
            else:
                negative_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                neg_count += 1
        negative_handle.close()
    handle.close()
    return pos_count, neg_count


def sff_filter(in_file, pos_file, neg_file, wanted):
    """SFF filter."""
    try:
        from Bio.SeqIO.SffIO import SffIterator, SffWriter
    except ImportError:
        sys.exit("SFF filtering requires Biopython 1.54 or later")

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
    if out_positive_file is not None:
        out_handle = open(out_positive_file, "wb")
        writer = SffWriter(out_handle, xml=manifest)
        in_handle.seek(0) #start again after getting manifest
        pos_count = writer.write_file(rec for rec in SffIterator(in_handle) if clean_name(rec.id) in ids)
        out_handle.close()
    if out_negative_file is not None:
        out_handle = open(out_negative_file, "wb")
        writer = SffWriter(out_handle, xml=manifest)
        in_handle.seek(0) #start again
        neg_count = writer.write_file(rec for rec in SffIterator(in_handle) if clean_name(rec.id) not in ids)
        out_handle.close()
    #And we're done
    in_handle.close()
    return pos_count, neg_count


if seq_format.lower()=="sff":
    sequence_filter = sff_filter
elif seq_format.lower()=="fasta":
    sequence_filter = fasta_filter
elif seq_format.lower().startswith("fastq"):
    sequence_filter = fastq_filter
else:
    sys.exit("Unsupported file type %r" % seq_format)

pos_count, neg_count = sequence_filter(in_file, out_positive_file, out_negative_file, ids)
print("%i mapped and %i unmapped reads." % (pos_count, neg_count))
fraction = float(pos_count) * 100.0 / float(pos_count + neg_count)
print("In total %i reads, of which %0.1f%% mapped." % (pos_count + neg_count, fraction))
