#!/usr/bin/env python
"""Find ORFs in a nucleotide sequence file.

For more details, see the help text and argument descriptions in the
accompanying get_orfs_or_cdss.xml file which defines a Galaxy interface.

This tool is a short Python script which requires Biopython. If you use
this tool in scientific work leading to a publication, please cite the
Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2011-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI), Dundee, UK. All rights reserved.

See accompanying text file for licence details (MIT licence).
"""

from __future__ import print_function

import re
import sys

from optparse import OptionParser

usage = r"""Use as follows:

$ python get_orfs_or_cdss.py -i genome.fa -f fasta --table 11 \
-t CDS -e open -m all -s both --on cds.nuc.fa --op cds.protein.fa \
--ob cds.bed --og cds.gff3
"""

try:
    from Bio.Seq import Seq, reverse_complement, translate
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    from Bio.Data import CodonTable
except ImportError:
    sys.exit("Missing Biopython library")


parser = OptionParser(usage=usage)
parser.add_option('-i', '--input', dest='input_file',
                  default=None, help='Input fasta file',
                  metavar='FILE')
parser.add_option('-f', '--format', dest='seq_format',
                  default='fasta', help='Sequence format (e.g. fasta, fastq, sff)')
parser.add_option('--table', dest='table',
                  default=1, help='NCBI Translation table', type='int')
parser.add_option('-t', '--ftype', dest='ftype', type='choice',
                  choices=['CDS', 'ORF'], default='ORF',
                  help='Find ORF or CDSs')
parser.add_option('-e', '--ends', dest='ends', type='choice',
                  choices=['open', 'closed'], default='closed',
                  help='Open or closed. Closed ensures start/stop codons are present')
parser.add_option('-m', '--mode', dest='mode', type='choice',
                  choices=['all', 'top', 'one'], default='all',
                  help='Output all ORFs/CDSs from sequence, all ORFs/CDSs '
                  'with max length, or first with maximum length')
parser.add_option('--min_len', dest='min_len',
                  default=10, help='Minimum ORF/CDS length', type='int')
parser.add_option('-s', '--strand', dest='strand', type='choice',
                  choices=['forward', 'reverse', 'both'], default='both',
                  help='Strand to search for features on')
parser.add_option('--on', dest='out_nuc_file',
                  default=None, help='Output nucleotide sequences, or - for STDOUT',
                  metavar='FILE')
parser.add_option('--op', dest='out_prot_file',
                  default=None, help='Output protein sequences, or - for STDOUT',
                  metavar='FILE')
parser.add_option('--ob', dest='out_bed_file',
                  default=None, help='Output BED file, or - for STDOUT',
                  metavar='FILE')
parser.add_option('--og', dest='out_gff3_file',
                  default=None, help='Output GFF3 file, or - for STDOUT',
                  metavar='FILE')
parser.add_option('-v', '--version', dest='version',
                  default=False, action='store_true',
                  help='Show version and quit')

options, args = parser.parse_args()

if options.version:
    print("v0.2.3")
    sys.exit(0)

if not options.input_file:
    sys.exit("Input file is required")

if not any((options.out_nuc_file, options.out_prot_file, options.out_bed_file, options.out_gff3_file)):
    sys.exit("At least one output file is required")

try:
    table_obj = CodonTable.ambiguous_generic_by_id[options.table]
except KeyError:
    sys.exit("Unknown codon table %i" % options.table)

if options.seq_format.lower() == "sff":
    seq_format = "sff-trim"
elif options.seq_format.lower() == "fasta":
    seq_format = "fasta"
elif options.seq_format.lower().startswith("fastq"):
    seq_format = "fastq"
else:
    sys.exit("Unsupported file type %r" % options.seq_format)

print("Genetic code table %i" % options.table)
print("Minimum length %i aa" % options.min_len)
# print "Taking %s ORF(s) from %s strand(s)" % (mode, strand)

starts = sorted(table_obj.start_codons)
assert "NNN" not in starts
re_starts = re.compile("|".join(starts))

stops = sorted(table_obj.stop_codons)
assert "NNN" not in stops
re_stops = re.compile("|".join(stops))


def start_chop_and_trans(s, strict=True):
    """Returns offset, trimmed nuc, protein."""
    if strict:
        assert s[-3:] in stops, s
    assert len(s) % 3 == 0
    for match in re_starts.finditer(s):
        # Must check the start is in frame
        start = match.start()
        if start % 3 == 0:
            n = s[start:]
            assert len(n) % 3 == 0, "%s is len %i" % (n, len(n))
            if strict:
                t = translate(n, options.table, cds=True)
            else:
                # Use when missing stop codon,
                t = "M" + translate(n[3:], options.table, to_stop=True)
            return start, n, t
    return None, None, None


def break_up_frame(s):
    """Returns offset, nuc, protein."""
    start = 0
    for match in re_stops.finditer(s):
        index = match.start() + 3
        if index % 3 != 0:
            continue
        n = s[start:index]
        if options.ftype == "CDS":
            offset, n, t = start_chop_and_trans(n)
        else:
            offset = 0
            t = translate(n, options.table, to_stop=True)
        if n and len(t) >= options.min_len:
            yield start + offset, n, t
        start = index
    if options.ends == "open":
        # No stop codon, Biopython's strict CDS translate will fail
        n = s[start:]
        # Ensure we have whole codons
        # TODO - Try appending N instead?
        # TODO - Do the next four lines more elegantly
        if len(n) % 3:
            n = n[:-1]
        if len(n) % 3:
            n = n[:-1]
        if options.ftype == "CDS":
            offset, n, t = start_chop_and_trans(n, strict=False)
        else:
            offset = 0
            t = translate(n, options.table, to_stop=True)
        if n and len(t) >= options.min_len:
            yield start + offset, n, t


def get_all_peptides(nuc_seq):
    """Returns start, end, strand, nucleotides, protein.

    Co-ordinates are Python style zero-based.
    """
    # TODO - Refactor to use a generator function (in start order)
    # rather than making a list and sorting?
    answer = []
    full_len = len(nuc_seq)
    if options.strand != "reverse":
        for frame in range(0, 3):
            for offset, n, t in break_up_frame(nuc_seq[frame:]):
                start = frame + offset  # zero based
                answer.append((start, start + len(n), +1, n, t))
    if options.strand != "forward":
        rc = reverse_complement(nuc_seq)
        for frame in range(0, 3):
            for offset, n, t in break_up_frame(rc[frame:]):
                start = full_len - frame - offset  # zero based
                answer.append((start - len(n), start, -1, n, t))
    answer.sort()
    return answer


def get_top_peptides(nuc_seq):
    """Returns all peptides of max length."""
    values = list(get_all_peptides(nuc_seq))
    if not values:
        raise StopIteration
    max_len = max(len(x[-1]) for x in values)
    for x in values:
        if len(x[-1]) == max_len:
            yield x


def get_one_peptide(nuc_seq):
    """Returns first (left most) peptide with max length."""
    values = list(get_top_peptides(nuc_seq))
    if not values:
        raise StopIteration
    yield values[0]


if options.mode == "all":
    get_peptides = get_all_peptides
elif options.mode == "top":
    get_peptides = get_top_peptides
elif options.mode == "one":
    get_peptides = get_one_peptide

in_count = 0
out_count = 0
if options.out_nuc_file == "-":
    out_nuc = sys.stdout
elif options.out_nuc_file:
    out_nuc = open(options.out_nuc_file, "w")
else:
    out_nuc = None

if options.out_prot_file == "-":
    out_prot = sys.stdout
elif options.out_prot_file:
    out_prot = open(options.out_prot_file, "w")
else:
    out_prot = None

if options.out_bed_file == "-":
    out_bed = sys.stdout
elif options.out_bed_file:
    out_bed = open(options.out_bed_file, "w")
else:
    out_bed = None

if options.out_gff3_file == "-":
    out_gff3 = sys.stdout
elif options.out_gff3_file:
    out_gff3 = open(options.out_gff3_file, "w")
else:
    out_gff3 = None

if out_gff3:
    out_gff3.write('##gff-version 3\n')

for record in SeqIO.parse(options.input_file, seq_format):
    for i, (f_start, f_end, f_strand, n, t) in enumerate(get_peptides(str(record.seq).upper())):
        out_count += 1
        if f_strand == +1:
            loc = "%i..%i" % (f_start + 1, f_end)
        else:
            loc = "complement(%i..%i)" % (f_start + 1, f_end)
        descr = "length %i aa, %i bp, from %s of %s" \
                % (len(t), len(n), loc, record.description)
        fid = record.id + "|%s%i" % (options.ftype, i + 1)
        r = SeqRecord(Seq(n), id=fid, name="", description=descr)
        t = SeqRecord(Seq(t), id=fid, name="", description=descr)
        if out_nuc:
            SeqIO.write(r, out_nuc, "fasta")
        if out_prot:
            SeqIO.write(t, out_prot, "fasta")
        nice_strand = '+' if f_strand == +1 else '-'
        if out_bed:
            out_bed.write('\t'.join(map(str, [record.id, f_start, f_end, fid, 0, nice_strand])) + '\n')
        if out_gff3:
            out_gff3.write('\t'.join(map(str, [record.id, 'getOrfsOrCds', 'CDS', f_start + 1, f_end, '.',
                                               nice_strand, 0, 'ID=%s%s' % (options.ftype, i + 1)])) + '\n')
    in_count += 1
if out_nuc and out_nuc is not sys.stdout:
    out_nuc.close()
if out_prot and out_prot is not sys.stdout:
    out_prot.close()
if out_bed and out_bed is not sys.stdout:
    out_bed.close()

print("Found %i %ss in %i sequences" % (out_count, options.ftype, in_count))
