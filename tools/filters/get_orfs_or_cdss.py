#!/usr/bin/env python
"""Find ORFs in a nucleotide sequence file.

get_orfs_or_cdss.py $input_fasta $input_format $table $ftype $ends $mode $min_len $strand $out_nuc_file $out_prot_file

Takes ten command line options, input sequence filename, format, genetic
code, CDS vs ORF, end type (open, closed), selection mode (all, top, one),
minimum length (in amino acids), strand (both, forward, reverse), output
nucleotide filename, and output protein filename.

This tool is a short Python script which requires Biopython. If you use
this tool in scientific work leading to a publication, please cite the
Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2011 by Peter Cock, The James Hutton Institute
(formerly SCRI), Dundee, UK. All rights reserved.

See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script.
"""
import sys
import re

def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

try:
    from Bio.Seq import Seq, reverse_complement, translate
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    from Bio.Data import CodonTable
except ImportError:
    stop_err("Missing Biopython library")

#Parse Command Line
try:
    input_file, seq_format, table, ftype, ends, mode, min_len, strand, out_nuc_file, out_prot_file = sys.argv[1:]
except ValueError:
    stop_err("Expected ten arguments, got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))

try:
    table = int(table)
except ValueError:
    stop_err("Expected integer for genetic code table, got %s" % table)

try:
    table_obj = CodonTable.ambiguous_generic_by_id[table]
except KeyError:
    stop_err("Unknown codon table %i" % table)

if ftype not in ["CDS", "ORF"]:
    stop_err("Expected CDS or ORF, got %s" % ftype)

if ends not in ["open", "closed"]:
    stop_err("Expected open or closed for end treatment, got %s" % ends)

try:
    min_len = int(min_len)
except ValueError:
    stop_err("Expected integer for min_len, got %s" % min_len)

if seq_format.lower()=="sff":
    seq_format = "sff-trim"
elif seq_format.lower()=="fasta":
    seq_format = "fasta"
elif seq_format.lower().startswith("fastq"):
    seq_format = "fastq"
else:
    stop_err("Unsupported file type %r" % seq_format)

print "Genetic code table %i" % table
print "Minimum length %i aa" % min_len
#print "Taking %s ORF(s) from %s strand(s)" % (mode, strand)

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
        #Must check the start is in frame
        start = match.start()
        if start % 3 == 0:
            n = s[start:]
            assert len(n) % 3 == 0, "%s is len %i" % (n, len(n))
            if strict:
                t = translate(n, table, cds=True)
            else:
                #Use when missing stop codon,
                t = "M" + translate(n[3:], table, to_stop=True)
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
        if ftype=="CDS":
            offset, n, t = start_chop_and_trans(n)
        else:
            offset = 0
            t = translate(n, table, to_stop=True)
        if n and len(t) >= min_len:
            yield start + offset, n, t
        start = index
    if ends == "open":
        #No stop codon, Biopython's strict CDS translate will fail
        n = s[start:]
        #Ensure we have whole codons
        #TODO - Try appending N instead?
        #TODO - Do the next four lines more elegantly
        if len(n) % 3:
            n = n[:-1]
        if len(n) % 3:
            n = n[:-1]
        if ftype=="CDS":
            offset, n, t = start_chop_and_trans(n, strict=False)
        else:
            offset = 0
            t = translate(n, table, to_stop=True)
        if n and len(t) >= min_len:
            yield start + offset, n, t
                        

def get_all_peptides(nuc_seq):
    """Returns start, end, strand, nucleotides, protein.

    Co-ordinates are Python style zero-based.
    """
    #TODO - Refactor to use a generator function (in start order)
    #rather than making a list and sorting?
    answer = []
    full_len = len(nuc_seq)
    if strand != "reverse":
        for frame in range(0,3):
            for offset, n, t in break_up_frame(nuc_seq[frame:]):
                start = frame + offset #zero based
                answer.append((start, start + len(n), +1, n, t))
    if strand != "forward":
        rc = reverse_complement(nuc_seq)
        for frame in range(0,3) :
            for offset, n, t in break_up_frame(rc[frame:]):
                start = full_len - frame - offset #zero based
                answer.append((start, start + len(n), -1, n ,t))
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

if mode == "all":
    get_peptides = get_all_peptides
elif mode == "top":
    get_peptides = get_top_peptides
elif mode == "one":
    get_peptides = get_one_peptide

in_count = 0
out_count = 0
if out_nuc_file == "-":
    out_nuc = sys.stdout
else:
    out_nuc = open(out_nuc_file, "w")
if out_prot_file == "-":
    out_prot = sys.stdout
else:
    out_prot = open(out_prot_file, "w")
for record in SeqIO.parse(input_file, seq_format):
    for i, (f_start, f_end, f_strand, n, t) in enumerate(get_peptides(str(record.seq).upper())):
        out_count += 1
        if f_strand == +1:
            loc = "%i..%i" % (f_start+1, f_end)
        else:
            loc = "complement(%i..%i)" % (f_start+1, f_end)
        descr = "length %i aa, %i bp, from %s" % (len(t), len(n), loc)
        r = SeqRecord(Seq(n), id = record.id + "|%s%i" % (ftype, i+1), name = "", description= descr)
        t = SeqRecord(Seq(t), id = record.id + "|%s%i" % (ftype, i+1), name = "", description= descr)
        SeqIO.write(r, out_nuc, "fasta")
        SeqIO.write(t, out_prot, "fasta")
    in_count += 1
if out_nuc is not sys.stdout:
    out_nuc.close()
if out_prot is not sys.stdout:
    out_prot.close()

print "Found %i %ss in %i sequences" % (out_count, ftype, in_count)
