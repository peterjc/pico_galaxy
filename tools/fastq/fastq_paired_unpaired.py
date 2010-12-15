#!/usr/bin/env python
"""Divides a FASTQ into paired and single (orphan reads) as separate files.

This tool is a short Python script which requires Biopython 1.50 or later.
If you use this tool in scientific work leading to a publication, please cite
the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

The input file should be a valid FASTQ file which has been sorted so that
any partner forward+reverse reads are consecutive. The output files all
preserve this sort order. Pairing are recognised based on standard name
suffices. See below or run the tool with no arguments for more details.

Note that the FASTQ variant is unimportant (Sanger, Solexa, Illumina, or even
Color Space should all work equally well).

This script is copyright 2010 by Peter Cock, SCRI, UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.3 of the script.
"""
import sys
import re

def stop_err(msg, err=1):
   sys.stderr.write(msg.rstrip() + "\n")
   sys.exit(err)

try:
    from Bio import SeqIO
except ImportError:
    stop_err("Biopython not installed")

msg = """Expect either 3 or 4 arguments, all FASTQ filenames.

If you want two output files, use three arguments:
 - Sorted input FASTQ filename,
 - Output paired FASTQ filename (forward then reverse interleaved),
 - Output singles FASTQ filename (orphan reads)

If you want three output files, use four arguments:
 - Sorted input FASTQ filename,
 - Output forward paired FASTQ filename,
 - Output reverse paired FASTQ filename,
 - Output singles FASTQ filename (orphan reads)

The input file should be a valid FASTQ file which has been sorted so that
any partner forward+reverse reads are consecutive. The output files all
preserve this sort order.

Any reads where the forward/reverse naming suffix used is not recognised
are treated as orphan reads. The tool supports the /1 and /2 convention
used by Illumina, the .f and .r convention, and the Sanger convention
(see http://staden.sourceforge.net/manual/pregap4_unix_50.html for details).

Note that this does support multiple forward and reverse reads per template
(which is quite common with Sanger sequencing), e.g. this which is sorted
alphabetically:

WTSI_1055_4p17.p1kapIBF
WTSI_1055_4p17.p1kpIBF
WTSI_1055_4p17.q1kapIBR
WTSI_1055_4p17.q1kpIBR

or this where the reads already come in pairs:

WTSI_1055_4p17.p1kapIBF
WTSI_1055_4p17.q1kapIBR
WTSI_1055_4p17.p1kpIBF
WTSI_1055_4p17.q1kpIBR

both become:

WTSI_1055_4p17.p1kapIBF paired with WTSI_1055_4p17.q1kapIBR
WTSI_1055_4p17.p1kpIBF paired with WTSI_1055_4p17.q1kpIBR
"""

if len(sys.argv) == 4:
    input_fastq, pairs_fastq, singles_fastq = sys.argv[1:]
elif len(sys.argv) == 5:
    pairs_fastq = None
    input_fastq, pairs_f_fastq, pairs_r_fastq, singles_fastq = sys.argv[1:]
else:
    stop_err(msg)

def f_match(name):
   if name.endswith("/1") or name.endswith(".f"):
      return True

#Cope with three widely used suffix naming convensions,
#Illumina: /1 or /2
#Forward/revered: .f or .r
#Sanger, e.g. .p1k and .q1k
#See http://staden.sourceforge.net/manual/pregap4_unix_50.html
re_f = re.compile(r"(/1|\.f|\.[sfp]\d\w*)$")
re_r = re.compile(r"(/2|\.r|\.[rq]\d\w*)$")

#assert re_f.match("demo/1")
assert re_f.search("demo.f")
assert re_f.search("demo.s1")
assert re_f.search("demo.f1k")
assert re_f.search("demo.p1")
assert re_f.search("demo.p1k")
assert re_f.search("demo.p1lk")
assert re_r.search("demo/2")
assert re_r.search("demo.r")
assert re_r.search("demo.q1")
assert re_r.search("demo.q1lk")
assert not re_r.search("demo/1")
assert not re_r.search("demo.f")
assert not re_r.search("demo.p")
assert not re_f.search("demo/2")
assert not re_f.search("demo.r")
assert not re_f.search("demo.q")

count, forward, reverse, neither, pairs, singles = 0, 0, 0, 0, 0, 0
in_handle = open(input_fastq)
if pairs_fastq:
    pairs_f_handle = open(pairs_fastq, "w")
    pairs_r_handle = pairs_f_handle
else:
    pairs_f_handle = open(pairs_f_fastq, "w")
    pairs_r_handle = open(pairs_r_fastq, "w")
singles_handle = open(singles_fastq, "w")
last_template, buffered_reads = None, []
for title, seq, qual in SeqIO.QualityIO.FastqGeneralIterator(in_handle):
    count += 1
    name = title.split(None,1)[0]
    suffix = re_f.search(name)
    if suffix:
        #============
        #Forward read
        #============
        template = name[:suffix.start()]
        #print name, "forward", template
        forward += 1
        if last_template == template:
            buffered_reads.append((title, seq, qual))
        else:
            #Any buffered reads are orphans
            for read in buffered_reads:
                singles_handle.write("@%s\n%s\n+\n%s\n" % read)
                singles += 1
            buffered_reads = [(title, seq, qual)]
            last_template = template  
    else:
        suffix = re_r.search(name)
        if suffix:
            #============
            #Reverse read
            #============
            template = name[:suffix.start()]
            #print name, "reverse", template
            reverse += 1
            if last_template == template and buffered_reads:
                #We have a pair!
                #If there are multiple buffered forward reads, want to pick
                #the first one (although we could try and do something more
                #clever looking at the suffix to match them up...)
                read = buffered_reads.pop(0)
                pairs_f_handle.write("@%s\n%s\n+\n%s\n" % read)
                pairs_r_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                pairs += 2
            else:
                #As this is a reverse read, this and any buffered read(s) are
                #all orphans
                for read in buffered_reads:
                    singles_handle.write("@%s\n%s\n+\n%s\n" % read)
                    singles += 1
                buffered_reads = []
                singles_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                singles += 1
                last_template = None
        else:
            #===========================
            #Neither forward nor reverse
            #===========================
            singles_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            singles += 1
            neither += 1
            for read in buffered_reads:
                singles_handle.write("@%s\n%s\n+\n%s\n" % read)
                singles += 1
            buffered_reads = []
            last_template = None
if last_template:
    #Left over singles...
    for read in buffered_reads:
        singles_handle.write("@%s\n%s\n+\n%s\n" % read)
        singles += 1
in_handle.close
singles_handle.close()
if pairs_fastq:
    pairs_f_handle.close()
    assert pairs_r_handle.closed
else:
    pairs_f_handle.close()
    pairs_r_handle.close()

if neither:
    print "%i reads (%i forward, %i reverse, %i neighter), %i in pairs, %i as singles" \
           % (count, forward, reverse, neither, pairs, singles)
else:
    print "%i reads (%i forward, %i reverse), %i in pairs, %i as singles" \
           % (count, forward, reverse, pairs, singles)

assert count == pairs + singles == forward + reverse + neither, \
       "%i vs %i+%i=%i vs %i+%i=%i" \
       % (count,pairs,singles,pairs+singles,forward,reverse,neither,forward+reverse+neither)
