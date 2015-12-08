#!/usr/bin/env python
"""Divides a FASTQ into paired and single (orphan reads) as separate files.

The input file should be a valid FASTQ file which has been sorted so that
any partner forward+reverse reads are consecutive. The output files all
preserve this sort order. Pairing are recognised based on standard name
suffices. See below or run the tool with no arguments for more details.

Note that the FASTQ variant is unimportant (Sanger, Solexa, Illumina, or even
Color Space should all work equally well).

This script is copyright 2010-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI), Scotland, UK. All rights reserved.

See accompanying text file for licence details (MIT license).
"""
import sys
import re

if "-v" in sys.argv or "--version" in sys.argv:
    print("Version 0.1.0")
    sys.exit(0)

try:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError:
    sys.exit("Biopython missing")

msg = """Expect either 3 or 4 arguments, all FASTQ filenames.

If you want two output files, use four arguments:
 - FASTQ variant (e.g. sanger, solexa, illumina or cssanger)
 - Sorted input FASTQ filename,
 - Output paired FASTQ filename (forward then reverse interleaved),
 - Output singles FASTQ filename (orphan reads)

If you want three output files, use five arguments:
 - FASTQ variant (e.g. sanger, solexa, illumina or cssanger)
 - Sorted input FASTQ filename,
 - Output forward paired FASTQ filename,
 - Output reverse paired FASTQ filename,
 - Output singles FASTQ filename (orphan reads)

The input file should be a valid FASTQ file which has been sorted so that
any partner forward+reverse reads are consecutive. The output files all
preserve this sort order.

Any reads where the forward/reverse naming suffix used is not recognised
are treated as orphan reads. The tool supports the /1 and /2 convention
originally used by Illumina, the .f and .r convention, and the Sanger
convention (see http://staden.sourceforge.net/manual/pregap4_unix_50.html
for details), and the new Illumina convention where the reads have the
same identifier with the fragment at the start of the description, e.g.

@HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 1:N:0:TGNCCA
@HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 2:N:0:TGNCCA

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

if len(sys.argv) == 5:
    format, input_fastq, pairs_fastq, singles_fastq = sys.argv[1:]
elif len(sys.argv) == 6:
    pairs_fastq = None
    format, input_fastq, pairs_f_fastq, pairs_r_fastq, singles_fastq = sys.argv[1:]
else:
    sys.exit(msg)

format = format.replace("fastq", "").lower()
if not format:
    format = "sanger"  # safe default
elif format not in ["sanger", "solexa", "illumina", "cssanger"]:
    sys.exit("Unrecognised format %s" % format)

# Cope with three widely used suffix naming convensions,
# Illumina: /1 or /2
# Forward/revered: .f or .r
# Sanger, e.g. .p1k and .q1k
# See http://staden.sourceforge.net/manual/pregap4_unix_50.html
re_f = re.compile(r"(/1|\.f|\.[sfp]\d\w*)$")
re_r = re.compile(r"(/2|\.r|\.[rq]\d\w*)$")

# assert re_f.match("demo/1")
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

re_illumina_f = re.compile(r"^[a-zA-Z0-9_:-]+ 1:.*$")
re_illumina_r = re.compile(r"^[a-zA-Z0-9_:-]+ 2:.*$")
assert re_illumina_f.match("HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 1:N:0:TGNCCA")
assert re_illumina_r.match("HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 2:N:0:TGNCCA")
assert not re_illumina_f.match("HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 2:N:0:TGNCCA")
assert not re_illumina_r.match("HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 1:N:0:TGNCCA")

FASTQ_TEMPLATE = "@%s\n%s\n+\n%s\n"

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

for title, seq, qual in FastqGeneralIterator(in_handle):
    count += 1
    name = title.split(None, 1)[0]
    is_forward = False
    suffix = re_f.search(name)
    if suffix:
        # ============
        # Forward read
        # ============
        template = name[:suffix.start()]
        is_forward = True
    elif re_illumina_f.match(title):
        template = name  # No suffix
        is_forward = True
    if is_forward:
        # print(name, "forward", template)
        forward += 1
        if last_template == template:
            buffered_reads.append((title, seq, qual))
        else:
            # Any old buffered reads are orphans
            for old in buffered_reads:
                singles_handle.write(FASTQ_TEMPLATE % old)
                singles += 1
            # Save this read in buffer
            buffered_reads = [(title, seq, qual)]
            last_template = template
    else:
        is_reverse = False
        suffix = re_r.search(name)
        if suffix:
            # ============
            # Reverse read
            # ============
            template = name[:suffix.start()]
            is_reverse = True
        elif re_illumina_r.match(title):
            template = name  # No suffix
            is_reverse = True
        if is_reverse:
            # print(name, "reverse", template)
            reverse += 1
            if last_template == template and buffered_reads:
                # We have a pair!
                # If there are multiple buffered forward reads, want to pick
                # the first one (although we could try and do something more
                # clever looking at the suffix to match them up...)
                old = buffered_reads.pop(0)
                pairs_f_handle.write(FASTQ_TEMPLATE % old)
                pairs_r_handle.write(FASTQ_TEMPLATE % (title, seq, qual))
                pairs += 2
            else:
                # As this is a reverse read, this and any buffered read(s) are
                # all orphans
                for old in buffered_reads:
                    singles_handle.write(FASTQ_TEMPLATE % old)
                    singles += 1
                buffered_reads = []
                singles_handle.write(FASTQ_TEMPLATE % (title, seq, qual))
                singles += 1
                last_template = None
        else:
            # ===========================
            # Neither forward nor reverse
            # ===========================
            singles_handle.write(FASTQ_TEMPLATE % (title, seq, qual))
            singles += 1
            neither += 1
            for old in buffered_reads:
                singles_handle.write(FASTQ_TEMPLATE % old)
                singles += 1
            buffered_reads = []
            last_template = None
if last_template:
    # Left over singles...
    for old in buffered_reads:
        singles_handle.write(FASTQ_TEMPLATE % old)
        singles += 1
in_handle.close()
singles_handle.close()
if pairs_fastq:
    pairs_f_handle.close()
    assert pairs_r_handle.closed
else:
    pairs_f_handle.close()
    pairs_r_handle.close()

if neither:
    print("%i reads (%i forward, %i reverse, %i neither), %i in pairs, %i as singles"
          % (count, forward, reverse, neither, pairs, singles))
else:
    print("%i reads (%i forward, %i reverse), %i in pairs, %i as singles"
          % (count, forward, reverse, pairs, singles))

assert count == pairs + singles == forward + reverse + neither, \
    "%i vs %i+%i=%i vs %i+%i+%i=%i" \
    % (count, pairs, singles, pairs + singles, forward, reverse, neither, forward + reverse + neither)
