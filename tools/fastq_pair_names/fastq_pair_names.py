#!/usr/bin/env python
"""Extract paired read names from FASTQ file(s).

The input file should be a valid FASTQ file(s), the output is two tabular
files - the paired read names (without suffixes), and unpaired read names
(including any unrecognised pair names).

Note that the FASTQ variant is unimportant (Sanger, Solexa, Illumina, or even
Color Space should all work equally well).

This script is copyright 2014-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI), Scotland, UK. All rights reserved.

See accompanying text file for licence details (MIT license).
"""

from __future__ import print_function

import os
import re
import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print("0.0.5")
    sys.exit(0)

from galaxy_utils.sequence.fastq import fastqReader

msg = """Expects at least 3 arguments:

 - Pair names tabular output filename
 - Non-pair names tabular output filename
 - Input FASTQ input filename(s)
"""

if len(sys.argv) < 3:
    sys.exit(msg)

output_pairs = sys.argv[1]
output_nonpairs = sys.argv[2]
input_fastq_filenames = sys.argv[3:]

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

re_illumina_f = re.compile(r"^@[a-zA-Z0-9_:-]+ 1:.*$")
re_illumina_r = re.compile(r"^@[a-zA-Z0-9_:-]+ 2:.*$")
assert re_illumina_f.match("@HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 1:N:0:TGNCCA")
assert re_illumina_r.match("@HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 2:N:0:TGNCCA")
assert not re_illumina_f.match(
    "@HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 2:N:0:TGNCCA"
)
assert not re_illumina_r.match(
    "@HWI-ST916:79:D04M5ACXX:1:1101:10000:100326 1:N:0:TGNCCA"
)

count = 0
pairs = set()  # Will this scale OK?
forward = 0
reverse = 0
neither = 0

out_pairs = open(output_pairs, "w")
out_nonpairs = open(output_nonpairs, "w")

for input_fastq in input_fastq_filenames:
    if not os.path.isfile(input_fastq):
        sys.exit("Missing input FASTQ file %r" % input_fastq)
    in_handle = open(input_fastq)

    # Don't care about the FASTQ type really...
    for record in fastqReader(in_handle, "sanger"):
        count += 1
        name = record.identifier.split(None, 1)[0]
        assert name[0] == "@", record.identifier  # Quirk of the Galaxy parser
        name = name[1:]
        is_forward = False
        suffix = re_f.search(name)
        if suffix:
            # ============
            # Forward read
            # ============
            template = name[: suffix.start()]
            is_forward = True
        elif re_illumina_f.match(record.identifier):
            template = name  # No suffix
            is_forward = True
        if is_forward:
            forward += 1
            if template not in pairs:
                pairs.add(template)
                out_pairs.write(template + "\n")
        else:
            is_reverse = False
            suffix = re_r.search(name)
            if suffix:
                # ============
                # Reverse read
                # ============
                template = name[: suffix.start()]
                is_reverse = True
            elif re_illumina_r.match(record.identifier):
                template = name  # No suffix
                is_reverse = True
            if is_reverse:
                reverse += 1
                if template not in pairs:
                    pairs.add(template)
                    out_pairs.write(template + "\n")
            else:
                # ===========================
                # Neither forward nor reverse
                # ===========================
                out_nonpairs.write(name + "\n")
                neither += 1
    in_handle.close()
out_pairs.close()
out_nonpairs.close()

print(
    "%i reads (%i forward, %i reverse, %i neither), %i pairs"
    % (count, forward, reverse, neither, len(pairs))
)
