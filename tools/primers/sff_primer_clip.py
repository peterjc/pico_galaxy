#!/usr/bin/env python
"""Looks for the given primer sequences and clips matching SFF reads.

Takes seven command line options, input SFF read filename, input primer FASTA
filename, type of primers (forward, reverse or reverse-complement), number of
mismatches (currently only 0, 1 and 2 are supported), minimum length to keep a
read (after primer trimming), should primer-less reads be kept, and finally the
output SFF filename.

This can also be used for stripping off (and optionally filtering on) barcodes.

Note that only the trim/clip values in the SFF file are changed, not the flow
information of the full read sequence.

This script is copyright 2011 by Peter Cock, SCRI, UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script.
"""
import sys
from core_primer_clip import stop_err, load_primers_as_re

try:
    from Bio.SeqIO.SffIO import SffIterator, SffWriter
except ImportError:
    stop_err("Requires Biopython 1.54 or later")

try:
    from Bio.SeqIO.SffIO import ReadRocheXmlManifest
except ImportError:
    #Prior to Biopython 1.56 this was a private function
    from Bio.SeqIO.SffIO import _sff_read_roche_index_xml as ReadRocheXmlManifest


#Parse Command Line
try:
    input_sff, primer_fasta, primer_type, mm, min_len, keep_negatives, output_sff = sys.argv[1:]
except ValueError:
    stop_err("Expected seven arguments, got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))

try:
    mm = int(mm)
except ValueError:
    stop_err("Expected non-negative integer number of mismatches (e.g. 0 or 1), not %r" % mm)
if mm < 0:
    stop_err("Expected non-negtive integer number of mismatches (e.g. 0 or 1), not %r" % mm)
if mm not in [0,1,2]:
    raise NotImplementedError

try:
    min_len = int(min_len)
except ValueError:
    stop_err("Expected non-negative integer min_len (e.g. 0 or 1), not %r" % min_len)
if min_len < 0:
    stop_err("Expected non-negtive integer min_len (e.g. 0 or 1), not %r" % min_len)


if keep_negatives.lower() in ["true", "yes", "on"]:
    keep_negatives = True
elif keep_negatives.lower() in ["false", "no", "off"]:
    keep_negatives = False
else:
    stop_err("Expected boolean for keep_negatives (e.g. true or false), not %r" % keep_negatives)


if primer_type.lower() == "forward":
    forward = True
    rc = False
elif primer_type.lower() == "reverse":
    forward = False
    rc = False
elif primer_type.lower() == "reverse-complement":
    forward = False
    rc = True
else:
    stop_err("Expected foward, reverse or reverse-complement not %r" % primer_type)


#Read primer file and record all specified sequences
count, primer = load_primers_as_re(primer_fasta, mm, rc)
print "%i primer sequences" % count


short = 0
clipped = 0
negs = 0
def process(records):
    global short, clipped, negs
    if forward:
        for record in records:
            left_clip = record.annotations["clip_qual_left"]
            right_clip = record.annotations["clip_qual_right"]
            seq = str(record.seq)[left_clip:right_clip].upper()
            result = primer.search(seq)
            if result:
                #Forward primer, take everything after it
                #so move the left clip along
                record.annotations["clip_qual_left"] = left_clip + result.end()
                if right_clip - left_clip - result.end() >= min_len:
                    clipped += 1
                    yield record
                else:
                    short += 1
            elif keep_negatives:
                negs += 1
                yield record
    else:
        for record in records:
            left_clip = record.annotations["clip_qual_left"]
            right_clip = record.annotations["clip_qual_right"]
            seq = str(record.seq)[left_clip:right_clip].upper()
            result = primer.search(seq)
            if result:
                #Reverse primer, take everything before it
                #so move the right clip back
                record.annotations["clip_qual_right"] = left_clip + result.start()
                if right_clip - left_clip - result.start() >= min_len:
                    clipped += 1
                    yield record
                else:
                    short += 1
            elif keep_negatives:
                negs += 1
                yield record

in_handle = open(input_sff, "rb")
try:
    manifest = ReadRocheXmlManifest(in_handle)
except ValueError:
    manifest = None
in_handle.seek(0)
out_handle = open(output_sff, "wb")
writer = SffWriter(out_handle, xml=manifest)
writer.write_file(process(SffIterator(in_handle)))
in_handle.close()
out_handle.close()

print "Kept %i clipped reads" % clipped
print "Discarded %i short reads" % short
if keep_negatives:
    print "Kept %i non-matching reads" % negs
