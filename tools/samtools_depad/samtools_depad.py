#!/usr/bin/env python
"""Wrapper for "samtools depad" for use in Galaxy.

This script takes exactly four command line arguments:
 * Input padded reference FASTA file
 * Input SAM/BAM filename
 * Input format ("sam" or "bam")
 * Output BAM filename

Runs "samtools depad" and captures the output to the desired BAM file.
"""
import sys
import os

if "-v" in sys.argv or "--version" in sys.argv:
    # Galaxy seems to invert the order of the two lines
    print "(Galaxy wrapper v0.0.2)"
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

if len(sys.argv) != 5:
    sys.exit("Require four arguments: padded FASTA, SAM/BAM file, format (SAM or BAM), output BAM filenames")

padded_ref, bam_filename, input_format, output_filename = sys.argv[1:]

if not os.path.isfile(padded_ref):
    sys.exit("Input padded reference FASTA file not found: %s" % padded_ref)
if not os.path.isfile(bam_filename):
    sys.exit("Input BAM file not found: %s" % bam_filename)
if input_format.lower() not in ["sam", "bam"]:
    sys.exit("Input format should be SAM or BAM, not %r" % input_format)

# Run samtools depad:
if input_format.lower() == "sam":
    cmd = "samtools depad -S -T %s %s > %s" % (padded_ref, bam_filename, output_filename)
else:
    cmd = "samtools depad -T %s %s > %s" % (padded_ref, bam_filename, output_filename)
return_code = os.system(cmd)

if return_code:
    sys.exit("Return code %i from command:\n%s" % (return_code, cmd))
