#!/usr/bin/env python
"""MUMmer wrapper calling nucmer/promer/mummerplot etc.

Takes the following command line options,
1. FASTA filename of species A
2. FASTA filename of species B
3. Algorithm, nucmer or promer
4. PNG output filename
5. PDF output filename
"""

import os
import shutil
import sys
import tempfile


def run(cmd):
    print(cmd)
    return_code = os.system(cmd)
    if return_code:
        sys.exit("Error %i from: %s" % (return_code, cmd))


if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print("MUMmer wrapper v0.0.8\n")
    # TODO - How to get a version string from the mummer binary?
    os.system("nucmer --version")
    os.system("promer --version")
    os.system("mummerplot --version")
    os.system("gnuplot --version")
    # TODO - Should we include "gs --version" as a proxy for ps2pdf?
    sys.exit(0)

# Parse Command Line
# TODO - optparse
try:
    fasta_a, fasta_b, algorithm, png_out, pdf_out = sys.argv[1:]
except ValueError:
    sys.exit(
        "Expect 5 arguments (FASTA, FASTA, algorithm, PNG out, PDF out), got %i"
        % (len(sys.argv) - 1)
    )


valid_algo = ["mummer", "nucmer", "promer"]
if algorithm not in valid_algo:
    sys.exit(
        "Invalid algorithm argument %r, should be: %s"
        % (algorithm, ", ".join(valid_algo))
    )

base_path = tempfile.mkdtemp()
prefix = os.path.join(base_path, "ref_qry")
coords = prefix + ".mums"
# gnuplot = prefix + ".gp"
ps_image = prefix + ".ps"
png_image = prefix + ".png"

if algorithm == "mummer":
    # Add -mum as per example to find maximal unique matches between ref and query.
    # Add the -b -c options to search both strands and report relative to forward strand
    # which then matches the default dual-strand approach in nucmer and promer
    cmd = '%s -mum -b -c "%s" "%s" > %s' % (algorithm, fasta_a, fasta_b, coords)
else:
    coords = "out.delta"
    cmd = '%s "%s" "%s"' % (algorithm, fasta_a, fasta_b)
run(cmd)

output_failed = False

# PNG
# ===
cmd = 'mummerplot -R "%s" -Q "%s" --png --large --prefix=%s %s' % (
    fasta_a,
    fasta_b,
    prefix,
    coords,
)
run(cmd)
if os.path.isfile(png_image):
    shutil.move(png_image, png_out)
else:
    sys.stderr.write("ERROR: PNG file not created.\n")
    output_failed = True

# PS --> PDF
# ==========
# Using --large, puts "set size 3,3" in the gnuplot - which seems to mess up.
# Problem here is the default bbox (BoundingBox) in the PS output is letter page size,
# and even if we override that, at least when view the PS output in Adobe Illustrator
# things don't seem to be lined up properly :(
#
# Using "set size 1,1" works better - which is what --small gives:
cmd = 'mummerplot -R "%s" -Q "%s" --postscript --small --prefix=%s %s' % (
    fasta_a,
    fasta_b,
    prefix,
    coords,
)
run(cmd)
if not os.path.isfile(ps_image):
    sys.stderr.write("ERROR: PostScript file needed for PDF output was not created.\n")
    output_failed = True
else:
    cmd = 'ps2pdf -dEPSCrop "%s" "%s"' % (ps_image, pdf_out)
    run(cmd)
    if not os.path.isfile(pdf_out):
        sys.stderr.write("ERROR: PDF file not created.\n")
        output_failed = True

# Remove temp files...
os.remove(coords)  # Might not be under the temp directory...
shutil.rmtree(base_path)

if output_failed:
    sys.exit("ERROR: Failed to produce output file(s).")
