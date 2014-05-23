#!/usr/bin/env python
"""MUMmer wrapper calling nucmer/promer/mummerplot etc,

Takes the following command line options,
1. FASTA filename of species A
2. FASTA filename of species B
3. Algorithm, nucmer or promer
4. PNG output filename
5. PDF output filename
"""

import os
import sys
import tempfile
import shutil

def stop_err( msg ):
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)

def run(cmd):
    print(cmd)
    return_code = os.system(cmd)
    if return_code:
        stop_err("Error %i from: %s" % (return_code, cmd))

if "-v" in sys.argv [1:]or "--version" in sys.argv[1:]:
    print("MUMmer wrapper v0.0.1\n")
    os.system("nucmer --version")
    os.system("promer --version")
    os.system("mummerplot --version")
    sys.exit(0)

#Parse Command Line
#TODO - optparse
try:
    fasta_a, fasta_b, algorithm, png_out, pdf_out = sys.argv[1:]
except:
    stop_err("Expect 5 arguments, got %i" % (len(sys.argv) - 1))


valid_algo = ["mummer", "nucmer", "promer"]
if algorithm not in valid_algo:
    stop_err("Invalid algorithm argument %r, should be: %s" % (algorithm, ", ".join(valid_algo)))

base_path = tempfile.mkdtemp()
prefix = os.path.join(base_path, "ref_qry")
coords = prefix + ".mums"
#gnuplot = prefix + ".gp"
ps_image = prefix + ".ps"
png_image = prefix + ".png"

if algorithm == "mummer":
    #Add the -b -c options to search both strands and report relative to forward strand
    #which then matches the default dual-strand approach in nucmer and promer
    cmd = '%s -b -c "%s" "%s" > %s' % (algorithm, fasta_a, fasta_b, coords)
else:
    coords = "out.delta"
    cmd = '%s "%s" "%s"' % (algorithm, fasta_a, fasta_b)
run(cmd)

# PNG
# ===
cmd = 'mummerplot -R "%s" -Q "%s" --png --large --prefix=%s %s' % (fasta_a, fasta_b, prefix, coords)
run(cmd)
shutil.move(png_image, png_out)

# PS --> PDF
# ==========
# Using --large, puts "set size 3,3" in the gnuplot - which seems to mess up.
# Problem here is the default bbox (BoundingBox) in the PS output is letter page size,
# and even if we override that, at least when view the PS output in Adobe Illustrator
# things don't seem to be lined up properly :(
#
# Using "set size 1,1" works better - which is what --small gives:
cmd = 'mummerplot -R "%s" -Q "%s" --postscript --small --prefix=%s %s' % (fasta_a, fasta_b, prefix, coords)
run(cmd)
cmd = 'ps2pdf -dEPSCrop "%s" "%s"' % (ps_image, pdf_out)
run(cmd)

#Remove temp files...
shutil.rmtree(base_path)
