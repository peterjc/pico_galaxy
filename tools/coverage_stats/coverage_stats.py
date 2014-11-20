#!/usr/bin/env python
"""BAM coverage statistics using samtools idxstats and depth.

This script takes exactly three command line arguments:
 * Input BAM filename
 * Input BAI filename (via Galaxy metadata)
 * Output tabular filename

This messes about with the filenames to make samtools happy, then
runs "samtools idxstats" and "samtools depth", captures and combines
the output to the desired output tabular file.
"""
import sys
import os
import subprocess
import tempfile

if "-v" in sys.argv or "--version" in sys.argv:
    #Galaxy seems to invert the order of the two lines
    print "BAM coverage statistics v0.0.1"
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

def stop_err(msg, error_level=1):
   """Print error message to stdout and quit with given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)

if len(sys.argv) != 4:
   stop_err("Require three arguments: BAM, BAI, tabular filenames")

bam_filename, bai_filename, tabular_filename = sys.argv[1:]

if not os.path.isfile(bam_filename):
    stop_err("Input BAM file not found: %s" % bam_filename)
if not os.path.isfile(bai_filename):
    if bai_filename == "None":
        stop_err("Error: Galaxy did not index your BAM file")
    stop_err("Input BAI file not found: %s" % bai_filename)

#Assign sensible names with real extensions, and setup symlinks:
tmp_dir = tempfile.mkdtemp()
bam_file = os.path.join(tmp_dir, "temp.bam")
bai_file = os.path.join(tmp_dir, "temp.bam.bai")
idxstats_filename = os.path.join(tmp_dir, "idxstats.tsv")
depth_filename = os.path.join(tmp_dir, "depth.tsv")
os.symlink(os.path.abspath(bam_filename), bam_file)
os.symlink(os.path.abspath(bai_filename), bai_file)
assert os.path.isfile(bam_file), bam_file
assert os.path.isfile(bai_file), bai_file
assert os.path.isfile(bam_file + ".bai"), bam_file

def clean_up():
    os.remove(idxstats_filename)
    os.remove(depth_filename)
    os.remove(bam_file)
    os.remove(bai_file)
    os.rmdir(tmp_dir)

# Run samtools idxstats:
cmd = 'samtools idxstats "%s" > "%s"' % (bam_file, idxstats_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    stop_err("Return code %i from command:\n%s" % (return_code, cmd))

# Run samtools depth:
# TODO - Parse stdout instead?
cmd = 'samtools depth "%s" > "%s"' % (bam_file, depth_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    stop_err("Return code %i from command:\n%s" % (return_code, cmd))

# Parse and combine the output
handle = open(tabular_filename, "w")
handle.write("#TODO\n")
handle.close()

# Remove the temp symlinks and files:
clean_up()
