#!/usr/bin/env python
"""A simple wrapper script to call MIRA and collect its output.

This focuses on the miraconvert binary.
"""
import os
import sys
import subprocess
import shutil
import time
import tempfile
from optparse import OptionParser

#Do we need any PYTHONPATH magic?
from mira4_make_bam import make_bam

WRAPPER_VER = "0.0.5" #Keep in sync with the XML file

def stop_err(msg, err=1):
    sys.stderr.write(msg+"\n")
    sys.exit(err)


def get_version(mira_binary):
    """Run MIRA to find its version number"""
    # At the commend line I would use: mira -v | head -n 1
    # however there is some pipe error when doing that here.
    cmd = [mira_binary, "-v"]
    try:
        child = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
    except Exception, err:
        sys.stderr.write("Error invoking command:\n%s\n\n%s\n" % (" ".join(cmd), err))
        sys.exit(1)
    ver, tmp = child.communicate()
    del child
    return ver.split("\n", 1)[0].strip()

#Parse Command Line
usage = """Galaxy MIRA4 wrapper script v%s - use as follows:

$ python mira4_convert.py ...

This will run the MIRA miraconvert binary and collect its output files as directed.
""" % WRAPPER_VER
parser = OptionParser(usage=usage)
parser.add_option("--input", dest="input",
                  default=None, metavar="FILE",
                  help="MIRA input filename")
parser.add_option("--min_length", dest="min_length",
                  default="0",
                  help="Minimum contig length")
parser.add_option("--min_reads", dest="min_reads",
                  default="0",
                  help="Minimum reads per contig")
parser.add_option("--maf", dest="maf",
                  default="-", metavar="FILE",
                  help="MIRA MAF output filename")
parser.add_option("--ace", dest="ace",
                  default="-", metavar="FILE",
                  help="ACE output filename")
parser.add_option("--sam", dest="sam",
                  default="-", metavar="FILE",
                  help="Unpadded SAM output filename")
parser.add_option("--bam", dest="bam",
                  default="-", metavar="FILE",
                  help="Unpadded BAM output filename")
parser.add_option("--fasta", dest="fasta",
                  default="-", metavar="FILE",
                  help="Unpadded FASTA output filename")
options, args = parser.parse_args()
input_maf = options.input
out_maf = options.maf
out_sam = options.sam
out_bam = options.bam
out_fasta = options.fasta
out_ace = options.ace

try:
    mira_path = os.environ["MIRA4"]
except KeyError:
    stop_err("Environment variable $MIRA4 not set")
mira_convert = os.path.join(mira_path, "miraconvert")
if not os.path.isfile(mira_convert):
    stop_err("Missing miraconvert under $MIRA4, %r\nFolder contained: %s"
             % (mira_convert, ", ".join(os.listdir(mira_path))))

mira_convert_ver = get_version(mira_convert)
if not mira_convert_ver.strip().startswith("4.0"):
    stop_err("This wrapper is for MIRA V4.0, not:\n%s\n%s" % (mira_ver, mira_convert))
if "-v" in sys.argv or "--version" in sys.argv:
    print "%s, MIRA wrapper version %s" % (mira_convert_ver, WRAPPER_VER)
    sys.exit(0)

if not input_maf:
    stop_err("Input MIRA file is required")
elif not os.path.isfile(input_maf):
    stop_err("Missing input MIRA file: %r" % input_maf)

print("TODO...")
sys.exit(1)
print("Done")
