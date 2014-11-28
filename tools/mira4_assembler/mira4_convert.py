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
try:
    from io import BytesIO
except ImportError:
    #Should we worry about Python 2.5 or older?
    from StringIO import StringIO as BytesIO

#Do we need any PYTHONPATH magic?
from mira4_make_bam import depad

WRAPPER_VER = "0.0.6" #Keep in sync with the XML file

def stop_err(msg, err=1):
    sys.stderr.write(msg+"\n")
    sys.exit(err)

def run(cmd):
    #Avoid using shell=True when we call subprocess to ensure if the Python
    #script is killed, so too is the child process.
    try:
        child = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception, err:
        stop_err("Error invoking command:\n%s\n\n%s\n" % (" ".join(cmd), err))
    #Use .communicate as can get deadlocks with .wait(),
    stdout, stderr = child.communicate()
    return_code = child.returncode
    if return_code:
        cmd_str = " ".join(cmd)  # doesn't quote spaces etc
        if stderr and stdout:
            stop_err("Return code %i from command:\n%s\n\n%s\n\n%s" % (return_code, cmd_str, stdout, stderr))
        else:
            stop_err("Return code %i from command:\n%s\n%s" % (return_code, cmd_str, stderr))

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
parser.add_option("-x", "--min_length", dest="min_length",
                  default="0",
                  help="Minimum contig length")
parser.add_option("-y", "--min_cover", dest="min_cover",
                  default="0",
                  help="Minimum average contig coverage")
parser.add_option("-z", "--min_reads", dest="min_reads",
                  default="0",
                  help="Minimum reads per contig")
parser.add_option("--maf", dest="maf",
                  default="", metavar="FILE",
                  help="MIRA MAF output filename")
parser.add_option("--ace", dest="ace",
                  default="", metavar="FILE",
                  help="ACE output filename")
parser.add_option("--bam", dest="bam",
                  default="", metavar="FILE",
                  help="Unpadded BAM output filename")
parser.add_option("--fasta", dest="fasta",
                  default="", metavar="FILE",
                  help="Unpadded FASTA output filename")
parser.add_option("--cstats", dest="cstats",
                  default="", metavar="FILE",
                  help="Contig statistics filename")
parser.add_option("-v", "--version", dest="version",
                  default=False, action="store_true",
                  help="Show version and quit")
options, args = parser.parse_args()
if args:
    stop_err("Expected options (e.g. --input example.maf), not arguments")

input_maf = options.input
out_maf = options.maf
out_bam = options.bam
out_fasta = options.fasta
out_ace = options.ace
out_cstats = options.cstats

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
if options.version:
    print("%s, MIRA wrapper version %s" % (mira_convert_ver, WRAPPER_VER))
    sys.exit(0)

if not input_maf:
    stop_err("Input MIRA file is required")
elif not os.path.isfile(input_maf):
    stop_err("Missing input MIRA file: %r" % input_maf)

if not (out_maf or out_bam or out_fasta or out_ace or out_cstats):
    stop_err("No output requested")


def check_min_int(value, name):
    try:
        i = int(value)
    except:
        stop_err("Bad %s setting, %r" % (name, value))
    if i < 0:
        stop_err("Negative %s setting, %r" % (name, value))
    return i

min_length = check_min_int(options.min_length, "minimum length")
min_cover = check_min_int(options.min_cover, "minimum cover")
min_reads = check_min_int(options.min_reads, "minimum reads")

#TODO - Run MIRA in /tmp or a configurable directory?
#Currently Galaxy puts us somewhere safe like:
#/opt/galaxy-dist/database/job_working_directory/846/
temp = "."


cmd_list = [mira_convert]
if min_length:
    cmd_list.extend(["-x", str(min_length)])
if min_cover:
    cmd_list.extend(["-y", str(min_cover)])
if min_reads:
    cmd_list.extend(["-z", str(min_reads)])
cmd_list.extend(["-f", "maf", input_maf, os.path.join(temp, "converted")])
if out_maf:
    cmd_list.append("maf")
if out_bam:
    cmd_list.append("samnbb")
    if not out_fasta:
        #Need this for samtools depad
        out_fasta = os.path.join(temp, "depadded.fasta")
if out_fasta:
    cmd_list.append("fasta")
if out_ace:
    cmd_list.append("ace")
if out_cstats:
    cmd_list.append("cstats")
run(cmd_list)

def collect(old, new):
    if not os.path.isfile(old):
        stop_err("Missing expected output file %s" % old)
    shutil.move(old, new)

if out_maf:
    collect(os.path.join(temp, "converted.maf"), out_maf)
if out_fasta:
    #Can we look at the MAF file to see if there are multiple strains?
    old = os.path.join(temp, "converted_AllStrains.unpadded.fasta")
    if os.path.isfile(old):
        collect(old, out_fasta)
    else:
        #Might the output be filtered down to zero contigs?
        old = os.path.join(temp, "converted.fasta")
        if not os.path.isfile(old):
            stop_err("Missing expected output FASTA file")
        elif os.path.getsize(old) == 0:
            print("Warning - no contigs (harsh filters?)")
            collect(old, out_fasta)
        else:
            stop_err("Missing expected output FASTA file (only generic file present)")
if out_ace:
    collect(os.path.join(temp, "converted.maf"), out_ace)
if out_cstats:
    collect(os.path.join(temp, "converted_info_contigstats.txt"), out_cstats)

if out_bam:
    assert os.path.isfile(out_fasta)
    old = os.path.join(temp, "converted.samnbb")
    if not os.path.isfile(old):
        old = os.path.join(temp, "converted.sam")
    if not os.path.isfile(old):
        stop_err("Missing expected intermediate file %s" % old)
    h = BytesIO()
    msg = depad(out_fasta, old, out_bam, h)
    if msg:
        print(msg)
        print(h.getvalue())
        h.close()
        sys.exit(1)
    h.close()
    if out_fasta == os.path.join(temp, "depadded.fasta"):
        #Not asked for by Galaxy, no longer needed
        os.remove(out_fasta)

if min_length or min_cover or min_reads:
    print("Filtered.")
else:
    print("Converted.")
