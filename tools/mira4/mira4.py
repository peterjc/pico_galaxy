#!/usr/bin/env python
"""A simple wrapper script to call MIRA and collect its output.
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

WRAPPER_VER = "0.0.4" #Keep in sync with the XML file

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

$ python mira4.py ...

This will run the MIRA binary and collect its output files as directed.
""" % WRAPPER_VER
parser = OptionParser(usage=usage)
parser.add_option("-m", "--manifest", dest="manifest",
                  default=None, metavar="FILE",
                  help="MIRA manifest filename")
parser.add_option("--maf", dest="maf",
                  default="-", metavar="FILE",
                  help="MIRA MAF output filename")
parser.add_option("--bam", dest="bam",
                  default="-", metavar="FILE",
                  help="Unpadded BAM output filename")
parser.add_option("--fasta", dest="fasta",
                  default="-", metavar="FILE",
                  help="Unpadded FASTA output filename")
parser.add_option("--log", dest="log",
                  default="-", metavar="FILE",
                  help="MIRA logging output filename")
options, args = parser.parse_args()
manifest = options.manifest
out_maf = options.maf
out_bam = options.bam
out_fasta = options.fasta
out_log = options.log

try:
    mira_path = os.environ["MIRA4"]
except KeyError:
    stop_err("Environment variable $MIRA4 not set")
mira_binary = os.path.join(mira_path, "mira")
if not os.path.isfile(mira_binary):
    stop_err("Missing mira under $MIRA4, %r\nFolder contained: %s"
             % (mira_binary, ", ".join(os.listdir(mira_path))))
mira_convert = os.path.join(mira_path, "miraconvert")
if not os.path.isfile(mira_convert):
    stop_err("Missing miraconvert under $MIRA4, %r\nFolder contained: %s"
             % (mira_convert, ", ".join(os.listdir(mira_path))))

mira_ver = get_version(mira_binary)
if not mira_ver.strip().startswith("4.0"):
    stop_err("This wrapper is for MIRA V4.0, not:\n%s\n%s" % (mira_ver, mira_binary))
mira_convert_ver = get_version(mira_convert)
if not mira_convert_ver.strip().startswith("4.0"):
    stop_err("This wrapper is for MIRA V4.0, not:\n%s\n%s" % (mira_ver, mira_convert))
if "-v" in sys.argv or "--version" in sys.argv:
    print "%s, MIRA wrapper version %s" % (mira_ver, WRAPPER_VER)
    if mira_ver != mira_convert_ver:
        print "WARNING: miraconvert %s" % mira_convert_ver
    sys.exit(0)

if not manifest:
    stop_err("Manifest is required")
elif not os.path.isfile(manifest):
    stop_err("Missing input MIRA manifest file: %r" % manifest)


try:
    threads = int(os.environ.get("GALAXY_SLOTS", "1"))
except ValueError:
    threads = 1
assert 1 <= threads, threads


def override_temp(manifest):
    """Override ``-DI:trt=/tmp`` in manifest with environment variable.

    Currently MIRA 4 does not allow envronment variables like ``$TMP``
    inside the manifest, which is a problem if you need to override
    the default at run time.

    The tool XML will ``/tmp`` and we replace that here with
    ``tempfile.gettempdir()`` which will respect $TMPDIR, $TEMP, $TMP
    as explained in the Python standard library documentation:
    http://docs.python.org/2/library/tempfile.html#tempfile.tempdir

    By default MIRA 4 would write its temporary files within the output
    folder, which is a problem if that is a network drive.
    """
    handle = open(manifest, "r")
    text = handle.read()
    handle.close()

    #At time of writing, this is at the end of a file,
    #but could be followed by a space in future...
    text = text.replace("-DI:trt=/tmp", "-DI:trt=" + tempfile.gettempdir())

    #Want to try to ensure this gets written to disk before MIRA attempts
    #to open it - any networked file system may impose a delay...
    handle = open(manifest, "w")
    handle.write(text)
    handle.flush()
    os.fsync(handle.fileno())
    handle.close()


def log_manifest(manifest):
    """Write the manifest file to stderr."""
    sys.stderr.write("\n%s\nManifest file\n%s\n" % ("="*60, "="*60))
    with open(manifest) as h:
        for line in h:
            sys.stderr.write(line)
    sys.stderr.write("\n%s\nEnd of manifest\n%s\n" % ("="*60, "="*60))


def collect_output(temp, name, handle):
    """Moves files to the output filenames (global variables)."""
    n3 = (temp, name, name, name)
    f = "%s/%s_assembly/%s_d_results" % (temp, name, name)
    if not os.path.isdir(f):
        log_manifest(manifest)
        stop_err("Missing output folder")
    if not os.listdir(f):
        log_manifest(manifest)
        stop_err("Empty output folder")
    missing = []

    old_maf = "%s/%s_out.maf" % (f, name)
    if not os.path.isfile(old_maf):
        #Triggered extractLargeContigs.sh?
        old_maf = "%s/%s_LargeContigs_out.maf" % (f, name)

    #De novo or single strain mapping,
    old_fasta = "%s/%s_out.unpadded.fasta" % (f, name)
    ref_fasta = "%s/%s_out.padded.fasta" % (f, name)
    if not os.path.isfile(old_fasta):
        #Mapping (StrainX versus reference) or de novo
        old_fasta = "%s/%s_out_StrainX.unpadded.fasta" % (f, name)
        ref_fasta = "%s/%s_out_StrainX.padded.fasta" % (f, name)
    if not os.path.isfile(old_fasta):
        old_fasta = "%s/%s_out_ReferenceStrain.unpadded.fasta" % (f, name)
        ref_fasta = "%s/%s_out_ReferenceStrain.padded.fasta" % (f, name)
        

    missing = False
    for old, new in [(old_maf, out_maf),
                     (old_fasta, out_fasta)]:
        if not os.path.isfile(old):
            missing = True
        elif not new or new == "-":
            handle.write("Ignoring %s\n" % old)
        else:
            handle.write("Capturing %s\n" % old)
            shutil.move(old, new)
    if missing:
        log_manifest(manifest)
        sys.stderr.write("Contents of %r:\n" % f)
        for filename in sorted(os.listdir(f)):
            sys.stderr.write("%s\n" % filename)

    #For mapping mode, probably most people would expect a BAM file
    #using the reference FASTA file...
    if out_bam and out_bam != "-":
        msg = make_bam(mira_convert, out_maf, ref_fasta, out_bam, handle)
        if msg:
            stop_err(msg)

def clean_up(temp, name):
    folder = "%s/%s_assembly" % (temp, name)
    if os.path.isdir(folder):
        shutil.rmtree(folder)

#TODO - Run MIRA in /tmp or a configurable directory?
#Currently Galaxy puts us somewhere safe like:
#/opt/galaxy-dist/database/job_working_directory/846/
temp = "."

name = "MIRA"

override_temp(manifest)

start_time = time.time()
cmd_list = [mira_binary, "-t", str(threads), manifest]
cmd = " ".join(cmd_list)

assert os.path.isdir(temp)
d = "%s_assembly" % name
assert not os.path.isdir(d), "Path %s already exists" % d
try:
    #Check path access
    os.mkdir(d)
except Exception, err:
    log_manifest(manifest)
    sys.stderr.write("Error making directory %s\n%s" % (d, err))
    sys.exit(1)

#print os.path.abspath(".")
#print cmd

if out_log and out_log != "-":
    handle = open(out_log, "w")
else:
    handle = open(os.devnull, "w")
handle.write("======================== MIRA manifest (instructions) ========================\n")
m = open(manifest, "rU")
for line in m:
    handle.write(line)
m.close()
del m
handle.write("\n")
handle.write("============================ Starting MIRA now ===============================\n")
handle.flush()
try:
    #Run MIRA
    child = subprocess.Popen(cmd_list,
                             stdout=handle,
                             stderr=subprocess.STDOUT)
except Exception, err:
    log_manifest(manifest)
    sys.stderr.write("Error invoking command:\n%s\n\n%s\n" % (cmd, err))
    #TODO - call clean up?
    handle.write("Error invoking command:\n%s\n\n%s\n" % (cmd, err))
    handle.close()
    sys.exit(1)
#Use .communicate as can get deadlocks with .wait(),
stdout, stderr = child.communicate()
assert not stdout and not stderr #Should be empty as sent to handle
run_time = time.time() - start_time
return_code = child.returncode
handle.write("\n")
handle.write("============================ MIRA has finished ===============================\n")
handle.write("MIRA took %0.2f hours\n" % (run_time / 3600.0))
if return_code:
    print "MIRA took %0.2f hours" % (run_time / 3600.0)
    handle.write("Return error code %i from command:\n" % return_code)
    handle.write(cmd + "\n")
    handle.close()
    clean_up(temp, name)
    log_manifest(manifest)
    stop_err("Return error code %i from command:\n%s" % (return_code, cmd),
             return_code)
handle.flush()

if os.path.isfile("MIRA_assembly/MIRA_d_results/ec.log"):
    handle.write("\n")
    handle.write("====================== Extract Large Contigs failed ==========================\n")
    e = open("MIRA_assembly/MIRA_d_results/ec.log", "rU")
    for line in e:
        handle.write(line)
    e.close()
    handle.write("============================ (end of ec.log) =================================\n")
    handle.flush()

#print "Collecting output..."
start_time = time.time()
collect_output(temp, name, handle)
collect_time = time.time() - start_time
handle.write("MIRA took %0.2f hours; collecting output %0.2f minutes\n" % (run_time / 3600.0, collect_time / 60.0))
print("MIRA took %0.2f hours; collecting output %0.2f minutes\n" % (run_time / 3600.0, collect_time / 60.0))

if os.path.isfile("MIRA_assembly/MIRA_d_results/ec.log"):
    #Treat as an error, but doing this AFTER collect_output
    sys.stderr.write("Extract Large Contigs failed\n")
    handle.write("Extract Large Contigs failed\n")
    handle.close()
    sys.exit(1)

#print "Cleaning up..."
clean_up(temp, name)

handle.write("\nDone\n")
handle.close()
print("Done")
