#!/usr/bin/env python
"""A simple wrapper script to call MIRA and collect its output.
"""
import os
import sys
import subprocess
import shutil
import time

WRAPPER_VER = "0.0.1" #Keep in sync with the XML file

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
    return ver.split("\n", 1)[0]


os.environ["PATH"] = "/mnt/galaxy/downloads/mira_4.0rc3_linux-gnu_x86_64_static/bin/:%s" % os.environ["PATH"]
mira_binary = "mira"
mira_ver = get_version(mira_binary)
if not mira_ver.strip().startswith("4.0"):
    stop_err("This wrapper is for MIRA V4.0, not:\n%s" % mira_ver)
if "-v" in sys.argv:
    print "MIRA wrapper version %s," % WRAPPER_VER
    print mira_ver
    sys.exit(0)


def log_manifest(manifest):
    """Write the manifest file to stderr."""
    sys.stderr.write("\n%s\nManifest file\n%s\n" % ("="*60, "="*60))
    with open(manifest) as h:
        for line in h:
            sys.stderr.write(line)
    sys.stderr.write("\n%s\nEnd of manifest\n%s\n" % ("="*60, "="*60))


def collect_output(temp, name, handle):
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
    if not os.path.isfile(old_fasta):
        #Mapping (currently StrainX versus reference)
        old_fasta = "%s/%s_out_StrainX.unpadded.fasta" % (f, name)
    if not os.path.isfile(old_fasta):
        #Triggered extractLargeContigs.sh?
        old_fasta = "%s/%s_LargeContigs_out.fasta" % (f, name)

    missing = False
    for old, new in [(old_maf, out_maf),
                     (old_fasta, out_fasta)]:
        if not os.path.isfile(old):
            missing = True
        else:
            handle.write("Capturing %s\n" % old)
            shutil.move(old, new)
    if missing:
        log_manifest(manifest)
        sys.stderr.write("Contents of %r:\n" % f)
        for filename in sorted(os.listdir(f)):
            sys.stderr.write("%s\n" % filename)

def clean_up(temp, name):
    folder = "%s/%s_assembly" % (temp, name)
    if os.path.isdir(folder):
        shutil.rmtree(folder)

#TODO - Run MIRA in /tmp or a configurable directory?
#Currently Galaxy puts us somewhere safe like:
#/opt/galaxy-dist/database/job_working_directory/846/
temp = "."
#name, out_fasta, out_qual, out_ace, out_caf, out_wig, out_log = sys.argv[1:8]
name = "MIRA"
manifest, out_maf, out_fasta, out_log = sys.argv[1:5]

start_time = time.time()
#cmd_list =sys.argv[8:]
cmd_list = [mira_binary, manifest]
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

handle = open(out_log, "w")
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
print "MIRA took %0.2f hours" % (run_time / 3600.0)
if return_code:
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
    handle.write("===================== Extract Large Contigs failed? ==========================\n")
    e = open("MIRA_assembly/MIRA_d_results/ec.log", "rU")
    for line in e:
        handle.write(line)
    e.close()
    handle.write("============================ (end of ec.log) =================================\n")
    handle.flush()

#print "Collecting output..."
collect_output(temp, name, handle)

#print "Cleaning up..."
clean_up(temp, name)

handle.close()
print "Done"
