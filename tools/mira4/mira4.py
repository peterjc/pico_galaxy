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


def massage_symlinks(manifest):
    """Create FASTQ aliases and edit the manifest to use them.

    Short term measure for MIRA 4.0RC2 which depends on data file
    extensions to decide the file format, and doesn't like *.dat
    as used in Galaxy.
    """
    base = os.path.split(manifest)[0]
    with open(manifest) as h:
        lines = h.readlines()
    f = 0
    for i, line in enumerate(lines):
         if not line.startswith("data ="):
             continue
         #Assumes no spaces in filename, would they have to be escaped?
         new_line = "data ="
         for filename in line[6:].strip().split():
             if not filename:
                 continue
             assert os.path.isfile(filename), filename
             f += 1
             alias = os.path.join(base, "input%i.fastq" % f)
             new_line += " " + alias
             cmd = "ln -s %s %s" % (filename, alias)
             if os.system(cmd):
                 stop_err("Problem creating FASTQ alias:\n%s" % cmd)
         lines[i] = new_line + "\n"
    with open(manifest, "w") as h:
        for line in lines:
            #sys.stderr.write(line)
            h.write(line)
    return True


def collect_output(temp, name):
    n3 = (temp, name, name, name)
    f = "%s/%s_assembly/%s_d_results" % (temp, name, name)
    if not os.path.isdir(f):
        log_manifest(manifest)
        stop_err("Missing output folder")
    if not os.listdir(f):
        log_manifest(manifest)
        stop_err("Empty output folder")
    missing = []
    for old, new in [("%s/%s_out.maf" % (f, name), out_maf),
                     ("%s/%s_out.unpadded.fasta" % (f, name), out_fasta)]:
        if not os.path.isfile(old):
            missing.append(os.path.splitext(old)[-1])
        else:
            shutil.move(old, new)
    if missing:
        log_manifest(manifest)
        sys.stderr.write("Contents of %r: %r\n" % (f, os.listdir(f)))
        stop_err("Missing output files: %s" % ", ".join(missing))

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

#Hack until MIRA v4 lets us specify file format explicitly,
massage_symlinks(manifest)

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
handle.write("\n\nMIRA took %0.2f minutes\n" % (run_time / 60.0))
print "MIRA took %0.2f minutes" % (run_time / 60.0)
if return_code:
    handle.write("Return error code %i from command:\n" % return_code)
    handle.write(cmd + "\n")
    handle.close()
    clean_up(temp, name)
    log_manifest(manifest)
    stop_err("Return error code %i from command:\n%s" % (return_code, cmd),
             return_code)
handle.close()

#print "Collecting output..."
collect_output(temp, name)

#print "Cleaning up..."
clean_up(temp, name)

print "Done"
