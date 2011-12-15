#!/usr/bin/env python
"""A simple wrapper script to call MIRA and collect its output.
"""
import os
import sys
import subprocess
import shutil
import time

def stop_err(msg, err=1):
    sys.stderr.write(msg+"\n")
    sys.exit(err)

def collect_output(temp, name):
    n3 = (temp, name, name, name)
    f = "%s/%s_assembly/%s_d_results" % (temp, name, name)
    if not os.path.isdir(f):
        stop_err("Missing output folder")
    if not os.listdir(f):
        stop_err("Empty output folder")
    for old, new in [("%s/%s_out.unpadded.fasta" % (f, name), out_fasta),
                     ("%s/%s_out.unpadded.fasta.qual" % (f, name), out_qual),
                     ("%s/%s_out.wig" % (f, name), out_wig),
                     ("%s/%s_out.caf" % (f, name), out_caf),
                     ("%s/%s_out.ace" % (f, name), out_ace)]:
        if not os.path.isfile(old):
            stop_err("Missing %s output file" % os.path.splitext(old)[-1])
        else:
            shutil.move(old, new)

def clean_up(temp, name):
    folder = "%s/%s_assembly" % (temp, name)
    if os.path.isdir(folder):
        shutil.rmtree(folder)

#TODO - Run MIRA in /tmp or a configurable directory?
#Currently Galaxy puts us somewhere safe like:
#/opt/galaxy-dist/database/job_working_directory/846/
temp = "."
name, out_fasta, out_qual, out_ace, out_caf, out_wig, out_log = sys.argv[1:8]

start_time = time.time()
cmd_list =sys.argv[8:]
cmd = " ".join(cmd_list)

assert os.path.isdir(temp)
d = "%s_assembly" % name
assert not os.path.isdir(d), "Path %s already exists" % d
try:
    #Check path access
    os.mkdir(d)
except Exception, err:
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
    stop_err("Return error code %i from command:\n%s" % (return_code, cmd),
             return_code)
handle.close()

#print "Collecting output..."
collect_output(temp, name)

#print "Cleaning up..."
clean_up(temp, name)

print "Done"
