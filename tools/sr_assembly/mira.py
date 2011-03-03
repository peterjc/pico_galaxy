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

def tcs_to_tabular(old, new):
    in_handle = open(old, "rU")
    out_handle = open(new, "w")
    assert in_handle.readline() == "#TCS V1.0\n"
    assert in_handle.readline() == "#\n"
    assert in_handle.readline() == "# contig name          padPos  upadPos | B  Q | tcov covA covC covG covT cov* | qA qC qG qT q* |  S | Tags\n"
    assert in_handle.readline() == "#\n"
    out_handle.write("#%s\n" % "\t".join(["contig", "pasPos", "upadPos", "B", "Q",
                                         "tcov", "covA", "covC", "covG", "covT", "cov*",
                                         "qA", "qC", "qG", "qT", "q*", "S", "Tags"]))
    for line in in_handle:
        parts = line.rstrip("\n").split(None,22)
        assert parts[3] == parts[6] == parts[13] == parts[19] == parts[21] == "|"
        wanted = parts[:3] + parts[4:6]+parts[7:13]+parts[14:19]+parts[20:21]+parts[22:]
        out_handle.write("%s\n" % "\t".join(wanted))
    out_handle.close()
    in_handle.close()

def collect_output(temp, name):
    n3 = (temp, name, name, name)
    tcs_to_tabular("%s/%s_assembly/%s_d_results/%s_out.tcs" % n3, out_tcs)
    for old, new in [("%s/%s_assembly/%s_d_results/%s_out.unpadded.fasta" % n3, out_fasta),
                     ("%s/%s_assembly/%s_d_results/%s_out.unpadded.fasta.qual" % n3, out_qual),
                     ("%s/%s_assembly/%s_d_results/%s_out.wig" % n3, out_wig),
                     ("%s/%s_assembly/%s_d_results/%s_out.caf" % n3, out_caf),
                     ("%s/%s_assembly/%s_d_results/%s_out.ace" % n3, out_ace)]:
        if not os.path.isfile(old):
            stop_err("Missing %s output file" % os.path.splitext(old)[-1])
        else:
            shutil.move(old, new)

def clean_up(temp, name):
    folder = "%s/%s_assembly" % (temp, name)
    if os.path.isdir(folder):
        shutil.rmtree(folder)

#TODO - Run MIRA in /tmp or a configurable directory?
temp = "."
name, out_fasta, out_qual, out_tcs, out_ace, out_caf, out_wig, out_log = sys.argv[1:9]
start_time = time.time()
try:
    cmd = " ".join(sys.argv[9:])
    child = subprocess.Popen(sys.argv[9:],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
except Exception, err:
    sys.stderr.write("Error invoking command:\n%s\n\n%s\n" % (cmd, err))
    sys.exit(1)
#Use .communicate as can get deadlocks with .wait(),
stdout, stderr = child.communicate()
run_time = time.time() - start_time
return_code = child.returncode
handle = open(out_log, "w")
handle.write(stdout)
handle.write("MIRA took %0.2f minutes" % (run_time / 60.0))
print "MIRA took %0.2f minutes" % (run_time / 60.0)
if return_code:
    handle.write("Return error code %i from command:\n" % return_code)
    handle.write(cmd + "\n")
    handle.close()
    clean_up(temp, name)
    stop_err("Return error code %i from command:\n%s" % (return_code, cmd),
             return_code)
handle.close()

collect_output(temp, name)
clean_up(temp, name)
print "Done"
