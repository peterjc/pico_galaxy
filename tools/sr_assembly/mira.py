#!/usr/bin/env python
"""A simple wrapper script to call MIRA and collect its output
"""
import os
import sys
import subprocess
import shutil
import time

def stop_err(msg, err=1):
    sys.stderr.write(msg+"\n")
    sys.exit(err)


name, out_fasta, out_ace, out_wig, out_log = sys.argv[1:6]
start_time = time.time()
try:
    cmd = " ".join(sys.argv[6:])
    child = subprocess.Popen(sys.argv[6:],
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
handle.write("Mira took %0.2f minutes" % (run_time / 60.0))
print "Mira took %0.2f minutes" % (run_time / 60.0)
if return_code:
    handle.write("Return error code %i from command:\n" % return_code)
    handle.write(cmd + "\n")
    handle.close()
    stop_err("Return error code %i from command:\n%s" % (return_code, cmd),
             return_code)
handle.close()

temp = "."
n3 = (temp, name, name, name)
for old, new in [("%s/%s_assembly/%s_d_results/%s_out.unpadded.fasta" % n3, out_fasta),
                 ("%s/%s_assembly/%s_d_results/%s_out.wig" % n3, out_wig),
                 ("%s/%s_assembly/%s_d_results/%s_out.ace" % n3, out_ace)]:
    if not os.path.isfile(old):
        stop_err("Missing %s output file" % os.path.splitext(old)[-1])
    else:
        shutil.move(old, new)
shutil.rmtree("%s/%s_assembly" % (temp, name))
print "Done"
