#!/usr/bin/env python
"""Wrapper for Promoter 2.0 for use in Galaxy.

PRThis script takes exactly two command line arguments:
 * an input DNA FASTA filename
 * output tabular filename.

It calls the Promoter 2.0 binary (e.g. .../promoter-2.0/bin/promoter_Linux,
bypassing the Perl wrapper script 'promoter' which imposes a significant
performace overhead for no benefit here (we don't need HTML output for
example).

The main feature is this Python wrapper script parsers the bespoke
tabular output from Promoter 2.0 and reformats it into a Galaxy friendly
tab separated table.

Note that this wrapper does not (currently) take advantage of multiple
cores - we recommend using the job-splitting available in Galaxy for this.

TODO - Automatically extract the sequence containing a promoter prediction?
"""
import sys
import os
import subprocess
import commands
from seq_analysis_utils import stop_err

FASTA_CHUNK = 500
MAX_LEN = 6000 #Found by trial and error

if len(sys.argv) != 3:
   stop_err("Require two arguments: input DNA FASTA file & output tabular file. "
            "Got %i arguments." % (len(sys.argv)-1))

fasta_file = sys.argv[1]
tabular_file = sys.argv[2]

def get_path_and_binary():
    platform = commands.getoutput("uname") #e.g. Linux
    shell_script = commands.getoutput("which promoter")
    if not os.path.isfile(shell_script):
        stop_err("ERROR: Missing promoter executable shell script")
    path = None
    for line in open(shell_script):
        if line.startswith("setenv"): #could then be tab or space!
            parts = line.rstrip().split(None, 2)
            if parts[0] == "setenv" and parts[1] == "PROM":
                path = parts[2]
    if not path:
        stop_err("ERROR: Could not find promoter path (PROM) in %r" % shell_script)
    if not os.path.isdir(path):
        stop_error("ERROR: %r is not a directory" % path)
    bin = "%s/bin/promoter_%s" % (path, platform)
    if not os.path.isfile(bin):
        stop_err("ERROR: Missing promoter binary %r" % bin)
    return path, bin

def run_promoter(bin, fasta_file, tabular_file):
    child = subprocess.Popen([bin, fasta_file],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    identifier = None
    descr = None
    queries = 0
    out = open(tabular_file, "w")
    out.write("Identifier\tDescription\tPosition\tScore\tLikelihood\n")
    for line in child.stdout:
        #print repr(line)
        if not line.strip() or line == "Promoter prediction:\n":
            pass
        elif line[0] != " ":
            identifier = line.strip().replace("\t", " ")
            if " " in identifier:
                identifier, descr = identifier.split(None, 1)
            else:
                descr = ""
            queries += 1
        elif line == "  Position  Score  Likelihood\n":
            assert identifier
        else:
            position, score, likelihood = line.strip().split(None,2)
            out.write("%s\t%s\t%s\t%s\t%s\n" % (identifier, descr, position, score, likelihood))
    out.close()
    print "Results for %i sequences" % queries
    return_code = child.wait()
    if return_code:
       stop_err("ERROR: Return code %i from promoter" % return_code)
    
working_dir, bin = get_path_and_binary()

if not os.path.isfile(fasta_file):
   stop_err("ERROR: Missing input FASTA file %r" % fasta_file)

cur_dir = os.path.abspath(os.curdir)
os.chdir(working_dir)
run_promoter(bin, fasta_file, tabular_file)
os.chdir(cur_dir)
