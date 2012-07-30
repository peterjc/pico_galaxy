#!/usr/bin/env python
"""Wrapper for Promoter 2.0 for use in Galaxy.

This script takes exactly three command line arguments:
 * number of threads
 * an input DNA FASTA filename
 * output tabular filename.

It calls the Promoter 2.0 binary (e.g. .../promoter-2.0/bin/promoter_Linux,
bypassing the Perl wrapper script 'promoter' which imposes a significant
performace overhead for no benefit here (we don't need HTML output for
example).

The main feature is this Python wrapper script parsers the bespoke
tabular output from Promoter 2.0 and reformats it into a Galaxy friendly
tab separated table.

Additionally, in order to take advantage of multiple cores the input FASTA
file is broken into chunks and multiple copies of promoter run at once.
This can be used in combination with the job-splitting available in Galaxy.

Note that rewriting the FASTA input file allows us to avoid a bug in
promoter 2 with long descriptions in the FASTA header line (over 200
characters) which produces stray fragements of the description in the
output file, making parsing non-trivial.

TODO - Automatically extract the sequence containing a promoter prediction?
"""
import sys
import os
import commands
import tempfile
from seq_analysis_utils import stop_err, split_fasta, run_jobs

FASTA_CHUNK = 500

if len(sys.argv) != 4:
    stop_err("Require three arguments, number of threads (int), input DNA FASTA file & output tabular file. "
             "Got %i arguments." % (len(sys.argv)-1))
try:
    num_threads = int(sys.argv[1])
except:
    num_threads = 1 #Default, e.g. used "$NSLOTS" and environment variable not defined
if num_threads < 1:
    stop_err("Threads argument %s is not a positive integer" % sys.argv[1])

fasta_file = os.path.abspath(sys.argv[2])
tabular_file = os.path.abspath(sys.argv[3])

tmp_dir = tempfile.mkdtemp()

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

def make_tabular(raw_handle, out_handle):
    """Parse text output into tabular, return query count."""
    identifier = None
    queries = 0
    #out.write("#Identifier\tDescription\tPosition\tScore\tLikelihood\n")
    for line in raw_handle:
        #print repr(line)
        if not line.strip() or line == "Promoter prediction:\n":
            pass
        elif line[0] != " ":
            identifier = line.strip().replace("\t", " ").split(None,1)[0]
            queries += 1
        elif line == "  No promoter predicted\n":
            #End of a record
            identifier = None
        elif line == "  Position  Score  Likelihood\n":
            assert identifier
        else:
            try:
                position, score, likelihood = line.strip().split(None,2)
            except ValueError:
                print "WARNING: Problem with line: %r" % line
                continue
                #stop_err("ERROR: Problem with line: %r" % line)
            if likelihood not in ["ignored",
                                  "Marginal prediction",
                                  "Medium likely prediction",
                                  "Highly likely prediction"]:
                stop_err("ERROR: Problem with line: %r" % line)
            out_handle.write("%s\t%s\t%s\t%s\n" % (identifier, position, score, likelihood))
    #out.close()
    return queries
    
working_dir, bin = get_path_and_binary()

if not os.path.isfile(fasta_file):
   stop_err("ERROR: Missing input FASTA file %r" % fasta_file)

#Note that if the input FASTA file contains no sequences,
#split_fasta returns an empty list (i.e. zero temp files).
#We deliberately omit the FASTA descriptions to avoid a
#bug in promoter2 with descriptions over 200 characters.
fasta_files = split_fasta(fasta_file, os.path.join(tmp_dir, "promoter"), FASTA_CHUNK, keep_descr=False)
temp_files = [f+".out" for f in fasta_files]
jobs = ["%s %s > %s" % (bin, fasta, temp)
        for fasta, temp in zip(fasta_files, temp_files)]

def clean_up(file_list):
    for f in file_list:
        if os.path.isfile(f):
            os.remove(f)
    try:
        os.rmdir(tmp_dir)
    except:
        pass

if len(jobs) > 1 and num_threads > 1:
    #A small "info" message for Galaxy to show the user.
    print "Using %i threads for %i tasks" % (min(num_threads, len(jobs)), len(jobs))
cur_dir = os.path.abspath(os.curdir)
os.chdir(working_dir)
results = run_jobs(jobs, num_threads)
os.chdir(cur_dir)
for fasta, temp, cmd in zip(fasta_files, temp_files, jobs):
    error_level = results[cmd]
    if error_level:
        try:
            output = open(temp).readline()
        except IOError:
            output = ""
        clean_up(fasta_files + temp_files)
        stop_err("One or more tasks failed, e.g. %i from %r gave:\n%s" % (error_level, cmd, output),
                 error_level)

del results
del jobs

out_handle = open(tabular_file, "w")
out_handle.write("#Identifier\tDescription\tPosition\tScore\tLikelihood\n")
queries = 0
for temp in temp_files:
    data_handle = open(temp)
    count = make_tabular(data_handle, out_handle)
    data_handle.close()
    if not count:
        clean_up(fasta_files + temp_files)
        stop_err("No output from promoter2")
    out_handle.close()
    queries += count

clean_up(fasta_files + temp_files)
print "Results for %i queries" % queries
