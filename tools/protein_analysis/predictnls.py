#!/usr/bin/env python
"""Wrapper for predictNLS v1.3 for use in Galaxy.

This script takes exactly three command line arguments:
 * number of threads to use (integer)
 * an input protein FASTA filename
 * output tabular filename.

It then calls the standalone predictNLS v1.3 program predictNLS (not the
webservice), and coverts the text output from something like this:

Results of Nuclear Localization Signal Prediction(NLS)
The NLS server can be directly accessed at: http://cubic.bioc.columbia.edu/predictNLS/
For help on interpretation of results visit the predictNLS help page:
http://cubic.bioc.columbia.edu/predictNLS/doc/help.html
-----------------------------------------------------------------
----------------------------------------------------------------------
Input sequence Id: >MYOD1_HUMAN
Sequence Length: 319
----------------------------------------------------------------------
List of NLS's found in sequence
----------------------------------------------------------------------
|                NLS|Position in sequence|
|      CKRKTTNADRRKA|                 100|
|        KRKTTNADRRK|                 101|
|      RRKAATMRERRRL|                 109|
|       RRKAATMRERRR|                 109|
|        VNEAFETLKRC|                 124|
----------------------------------------------------------------------
...

In order to make it easier to use in Galaxy, this wrapper script reformats
this to use tab separators, with one line per NLS prediction:

#ID	Start	NLS
MYOD1_HUMAN	100	CKRKTTNADRRKA
MYOD1_HUMAN	101	KRKTTNADRRK
MYOD1_HUMAN	109	RRKAATMRERRRL
MYOD1_HUMAN	109	RRKAATMRERRR
MYOD1_HUMAN	124	VNEAFETLKRC

Additionally in order to take full advantage of multiple cores, multiple copies
of predictNLS are run in parallel. I would normally use Python's multiprocessing
library in this situation but it requires at least Python 2.6 and at the time of
writing Galaxy still supports Python 2.4.
"""
import sys
import os
from seq_analysis_utils import stop_err, split_fasta, run_jobs

FASTA_CHUNK = 500
exe = "predictNLS"


if len(sys.argv) != 4:
   stop_err("Require three arguments: threads, input protein FASTA file & output tabular file")

try:
   num_threads = int(sys.argv[1])
except:
   num_threads = 0
if num_threads < 1:
   stop_err("Threads argument %s is not a positive integer" % sys.argv[1])

fasta_file = sys.argv[2]

tabular_file = sys.argv[3]


def clean_tabular(raw_handle, out_handle):
    """Clean up predictNLS output to make it tabular.

    This should cope with concatenated output, or the raw single
    query output direct from predictNLS.
    """
    id = None
    for line in raw_handle:
        if line.startswith("Input sequence Id:"):
            id = line[19:].strip()
            if id.startswith(">"):
                #predictNLS seems to leave the FASTA > sign in
                id = id[1:]
        elif line.startswith("|") and line.count("|")==3:
            assert id is not None
            nls, start = [s.strip() for s in line[1:].split("|")[0:2]]
            if nls != "NLS" and start != "Position in sequence":
                out_handle.write("%s\t%s\t%s\n" % (id, start, nls))


#Due to the way predictNLS works, require one sequence per FASTA (!!!)
fasta_files = split_fasta(fasta_file, tabular_file, n=1)
temp_files = [f+".out" for f in fasta_files]
assert len(fasta_files) == len(temp_files)
#Send stdout and stderr to /dev/null, its very noisy!
jobs = ["%s in=%s out=%s html=0 > /dev/null 2> /dev/null" \
        % (exe, os.path.abspath(fasta), os.path.abspath(temp))
        for (fasta, temp) in zip(fasta_files, temp_files)]
assert len(fasta_files) == len(temp_files) == len(jobs)

def clean_up(file_list):
    for f in file_list:
        if os.path.isfile(f):
            os.remove(f)

if len(jobs) > 1 and num_threads > 1:
    #A small "info" message for Galaxy to show the user.
    print "Using %i threads for %i tasks" % (min(num_threads, len(jobs)), len(jobs))
#Note predictNLS will create temp files of its own in the current directory,
results = run_jobs(jobs, num_threads, pause=0.2)
assert len(fasta_files) == len(temp_files) == len(jobs)
for fasta, temp, cmd in zip(fasta_files, temp_files, jobs):
    error_level = results[cmd]
    if not os.path.isfile(temp):
        clean_up(fasta_files)
        clean_up(temp_files)
        stop_err("One or more tasks failed, e.g. no output file (return %i) from %r" \
                 % (error_level, cmd))
    try:
        output = open(temp).readline()
    except IOError, err:
        clean_up(fasta_files)
        clean_up(temp_files)
        stop_err("One or more tasks failed, e.g. error opening output (return %i) from %r" \
                                 % (error_level, cmd))
    if not output.strip():
       clean_up(fasta_files)
       clean_up(temp_files)
       stop_err("One or more tasks failed, e.g. empty output (return %i) from %r" \
                % (error_level, cmd))
    if error_level:
        clean_up(fasta_files)
        clean_up(temp_files)
        stop_err("One or more tasks failed, e.g. %i from %r gave:\n%s" \
                 % (error_level, cmd, output),
                 error_level)
del results

out_handle = open(tabular_file, "w")
out_handle.write("#ID\tStart\tNLS\n")
for temp in temp_files:
    data_handle = open(temp)
    clean_tabular(data_handle, out_handle)
    data_handle.close()
out_handle.close()

clean_up(fasta_files)
clean_up(temp_files)
