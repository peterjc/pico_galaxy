#!/usr/bin/env python
"""Wrapper for TMHMM v2.0 for use in Galaxy.

This script takes exactly three command line arguments - number of threads,
an input protein FASTA filename, and an output tabular filename. It then
calls the standalone TMHMM v2.0 program (not the webservice) requesting
the short output (one line per protein).

The first major feature is cleaning up the tabular output. The short form raw
output from TMHMM v2.0 looks like this (six columns tab separated):

 gi|2781234|pdb|1JLY|B	len=304 ExpAA=0.01	First60=0.00	PredHel=0	Topology=o
 gi|4959044|gb|AAD34209.1|AF069992_1	len=600	ExpAA=0.00	First60=0.00	PredHel=0	Topology=o
 gi|671626|emb|CAA85685.1|	len=473 ExpAA=0.19	First60=0.00 PredHel=0	Topology=o
 gi|3298468|dbj|BAA31520.1|	len=107	ExpAA=59.37	First60=31.17	PredHel=3	Topology=o23-45i52-74o89-106i

If there are any additional 'comment' lines starting with the hash (#)
character these are ignored by this script.

In order to make it easier to use in Galaxy, this wrapper script simplifies
this to remove the redundant tags, and instead adds a comment line at the
top with the column names:

 #ID	len	ExpAA	First60	PredHel	Topology 
 gi|2781234|pdb|1JLY|B	304	0.01	60	0.00	0	o
 gi|4959044|gb|AAD34209.1|AF069992_1	600	0.00	0	0.00	0	o
 gi|671626|emb|CAA85685.1|	473	0.19	0.00	0	o
 gi|3298468|dbj|BAA31520.1|	107	59.37	31.17	3	o23-45i52-74o89-106i

The second major potential feature is taking advantage of multiple cores
(since TMHMM v2.0 itself is single threaded) by dividing the input FASTA file
into chunks and running multiple copies of TMHMM in parallel. I would normally
use Python's multiprocessing library in this situation but it requires at
least Python 2.6 and at the time of writing Galaxy still supports Python 2.4.

Note that this is somewhat redundant with job-splitting available in Galaxy
itself (see the SignalP XML file for settings).

Also tmhmm2 can fail without returning an error code, for example if run on a
64 bit machine with only the 32 bit binaries installed. This script will spot
when there is no output from tmhmm2, and raise an error.
"""
import sys
import os
import tempfile
from seq_analysis_utils import stop_err, split_fasta, run_jobs, thread_count

FASTA_CHUNK = 500

if len(sys.argv) != 4:
   stop_err("Require three arguments, number of threads (int), input protein FASTA file & output tabular file")

num_threads = thread_count(sys.argv[1], default=4)
fasta_file = sys.argv[2]
tabular_file = sys.argv[3]

tmp_dir = tempfile.mkdtemp()

def clean_tabular(raw_handle, out_handle):
    """Clean up tabular TMHMM output, returns output line count."""
    count = 0
    for line in raw_handle:
        if not line.strip() or line.startswith("#"):
            #Ignore any blank lines or comment lines
            continue
        parts = line.rstrip("\r\n").split("\t")
        try:
            identifier, length, expAA, first60, predhel, topology = parts
        except:
            assert len(parts)!=6
            stop_err("Bad line: %r" % line)
        assert length.startswith("len="), line
        length = length[4:]
        assert expAA.startswith("ExpAA="), line
        expAA = expAA[6:]
        assert first60.startswith("First60="), line
        first60 = first60[8:]
        assert predhel.startswith("PredHel="), line
        predhel = predhel[8:]
        assert topology.startswith("Topology="), line
        topology = topology[9:]
        out_handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" \
                   % (identifier, length, expAA, first60, predhel, topology))
        count += 1
    return count

#Note that if the input FASTA file contains no sequences,
#split_fasta returns an empty list (i.e. zero temp files).
fasta_files = split_fasta(fasta_file, os.path.join(tmp_dir, "tmhmm"), FASTA_CHUNK)
temp_files = [f+".out" for f in fasta_files]
jobs = ["tmhmm -short %s > %s" % (fasta, temp)
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
results = run_jobs(jobs, num_threads)
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
out_handle.write("#ID\tlen\tExpAA\tFirst60\tPredHel\tTopology\n")
for temp in temp_files:
    data_handle = open(temp)
    count = clean_tabular(data_handle, out_handle)
    data_handle.close()
    if not count:
        clean_up(fasta_files + temp_files)
        stop_err("No output from tmhmm2")
out_handle.close()

clean_up(fasta_files + temp_files)
