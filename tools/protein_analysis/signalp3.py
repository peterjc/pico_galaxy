#!/usr/bin/env python
"""Wrapper for SignalP v3.0 for use in Galaxy.

This script takes exactly fives command line arguments:
 * the organism type (euk, gram+ or gram-)
 * length to truncate sequences to (integer)
 * number of threads to use (integer)
 * an input protein FASTA filename
 * output tabular filename.

It then calls the standalone SignalP v3.0 program (not the webservice)
requesting the short output (one line per protein) using both NN and HMM
for predictions.

First major feature is cleaning up the output. The raw output from SignalP
v3.0 looks like this (21 columns space separated):

# SignalP-NN euk predictions                                   	                # SignalP-HMM euk predictions
# name                Cmax  pos ?  Ymax  pos ?  Smax  pos ?  Smean ?  D     ? 	# name      !  Cmax  pos ?  Sprob ?
gi|2781234|pdb|1JLY|  0.061  17 N  0.043  17 N  0.199   1 N  0.067 N  0.055 N	gi|2781234|pdb|1JLY|B  Q  0.000  17 N  0.000 N  
gi|4959044|gb|AAD342  0.099 191 N  0.012  38 N  0.023  12 N  0.014 N  0.013 N	gi|4959044|gb|AAD34209.1|AF069992_1  Q  0.000   0 N  0.000 N  
gi|671626|emb|CAA856  0.139 381 N  0.020   8 N  0.121   4 N  0.067 N  0.044 N	gi|671626|emb|CAA85685.1|  Q  0.000   0 N  0.000 N  
gi|3298468|dbj|BAA31  0.208  24 N  0.184  38 N  0.980  32 Y  0.613 Y  0.398 N	gi|3298468|dbj|BAA31520.1|  Q  0.066  24 N  0.139 N

In order to make it easier to use in Galaxy, this wrapper script reformats
this to use tab separators. Also it removes the redundant truncated name
column, and assigns unique column names in the header:

#ID	NN_Cmax_score	NN_Cmax_pos	NN_Cmax_pred	NN_Ymax_score	NN_Ymax_pos	NN_Ymax_pred	NN_Smax_score	NN_Smax_pos	NN_Smax_pred	NN_Smean_score	NN_Smean_pred	NN_D_score	NN_D_pred	HMM_bang	HMM_Cmax_score	HMM_Cmax_pos	HMM_Cmax_pred	HMM_Sprob_score	HMM_Sprob_pred
gi|2781234|pdb|1JLY|B	0.061	17	N	0.043	17	N	0.199	1	N	0.067	N	0.055	N	Q	0.000	17	N	0.000	N
gi|4959044|gb|AAD34209.1|AF069992_1	0.099	191	N	0.012	38	N	0.023	12	N	0.014	N	0.013	N	Q	0.000	0	N	0.000	N
gi|671626|emb|CAA85685.1|	0.139	381	N	0.020	8	N	0.121	4	N	0.067	N	0.044	N	Q	0.000	0	N	0.000	N
gi|3298468|dbj|BAA31520.1|	0.208	24	N	0.184	38	N	0.980	32	Y	0.613	Y	0.398	N	Q	0.066	24	N	0.139	N

The second major feature is overcoming SignalP's built in limit of 4000
sequences by breaking up the input FASTA file into chunks. This also allows
us to pre-trim the sequences since SignalP only needs their starts.

The third major feature is taking advantage of multiple cores (since SignalP
v3.0 itself is single threaded) by using the individual FASTA input files to
run multiple copies of TMHMM in parallel. I would normally use Python's
multiprocessing library in this situation but it requires at least Python 2.6
and at the time of writing Galaxy still supports Python 2.4.
"""
import sys
import os
from seq_analysis_utils import stop_err, split_fasta, run_jobs

FASTA_CHUNK = 500
MAX_LEN = 6000 #Found by trial and error

if len(sys.argv) != 6:
   stop_err("Require five arguments, organism, truncate, threads, input protein FASTA file & output tabular file")

organism = sys.argv[1]
if organism not in ["euk", "gram+", "gram-"]:
   stop_err("Organism argument %s is not one of euk, gram+ or gram-" % organism)

try:
   truncate = int(sys.argv[2])
except:
   truncate = 0
if truncate < 0:
   stop_err("Truncate argument %s is not a positive integer (or zero)" % sys.argv[2])

try:
   num_threads = int(sys.argv[3])
except:
   num_threads = 0
if num_threads < 1:
   stop_err("Threads argument %s is not a positive integer" % sys.argv[3])

fasta_file = sys.argv[4]

tabular_file = sys.argv[5]

def clean_tabular(raw_handle, out_handle):
    """Clean up SignalP output to make it tabular."""
    for line in raw_handle:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\r\n").split()
        assert len(parts)==21, repr(line)
        assert parts[14].startswith(parts[0])
        #Remove redundant truncated name column (col 0)
        #and put full name at start (col 14)
        parts = parts[14:15] + parts[1:14] + parts[15:]
        out_handle.write("\t".join(parts) + "\n")

fasta_files = split_fasta(fasta_file, tabular_file, n=FASTA_CHUNK, truncate=truncate, max_len=MAX_LEN)
temp_files = [f+".out" for f in fasta_files]
assert len(fasta_files) == len(temp_files)
jobs = ["signalp -short -t %s %s > %s" % (organism, fasta, temp)
        for (fasta, temp) in zip(fasta_files, temp_files)]
assert len(fasta_files) == len(temp_files) == len(jobs)

def clean_up(file_list):
    for f in file_list:
        if os.path.isfile(f):
            os.remove(f)

if len(jobs) > 1 and num_threads > 1:
    #A small "info" message for Galaxy to show the user.
    print "Using %i threads for %i tasks" % (min(num_threads, len(jobs)), len(jobs))
results = run_jobs(jobs, num_threads)
assert len(fasta_files) == len(temp_files) == len(jobs)
for fasta, temp, cmd in zip(fasta_files, temp_files, jobs):
    error_level = results[cmd]
    output = open(temp).readline()
    if error_level or output.lower().startswith("error running"):
        clean_up(fasta_files)
        clean_up(temp_files)
        stop_err("One or more tasks failed, e.g. %i from %r gave:\n%s" % (error_level, cmd, output),
                 error_level)
del results

out_handle = open(tabular_file, "w")
fields = ["ID"]
#NN results:
for name in ["Cmax", "Ymax", "Smax"]:
    fields.extend(["NN_%s_score"%name, "NN_%s_pos"%name, "NN_%s_pred"%name])
fields.extend(["NN_Smean_score", "NN_Smean_pred", "NN_D_score", "NN_D_pred"])
#HMM results:
fields.extend(["HMM_type", "HMM_Cmax_score", "HMM_Cmax_pos", "HMM_Cmax_pred",
               "HMM_Sprob_score", "HMM_Sprob_pred"])
out_handle.write("#" + "\t".join(fields) + "\n")
for temp in temp_files:
    data_handle = open(temp)
    clean_tabular(data_handle, out_handle)
    data_handle.close()
out_handle.close()

clean_up(fasta_files)
clean_up(temp_files)
