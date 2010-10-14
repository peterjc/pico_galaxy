#!/usr/bin/env python
"""Wrapper for TMHMM v2.0 for use in Galaxy.

This script takes exactly two command line arguments - an input protein FASTA
filename and an output tabular filename. It then calls the standalone TMHMM
v2.0 program (not the webservice) requesting the short output (one line per
protein).

First major feature is cleaning up the tabular output. The raw output from
TMHMM v2.0 looks like this (six columns tab separated):

 gi|2781234|pdb|1JLY|B	len=304 ExpAA=0.01	First60=0.00	PredHel=0	Topology=o
 gi|4959044|gb|AAD34209.1|AF069992_1	len=600	ExpAA=0.00	First60=0.00	PredHel=0	Topology=o
 gi|671626|emb|CAA85685.1|	len=473 ExpAA=0.19	First60=0.00 PredHel=0	Topology=o
 gi|3298468|dbj|BAA31520.1|	len=107	ExpAA=59.37	First60=31.17	PredHel=3	Topology=o23-45i52-74o89-106i

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
"""
import sys
import os

def stop_err(msg, error_level=1):
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

if len(sys.argv) != 3:
   stop_err("Require two arguments, input protein FASTA file & output tabular file")
fasta_file = sys.argv[1]
tabular_file = sys.argv[2]
temp_file = tabular_file + ".tmp"

def clean_tabular(raw_handle, out_handle):
    for line in raw_handle:
        if not line:
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

cmd = "tmhmm %s > %s" % (fasta_file, temp_file)
error_level = os.system(cmd)
if error_level:
    if os.path.isfile(temp_file):
       os.remove(temp_file)
    stop_err("Failed with error level %i" % error_level, error_level)

data_handle = open(temp_file)
out_handle = open(tabular_file, "w")
out_handle.write("#ID\tlen\tExpAA\tFirst60\tPredHel\tTopology\n")
clean_tabular(data_handle, out_handle)
out_handle.close()
data_handle.close()
os.remove(temp_file)
