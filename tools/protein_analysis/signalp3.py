#!/usr/bin/env python
"""Wrapper for SignalP v3.0 for use in Galaxy.

This script takes exactly five command line arguments:
 * the organism type (euk, gram+ or gram-)
 * length to truncate sequences to (integer)
 * number of threads to use (integer, defaults to one)
 * an input protein FASTA filename
 * output tabular filename.

There are two further optional arguments
 * cut type (NN_Cmax, NN_Ymax, NN_Smax or HMM_Cmax)
 * output GFF3 filename

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

Note that this is somewhat redundant with job-splitting available in Galaxy
itself (see the SignalP XML file for settings).

Finally, you can opt to have a GFF3 file produced which will describe the
predicted signal peptide and mature peptide for each protein (using one of
the predictors which gives a cleavage site). *WORK IN PROGRESS*
"""  # noqa: E501

from __future__ import print_function

import os
import sys
import tempfile

from seq_analysis_utils import fasta_iterator, split_fasta
from seq_analysis_utils import run_jobs, thread_count

FASTA_CHUNK = 500
MAX_LEN = 6000  # Found by trial and error

if "-v" in sys.argv or "--version" in sys.argv:
    print("SignalP Galaxy wrapper version 0.0.19")
    sys.exit(os.system("signalp -version"))

if len(sys.argv) not in [6, 8]:
    sys.exit("Require five (or 7) arguments, organism, truncate, threads, "
             "input protein FASTA file & output tabular file (plus "
             "optionally cut method and GFF3 output file). "
             "Got %i arguments." % (len(sys.argv) - 1))

organism = sys.argv[1]
if organism not in ["euk", "gram+", "gram-"]:
    sys.exit("Organism argument %s is not one of euk, gram+ or gram-" % organism)

try:
    truncate = int(sys.argv[2])
except ValueError:
    truncate = 0
if truncate < 0:
    sys.exit("Truncate argument %s is not a positive integer (or zero)" % sys.argv[2])

num_threads = thread_count(sys.argv[3], default=4)
fasta_file = sys.argv[4]
tabular_file = sys.argv[5]

if len(sys.argv) == 8:
    cut_method = sys.argv[6]
    if cut_method not in ["NN_Cmax", "NN_Ymax", "NN_Smax", "HMM_Cmax"]:
        sys.exit("Invalid cut method %r" % cut_method)
    gff3_file = sys.argv[7]
else:
    cut_method = None
    gff3_file = None


tmp_dir = tempfile.mkdtemp()


def clean_tabular(raw_handle, out_handle, gff_handle=None):
    """Clean up SignalP output to make it tabular."""
    for line in raw_handle:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\r\n").split()
        out_handle.write("\t".join(parts) + "\n")


def make_gff(fasta_file, tabular_file, gff_file, cut_method):
    """Make a GFF file."""
    cut_col, score_col = {"NN_Cmax": (2, 1),
                          "NN_Ymax": (5, 4),
                          "NN_Smax": (8, 7),
                          "HMM_Cmax": (16, 15),
                          }[cut_method]

    source = "SignalP"
    strand = "."  # not stranded
    phase = "."  # not phased
    tags = "Note=%s" % cut_method

    tab_handle = open(tabular_file)
    line = tab_handle.readline()
    assert line.startswith("#ID\t"), line

    gff_handle = open(gff_file, "w")
    gff_handle.write("##gff-version 3\n")

    for (title, seq), line in zip(fasta_iterator(fasta_file), tab_handle):
        parts = line.rstrip("\n").split("\t")
        seqid = parts[0]
        assert title.startswith(seqid), "%s vs %s" % (seqid, title)
        if not seq:
            # Is it possible to have a zero length reference in GFF3?
            continue
        cut = int(parts[cut_col])
        if cut == 0:
            assert cut_method == "HMM_Cmax", cut_method
            # TODO - Why does it do this?
            cut = 1
        assert 1 <= cut <= len(seq), "%i for %s len %i" % (cut, seqid, len(seq))
        score = parts[score_col]
        gff_handle.write("##sequence-region %s %i %i\n"
                         % (seqid, 1, len(seq)))
        # If the cut is at the very begining, there is no signal peptide!
        if cut > 1:
            # signal_peptide = SO:0000418
            gff_handle.write("%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n"
                             % (seqid, source,
                                "signal_peptide", 1, cut - 1,
                                score, strand, phase, tags))
        # mature_protein_region = SO:0000419
        gff_handle.write("%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n"
                         % (seqid, source,
                            "mature_protein_region", cut, len(seq),
                            score, strand, phase, tags))
    tab_handle.close()
    gff_handle.close()


if num_threads == 1:
    # Still want to call split_fasta to apply truncation, but
    # no reason to make multiple files - and more chance of
    # hitting file system glitches if we do. So,
    FASTA_CHUNK = sys.maxsize

fasta_files = split_fasta(fasta_file, os.path.join(tmp_dir, "signalp"),
                          n=FASTA_CHUNK, truncate=truncate, max_len=MAX_LEN)
temp_files = [f + ".out" for f in fasta_files]
assert len(fasta_files) == len(temp_files)
jobs = ["signalp -f short -t %s %s > %s" % (organism, fasta, temp)
        for (fasta, temp) in zip(fasta_files, temp_files)]
assert len(fasta_files) == len(temp_files) == len(jobs)


def clean_up(file_list):
    """Remove temp files, and if possible the temp directory."""
    for f in file_list:
        if os.path.isfile(f):
            os.remove(f)
    try:
        os.rmdir(tmp_dir)
    except Exception:
        pass


if len(jobs) > 1 and num_threads > 1:
    # A small "info" message for Galaxy to show the user.
    print("Using %i threads for %i tasks" % (min(num_threads, len(jobs)), len(jobs)))
results = run_jobs(jobs, num_threads)
assert len(fasta_files) == len(temp_files) == len(jobs)
for fasta, temp, cmd in zip(fasta_files, temp_files, jobs):
    error_level = results[cmd]
    try:
        output = open(temp).readline()
    except IOError:
        output = "(no output)"
    if error_level or output.lower().startswith("error running"):
        clean_up(fasta_files + temp_files)
        if output:
            sys.stderr.write("One or more tasks failed, e.g. %i from %r gave:\n%s" % (error_level, cmd, output))
        else:
            sys.stderr.write("One or more tasks failed, e.g. %i from %r with no output\n" % (error_level, cmd))
        sys.exit(error_level)
del results

out_handle = open(tabular_file, "w")
fields = ["name", "Cmax", "pos", "Ymax", "pos", "Smax", "pos", "Smean", "D", "?", "Dmaxcut", "Networks-used"]
out_handle.write("#" + "\t".join(fields) + "\n")
for temp in temp_files:
    data_handle = open(temp)
    clean_tabular(data_handle, out_handle)
    data_handle.close()
out_handle.close()

# GFF3:
if cut_method:
    make_gff(fasta_file, tabular_file, gff3_file, cut_method)

clean_up(fasta_files + temp_files)
