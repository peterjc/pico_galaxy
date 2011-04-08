#!/usr/bin/env python
"""Implements the RXLR-EER motif regular expression method of Whisson et al. (2007)

This script takes exactly three command line arguments:
 * Protein FASTA filename
 * Number of threads
 * Output tabular filename

Note: This tool requires SignalP v3.0, and the provided wrapper script which takes
care of converting the output to tabular format and running it on multiple cores.
"""
import os
import sys
import re
from seq_analysis_utils import stop_err, fasta_iterator

re_rxlr_eer = re.compile("R.LR.{,40}[ED][ED][KR]")

if len(sys.argv) != 4:
   stop_err("Requires three arguments: protein FASTA filename, threads, and output filename")

fasta_file, threads, tabular_file = sys.argv[1:]
signalp_input_file = tabular_file + ".fasta.tmp"
signalp_output_file = tabular_file + ".tabular.tmp"
signalp_trunc = 70


#Prepare short list of candidates containing RXLR to pass to SignalP
count = 0
total = 0
handle = open(signalp_input_file, "w")
for title, seq in fasta_iterator(fasta_file):
    total += 1
    name = title.split(None,1)[0]
    match = re_rxlr_eer.search(seq.upper())
    if match and match.start() + 1 <= 140: #140 using 1-based counting
       #Allow up to 40aa for the signal peptide, and the RXLR must start
       #within 100aa after that - so at most 140aa from the protein start.
       #This is a potential RXLR, depending on the SignalP results.
       #Might as well truncate the sequence now, makes the temp file smaller
       handle.write(">%s (truncated)\n%s\n" % (name, seq[:signalp_trunc]))
       count += 1
handle.close()
#print "Running SignalP on %i/%i potentials." % (count, total)


#Run SignalP (using our wrapper script to get multi-core support etc)
signalp_script = os.path.join(os.path.split(sys.argv[0])[0], "signalp3.py")
if not os.path.isfile(signalp_script):
    stop_err("Error - missing signalp3.py script")
cmd = "python %s euk %i %s %s %s" % (signalp_script, signalp_trunc, threads, signalp_input_file, signalp_output_file)
return_code = os.system(cmd)
if return_code:
    stop_err("Error %i from SignalP:\n%s" % (return_code, cmd))


def parse_signalp(filename):
    """Parse SignalP output, yield tuples of ID, HMM_Sprob_score and NN predicted signal peptide length.

    For signal peptide length we use NN_Ymax_pos (minus one).
    """
    handle = open(filename)
    line = handle.readline()
    assert line.startswith("#ID\t"), line
    for line in handle:
        parts = line.rstrip("\t").split("\t")
        assert len(parts)==20, repr(line)
        yield parts[0], float(parts[18]), int(parts[5])-1
    handle.close()


#Parse SignalP results and apply the strict RXLR-EER criteria
count = 0
total = 0
handle = open(tabular_file, "w")
handle.write("#ID\tRXLR-EER\n")
signalp_results = parse_signalp(signalp_output_file)
for title, seq in fasta_iterator(fasta_file):
    total += 1
    rxlr_eer = "N"
    name = title.split(None,1)[0]
    match = re_rxlr_eer.search(seq.upper())
    if match and match.start() + 1 <= 140: #140 using 1-based counting
        del match
        #This was the criteria for calling SignalP,
        #so it will be in the SignalP results.
        sp_id, sp_hmm_score, sp_nn_len = signalp_results.next()
        assert name == sp_id, "%s vs %s" % (name, sp_id)
        if sp_hmm_score > 0.9 and 10 <= sp_nn_len <= 40:
            match = re_rxlr_eer.search(seq[sp_nn_len:].upper())
            if match and match.start() + 1 <= 100: #1-based counting
                #Within 100aa of the signal peptide cleavage point
                rxlr_eer = "Y"
                count += 1
    handle.write("%s\t%s\n" % (name, rxlr_eer))
handle.close()

#Check the iterator is finished
try:
    signalp_results.next()
    assert False, "Unexpected data in SignalP output"
except StopIteration:
    pass

#Cleanup
os.remove(signalp_input_file)
os.remove(signalp_output_file)

#Short summary to stdout for Galaxy's info display
print "%i out of %i have RXLR-EER motif" % (count, total)
