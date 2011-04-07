#!/usr/bin/env python
"""Python tool to implement the RXLR motif method of Win et al. (2007)

This script takes exactly three command line arguments:
 * Protein FASTA filename
 * Number of threads
 * Output tabular filename

Note:

Win et al. (2007) used SignalP v2.0, which is no longer available.
The current release is SignalP v3.0 (Mar 5, 2007). We have therefore
opted to use the NN Ymax position (column NN_Ymax_pos) for the
predicted cleavage site, as this is expected to be more accurate.
Also note that the HMM score values have changed from v2.0 to v3.0.
"""
import os
import sys
import re
from seq_analysis_utils import stop_err, fasta_iterator

re_rxlr = re.compile("R.LR")

if len(sys.argv) != 4:
   stop_err("Requires three arguments: protein FASTA filename, threads, and output filename")

fasta_file, threads, tabular_file = sys.argv[1:]

#Run SignalP (using our wrapper script to get multi-core support etc)
signalp_file = tabular_file + ".tmp"
signalp_script = os.path.join(os.path.split(sys.argv[0])[0], "signalp3.py")
if not os.path.isfile(signalp_script):
   stop_err("Error - missing signalp3.py script")
cmd = "python %s euk 70 %s %s %s" % (signalp_script, threads, fasta_file, signalp_file)
return_code = os.system(cmd)
if return_code:
   stop_err("Error %i from SignalP:\n%s" % (return_code, cmd))

def parse_signalp(filename):
    """Parse SignalP output, yield tuples of ID, HMM_Sprob_score and NN predicted signal peptide length.

    Note I leave the HHM probability as a string for loss-less output.

    For signal peptide length we use NN_Ymax_pos (minus one).
    """
    handle = open(filename)
    line = handle.readline()
    assert line.startswith("#ID\t"), line
    #Tabular
    for line in handle:
        parts = line.rstrip("\t").split("\t")
        assert len(parts)==20, repr(line)
        yield parts[0], parts[18], int(parts[5])-1
    handle.close()

handle = open(tabular_file, "w")
handle.write("#ID\tHMM_Sprob_score\tSP_len\tRXLR_start\tEER_start\tRXLR?\n")
for (title, seq),(sp_id, sp_hmm_score, sp_nn_len) \
in zip(fasta_iterator(fasta_file), parse_signalp(signalp_file)):
    assert title.split(None,1)[0] == sp_id
    rxlr = "N"
    rxlr_match = re_rxlr.search(seq[sp_nn_len:].upper())
    if rxlr_match:
       rxlr_pos = sp_nn_len + rxlr_match.start() + 1 #one based counting
       eer_pos = seq.upper().find("EER", rxlr_pos+3) #length four, zero based counting
       if eer_pos == -1:
          eer_pos = ""
       else:
          eer_pos = str(eer_pos + 1) #one based counting
       if float(sp_hmm_score) > 0.9 and 10 <= sp_nn_len <= 30 \
       and 30 <= rxlr_pos <= 60:
          rxlr = "Y"
       rxlr_pos = str(rxlr_pos) #already one based
    else:
       rxlr_pos = ""
    handle.write("%s\n" % "\t".join([sp_id, sp_hmm_score, str(sp_nn_len), rxlr_pos, eer_pos, rxlr]))
handle.close()

#Cleanup
os.remove(signalp_file)
print "Done"
