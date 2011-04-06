#!/usr/bin/env python
"""Python tool to implement the RXLR motif method of Win et al. (2007)

This script takes exactly three command line arguments:
 * Protein FASTA file
 * SignalP v3.0 predictions for this FASTA file (tabular)
 * Output tabular filename.

Note:

Win et al. (2007) used SignalP v2.0, which is no longer available.
The current release is SignalP v3.0 (Mar 5, 2007). We have therefore
opted to use the NN Ymax position (column NN_Ymax_pos) for the
predicted cleavage site, as this is expected to be more accurate.
Also note that the HMM score values have changed from v2.0 to v3.0.

Note:

Raw SignalP v3.0 output looks like this (21 columns space separated):

# SignalP-NN euk predictions                                   	                # SignalP-HMM euk predictions
# name                Cmax  pos ?  Ymax  pos ?  Smax  pos ?  Smean ?  D     ? 	# name      !  Cmax  pos ?  Sprob ?
gi|2781234|pdb|1JLY|  0.061  17 N  0.043  17 N  0.199   1 N  0.067 N  0.055 N	gi|2781234|pdb|1JLY|B  Q  0.000  17 N  0.000 N  
gi|4959044|gb|AAD342  0.099 191 N  0.012  38 N  0.023  12 N  0.014 N  0.013 N	gi|4959044|gb|AAD34209.1|AF069992_1  Q  0.000   0 N  0.000 N  
gi|671626|emb|CAA856  0.139 381 N  0.020   8 N  0.121   4 N  0.067 N  0.044 N	gi|671626|emb|CAA85685.1|  Q  0.000   0 N  0.000 N  
gi|3298468|dbj|BAA31  0.208  24 N  0.184  38 N  0.980  32 Y  0.613 Y  0.398 N	gi|3298468|dbj|BAA31520.1|  Q  0.066  24 N  0.139 N

In order to make it easier to use in Galaxy, our wrapper script reformats
this to use tab separators. Also it removes the redundant truncated name
column, and assigns unique column names in the header:

#ID	NN_Cmax_score	NN_Cmax_pos	NN_Cmax_pred	NN_Ymax_score	NN_Ymax_pos	NN_Ymax_pred	NN_Smax_score	NN_Smax_pos	NN_Smax_pred	NN_Smean_score	NN_Smean_pred	NN_D_score	NN_D_pred	HMM_bang	HMM_Cmax_score	HMM_Cmax_pos	HMM_Cmax_pred	HMM_Sprob_score	HMM_Sprob_pred
gi|2781234|pdb|1JLY|B	0.061	17	N	0.043	17	N	0.199	1	N	0.067	N	0.055	N	Q	0.000	17	N	0.000	N
gi|4959044|gb|AAD34209.1|AF069992_1	0.099	191	N	0.012	38	N	0.023	12	N	0.014	N	0.013	N	Q	0.000	0	N	0.000	N
gi|671626|emb|CAA85685.1|	0.139	381	N	0.020	8	N	0.121	4	N	0.067	N	0.044	N	Q	0.000	0	N	0.000	N
gi|3298468|dbj|BAA31520.1|	0.208	24	N	0.184	38	N	0.980	32	Y	0.613	Y	0.398	N	Q	0.066	24	N	0.139	N

This tool expects either style.
"""
import sys
import re
from seq_analysis_utils import stop_err, fasta_iterator

re_rxlr = re.compile("R.LR")

if len(sys.argv) != 4:
   stop_err("Requires three arguments: protein FASTA, SignalP predictions, and output filename")

fasta_file, signalp_file, tabular_file = sys.argv[1:]

def parse_signalp(filename):
    """Parse SignalP output, yield tuples of ID, HMM_Sprob_score and NN predicted signal peptide length.

    Note I leave the HHM probability as a string for loss-less output.

    For signal peptide length we use NN_Ymax_pos (minus one).
    """
    handle = open(filename)
    line = handle.readline()
    if line.startswith("#ID\t"):
        #Tabular
        for line in handle:
           parts = line.rstrip("\t").split("\t")
           assert len(parts)==20, repr(line)
           yield parts[0], parts[18], int(parts[5])-1
    else:
        #Native
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split()
            assert len(parts)==21, repr(line)
            assert parts[14].startswith(parts[0])
            #Use the full name (col 14)
            yield parts[14], parts[19], int(parts[5])-1
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
print "Done"
