#!/usr/bin/env python
"""Implements assorted RXLR motif methods from the literature

This script takes exactly four command line arguments:
 * Protein FASTA filename
 * Number of threads
 * Model name (Bhattacharjee2006, Win2007, Whisson2007)
 * Output tabular filename

The model names are:

Bhattacharjee2006: Simple regular expression search for RXLR
with additional requirements for positioning and signal peptide.

Win2007: Simple regular expression search for RXLR, but with
different positional requirements.

Whisson2007: As Bhattacharjee2006 but with a more complex regular
expression to look for RXLR-EER domain, and additionally calls HMMER.

See the help text in the accompanying Galaxy tool XML file for more
details including the full references.

Note:

Bhattacharjee et al. (2006) and Win et al. (2007) used SignalP v2.0,
which is no longer available. The current release is SignalP v3.0
(Mar 5, 2007). We have therefore opted to use the NN Ymax position for
the predicted cleavage site, as this is expected to be more accurate.
Also note that the HMM score values have changed from v2.0 to v3.0.
Whisson et al. (2007) used SignalP v3.0 anyway.

Whisson et al. (2007) used HMMER 2.3.2, and althought their HMM model
can still be used with hmmsearch from HMMER 3 this this does give
slightly different results. We expect the hmmsearch from HMMER 2.3.2
(the last stable release of HMMER 2) to be present on the path under
the name hmmsearch2 (allowing it to co-exist with HMMER 3).
"""
import os
import sys
import re
import subprocess
from seq_analysis_utils import stop_err, fasta_iterator

if len(sys.argv) != 5:
   stop_err("Requires four arguments: protein FASTA filename, threads, model, and output filename")

fasta_file, threads, model, tabular_file = sys.argv[1:]
hmm_output_file = tabular_file + ".hmm.tmp"
signalp_input_file = tabular_file + ".fasta.tmp"
signalp_output_file = tabular_file + ".tabular.tmp"
min_signalp_hmm = 0.9
hmmer_search = "hmmsearch2"

if model == "Bhattacharjee2006":
   signalp_trunc = 70
   re_rxlr = re.compile("R.LR")
   min_sp = 10
   max_sp = 40
   max_sp_rxlr = 100
   min_rxlr_start = 1
   #Allow signal peptide to be at most 40aa, and want RXLR to be
   #within 100aa, therefore for the prescreen the max start is 140:
   max_rxlr_start = max_sp + max_sp_rxlr
elif model == "Win2007":
   signalp_trunc = 70
   re_rxlr = re.compile("R.LR")
   min_sp = 10
   max_sp = 40
   min_rxlr_start = 30
   max_rxlr_start = 60
   #No explicit limit on separation of signal peptide clevage
   #and RXLR, but shortest signal peptide is 10, and furthest
   #away RXLR is 60, so effectively limit is 50.
   max_sp_rxlr = max_rxlr_start - min_sp + 1
elif model == "Whisson2007":
   signalp_trunc = 0 #zero for no truncation
   re_rxlr = re.compile("R.LR.{,40}[ED][ED][KR]")
   min_sp = 10
   max_sp = 40
   max_sp_rxlr = 100
   min_rxlr_start = 1
   max_rxlr_start = max_sp + max_sp_rxlr
else:
   stop_err("Did not recognise the model name %r\n"
            "Use Bhattacharjee2006, Win2007, or Whisson2007" % model)


def get_hmmer_version(exe, required=None):
    cmd = "%s -h" % exe
    try:
        child = subprocess.Popen([exe, "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        raise ValueError("Could not run %s" % exe)
    stdout, stderr = child.communicate()
    if required:
        return required in stdout
    elif "HMMER 2" in stdout:
        return 2
    elif "HMMER 3" in stdout:
        return 3
    else:
        raise ValueError("Could not determine version of %s" % exe)
    

#Run hmmsearch for Whisson et al. (2007)
if model == "Whisson2007":
    hmm_file = os.path.join(os.path.split(sys.argv[0])[0],
                       "whisson_et_al_rxlr_eer_cropped.hmm")
    if not os.path.isfile(hmm_file):
        stop_err("Missing HMM file for Whisson et al. (2007)")
    if not get_hmmer_version(hmmer_search, "HMMER 2.3.2 (Oct 2003)"):
        stop_err("Missing HMMER 2.3.2 (Oct 2003) binary, %s" % hmmer_searcher)
    #I've left the code to handle HMMER 3 in situ, in case
    #we revisit the choice to insist on HMMER 2.
    hmmer3 = (3 == get_hmmer_version(hmmer_search))
    #Using zero (or 5.6?) for bitscore threshold
    if hmmer3:
        #The HMMER3 table output is easy to parse
        #In HMMER3 can't use both -T and -E
        cmd = "%s -T 0 --tblout %s --noali %s %s > /dev/null" \
              % (hmmer_search, hmm_output_file, hmm_file, fasta_file)
    else:
        #For HMMER2 we are stuck with parsing stdout
        #Put 1e6 to effectively have no expectation threshold (otherwise
        #HMMER defaults to 10 and the calculated e-value depends on the
        #input FASTA file, and we can loose hits of interest).
        cmd = "%s -T 0 -E 1e6 %s %s > %s" \
              % (hmmer_search, hmm_file, fasta_file, hmm_output_file)
    return_code = os.system(cmd)
    if return_code:
        stop_err("Error %i from hmmsearch:\n%s" % (return_code, cmd))
    hmm_hits = set()
    valid_ids = set()
    for title, seq in fasta_iterator(fasta_file):
        name = title.split(None,1)[0]
        if name in valid_ids:
            stop_err("Duplicated identifier %r" % name)
        else:
            valid_ids.add(name)
    handle = open(hmm_output_file)
    for line in handle:
        if not line.strip():
            #We expect blank lines in the HMMER2 stdout
            continue
        elif line.startswith("#"):
            #Header
            continue
        else:
            name = line.split(None,1)[0]
            #Should be a sequence name in the HMMER3 table output.
            #Could be anything in the HMMER2 stdout.
            if name in valid_ids:
                hmm_hits.add(name)
            elif hmmer3:
                stop_err("Unexpected identifer %r in hmmsearch output" % name)
    handle.close()
    #if hmmer3:
    #    print "HMMER3 hits for %i/%i" % (len(hmm_hits), len(valid_ids))
    #else:
    #    print "HMMER2 hits for %i/%i" % (len(hmm_hits), len(valid_ids))  
    #print "%i/%i matched HMM" % (len(hmm_hits), len(valid_ids))
    os.remove(hmm_output_file)
    del valid_ids


#Prepare short list of candidates containing RXLR to pass to SignalP
assert min_rxlr_start > 0, "Min value one, since zero based counting"
count = 0
total = 0
handle = open(signalp_input_file, "w")
for title, seq in fasta_iterator(fasta_file):
    total += 1
    name = title.split(None,1)[0]
    match = re_rxlr.search(seq[min_rxlr_start-1:].upper())
    if match and min_rxlr_start - 1 + match.start() + 1 <= max_rxlr_start:
        #This is a potential RXLR, depending on the SignalP results.
        #Might as well truncate the sequence now, makes the temp file smaller
        if signalp_trunc:
            handle.write(">%s (truncated)\n%s\n" % (name, seq[:signalp_trunc]))
        else:
            #Does it matter we don't line wrap?
            handle.write(">%s\n%s\n" % (name, seq))
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
#print "SignalP done"

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


#Parse SignalP results and apply the strict RXLR criteria
total = 0
tally = dict()
handle = open(tabular_file, "w")
handle.write("#ID\t%s\n" % model)
signalp_results = parse_signalp(signalp_output_file)
for title, seq in fasta_iterator(fasta_file):
    total += 1
    rxlr = "N"
    name = title.split(None,1)[0]
    match = re_rxlr.search(seq[min_rxlr_start-1:].upper())
    if match and min_rxlr_start - 1 + match.start() + 1 <= max_rxlr_start:
        del match
        #This was the criteria for calling SignalP,
        #so it will be in the SignalP results.
        sp_id, sp_hmm_score, sp_nn_len = signalp_results.next()
        assert name == sp_id, "%s vs %s" % (name, sp_id)
        if sp_hmm_score >= min_signalp_hmm and min_sp <= sp_nn_len <= max_sp:
            match = re_rxlr.search(seq[sp_nn_len:].upper())
            if match and match.start() + 1 <= max_sp_rxlr: #1-based counting
                rxlr_start = sp_nn_len + match.start() + 1
                if min_rxlr_start <= rxlr_start <= max_rxlr_start:
                    rxlr = "Y"
    if model == "Whisson2007":
        #Combine the signalp with regular expression heuristic and the HMM
        if name in hmm_hits and rxlr == "N":
            rxlr = "hmm" #HMM only
        elif name not in hmm_hits and rxlr == "Y":
            rxlr = "re" #Heuristic only
        #Now have a four way classifier: Y, hmm, re, N
        #and count is the number of Y results (both HMM and heuristic)
    handle.write("%s\t%s\n" % (name, rxlr))
    try:
        tally[rxlr] += 1
    except KeyError:
        tally[rxlr] = 1
handle.close()
assert sum(tally.values()) == total

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
print "%s for %i sequences:" % (model, total)
print ", ".join("%s = %i" % kv for kv in sorted(tally.iteritems()))
