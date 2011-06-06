#!/usr/bin/env python

#Copyright 2011 by Peter Cock, James Hutton Institute (formerly SCRI), UK
#
#Licenced under the GPL (GNU General Public Licence) version 3.
#
#Based on Perl script predictNLS v1.3, copyright 2001-2005 and the later
#versions up to predictnls v1.0.17 (copright 2011), by Rajesh Nair
#(nair@rostlab.org) and Burkhard Rost (rost@rostlab.org), Rost Lab,
#Columbia University http://rostlab.org/

"""Batch mode predictNLS, for finding nuclear localization signals

This is a Python script implementing the predictNLS method, originally
written in Perl, described here:

Murat Cokol, Rajesh Nair, and Burkhard Rost.
Finding nuclear localization signals.
EMBO reports 1(5), 411-415, 2000

http://dx.doi.org/10.1093/embo-reports/kvd092

The original Perl script was designed to work on a single sequence at a time,
but offers quite detailed output, including HTML (webpage).

This Python version is designed to work on a single FASTA file containing
multiple sequences, and produces a single tabular output file, with one line
per NLS found (i.e. zero or more rows per query sequence).

It takes either two or three command line arguments:

predictNLS_batch input_file output_file [nls_motif_file]

The input file should be protein sequences in FASTA format, the output file
is tab separated plain text, and the NLS motif file defaults to using the
plain text My_NLS_list file located next to the script file, or in a data
subdirectory.

Tested with the My_NLS_list file included with predictnls-1.0.17.tar.gz
"""

import os
import sys
import re

def stop_err(msg, return_code=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(return_code)

if len(sys.argv) == 4:
    fasta_filename, tabular_filename, re_filename = sys.argv[1:]
elif len(sys.argv) == 3:
    fasta_filename, tabular_filename = sys.argv[1:]
    #Use os.path.realpath(...) to handle being called via a symlink
    #Try under subdirectory data:
    re_filename = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                               "data", "My_NLS_list")
    if not os.path.isfile(re_filename):
        #Try in same directory as this script:
        re_filename = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),
                                                   "My_NLS_list")
else:
    stop_err("Expect 2 or 3 arguments: input FASTA file, output tabular file, and NLS motif file")

if not os.path.isfile(fasta_filename):
    stop_err("Could not find FASTA input file: %s" % fasta_filename)

if not os.path.isfile(re_filename):
    stop_err("Could not find NLS motif file: %s" % re_filename)

def load_re(filename):
    """Parse the 5+ column tabular NLS motif file."""
    handle = open(filename, "rU")
    for line in handle:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        assert 5 <= len(parts), parts
        regex, evidence, p_count, percent_nuc, precent_non_nuc = parts[0:5]
        try:
            regex = re.compile(regex)
            p_count = int(p_count)
        except ValueError:
            stop_err("Bad data in line: %s" % line)
        if 6 <= len(parts):
            proteins = parts[5]
            assert p_count == len(proteins.split(",")), line
        else:
            proteins = ""
            assert p_count == 0
        if 7 <= len(parts):
            domains = parts[6]
            assert int(p_count) == len(domains.split(",")), line
        else:
            domains = ""
            assert p_count == 0
        #There can be further columns (DNA binding?), but we don't use them.
        yield regex, evidence, p_count, percent_nuc, proteins, domains
    handle.close()

def fasta_iterator(filename):
    """Simple FASTA parser yielding tuples of (name, upper case sequence)."""
    handle = open(filename)
    name, seq = "", ""
    for line in handle:
        if line.startswith(">"):
            if name:
                yield name, seq
            #Take the first word only as the name:
            name = line[1:].rstrip().split(None,1)[0]
            seq = ""
        elif name:
            #Simple way would leave in any internal white space,
            #seq += line.strip().upper()
            seq += "".join(line.strip().upper().split())
        elif not line.strip():
            #Ignore blank lines before first record
            pass
        else:
            raise ValueError("Bad FASTA line %r" % line)
    handle.close()
    if name:
        yield name, seq
    raise StopIteration

motifs = list(load_re(re_filename))
print "Looking for %i NLS motifs" % len(motifs)

out_handle = open(tabular_filename, "w")
out_handle.write("#ID\tNLS start\tNLS seq\tNLS pattern\tType\tProtCount\t%NucProt\tProtList\tProtLoci\n")
count = 0
nls = 0
for idn, seq in fasta_iterator(fasta_filename):
    for regex, evidence, p_count, percent_nuc_prot, proteins, domains in motifs:
        #Perl predictnls v1.0.17 (and older) take right most hit only, Bug #40
        #This has been fixed now, so we return all the matches
        for match in regex.finditer(seq):
            #Perl predictnls v1.0.17 (and older) return NLS start
            #position with zero based counting, Bug #38 (fixed)
            #We use one based couting, hence the start+1 here:
            out_handle.write("%s\t%i\t%s\t%s\t%s\t%i\t%s\t%s\t%s\n" \
                             % (idn, match.start()+1, match.group(),
                                regex.pattern, evidence, p_count,
                                percent_nuc_prot, proteins, domains))
            nls += 1
    count += 1
out_handle.close()
print "Found %i NLS motifs in %i sequences" % (nls, count)
