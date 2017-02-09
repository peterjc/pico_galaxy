#!/usr/bin/env python
"""Wrapper for psortb for use in Galaxy.

This script takes exactly six command line arguments - which includes the
number of threads, and the input protein FASTA filename and output
tabular filename. It then splits up the FASTA input and calls multiple
copies of the standalone psortb v3 program, then collates the output.
e.g. Rather than this,

psort $type -c $cutoff -d $divergent -o long $sequence > $outfile

Call this:

psort $threads $type $cutoff $divergent $sequence $outfile

If ommitting -c or -d options, set $cutoff and $divergent to zero or blank.

Note that this is somewhat redundant with job-splitting available in Galaxy
itself (see the SignalP XML file for settings), but both can be applied.

Additionally it ensures the header line (with the column names) starts
with a # character as used elsewhere in Galaxy.
"""

import os
import sys
import tempfile

from seq_analysis_utils import run_jobs, split_fasta, thread_count

FASTA_CHUNK = 500

if "-v" in sys.argv or "--version" in sys.argv:
    """Return underlying PSORTb's version"""
    sys.exit(os.system("psort --version"))

if len(sys.argv) != 8:
    sys.exit("Require 7 arguments, number of threads (int), type (e.g. archaea), "
             "output (e.g. terse/normal/long), cutoff, divergent, input protein "
             "FASTA file & output tabular file")

num_threads = thread_count(sys.argv[1], default=4)
org_type = sys.argv[2]
out_type = sys.argv[3]
cutoff = sys.argv[4]
if cutoff.strip() and float(cutoff.strip()) != 0.0:
    cutoff = "-c %s" % cutoff
else:
    cutoff = ""
divergent = sys.argv[5]
if divergent.strip() and float(divergent.strip()) != 0.0:
    divergent = "-d %s" % divergent
else:
    divergent = ""
fasta_file = sys.argv[6]
tabular_file = sys.argv[7]

if out_type == "terse":
    header = ['SeqID', 'Localization', 'Score']
elif out_type == "normal":
    sys.exit("Normal output not implemented yet, sorry.")
elif out_type == "long":
    if org_type == "-n":
        # Gram negative bacteria
        header = ['SeqID', 'CMSVM-_Localization', 'CMSVM-_Details', 'CytoSVM-_Localization', 'CytoSVM-_Details',
                  'ECSVM-_Localization', 'ECSVM-_Details', 'ModHMM-_Localization', 'ModHMM-_Details',
                  'Motif-_Localization', 'Motif-_Details', 'OMPMotif-_Localization', 'OMPMotif-_Details',
                  'OMSVM-_Localization', 'OMSVM-_Details', 'PPSVM-_Localization', 'PPSVM-_Details',
                  'Profile-_Localization', 'Profile-_Details',
                  'SCL-BLAST-_Localization', 'SCL-BLAST-_Details', 'SCL-BLASTe-_Localization', 'SCL-BLASTe-_Details',
                  'Signal-_Localization', 'Signal-_Details',
                  'Cytoplasmic_Score', 'CytoplasmicMembrane_Score', 'Periplasmic_Score', 'OuterMembrane_Score',
                  'Extracellular_Score', 'Final_Localization', 'Final_Localization_Details', 'Final_Score',
                  'Secondary_Localization', 'PSortb_Version']
    elif org_type == "-p":
        # Gram positive bacteria
        header = ['SeqID', 'CMSVM+_Localization', 'CMSVM+_Details', 'CWSVM+_Localization', 'CWSVM+_Details',
                  'CytoSVM+_Localization', 'CytoSVM+_Details', 'ECSVM+_Localization', 'ECSVM+_Details',
                  'ModHMM+_Localization', 'ModHMM+_Details', 'Motif+_Localization', 'Motif+_Details',
                  'Profile+_Localization', 'Profile+_Details',
                  'SCL-BLAST+_Localization', 'SCL-BLAST+_Details', 'SCL-BLASTe+_Localization', 'SCL-BLASTe+_Details',
                  'Signal+_Localization', 'Signal+_Details',
                  'Cytoplasmic_Score', 'CytoplasmicMembrane_Score', 'Cellwall_Score',
                  'Extracellular_Score', 'Final_Localization', 'Final_Localization_Details', 'Final_Score',
                  'Secondary_Localization', 'PSortb_Version']
    elif org_type == "-a":
        # Archaea
        header = ['SeqID', 'CMSVM_a_Localization', 'CMSVM_a_Details', 'CWSVM_a_Localization', 'CWSVM_a_Details',
                  'CytoSVM_a_Localization', 'CytoSVM_a_Details', 'ECSVM_a_Localization', 'ECSVM_a_Details',
                  'ModHMM_a_Localization', 'ModHMM_a_Details', 'Motif_a_Localization', 'Motif_a_Details',
                  'Profile_a_Localization', 'Profile_a_Details',
                  'SCL-BLAST_a_Localization', 'SCL-BLAST_a_Details', 'SCL-BLASTe_a_Localization', 'SCL-BLASTe_a_Details',
                  'Signal_a_Localization', 'Signal_a_Details',
                  'Cytoplasmic_Score', 'CytoplasmicMembrane_Score', 'Cellwall_Score',
                  'Extracellular_Score', 'Final_Localization', 'Final_Localization_Details', 'Final_Score',
                  'Secondary_Localization', 'PSortb_Version']
    else:
        sys.exit("Expected -n, -p or -a for the organism type, not %r" % org_type)
else:
    sys.exit("Expected terse, normal or long for the output type, not %r" % out_type)

tmp_dir = tempfile.mkdtemp()


def clean_tabular(raw_handle, out_handle):
    """Clean up tabular TMHMM output, returns output line count."""
    global header
    count = 0
    for line in raw_handle:
        if not line.strip() or line.startswith("#"):
            # Ignore any blank lines or comment lines
            continue
        parts = [x.strip() for x in line.rstrip("\r\n").split("\t")]
        if parts == header:
            # Ignore the header line
            continue
        if not parts[-1] and len(parts) == len(header) + 1:
            # Ignore dummy blank extra column, e.g.
            # "...2.0\t\tPSORTb version 3.0\t\n"
            parts = parts[:-1]
        assert len(parts) == len(header), \
            "%i fields, not %i, in line:\n%r" % (len(line), len(header), line)
        out_handle.write(line)
        count += 1
    return count


# Note that if the input FASTA file contains no sequences,
# split_fasta returns an empty list (i.e. zero temp files).
fasta_files = split_fasta(fasta_file, os.path.join(tmp_dir, "tmhmm"), FASTA_CHUNK)
temp_files = [f + ".out" for f in fasta_files]
jobs = ["psort %s %s %s -o %s %s > %s" % (org_type, cutoff, divergent, out_type, fasta, temp)
        for fasta, temp in zip(fasta_files, temp_files)]


def clean_up(file_list):
    for f in file_list:
        if os.path.isfile(f):
            os.remove(f)
    try:
        os.rmdir(tmp_dir)
    except Exception:
        pass


if len(jobs) > 1 and num_threads > 1:
    # A small "info" message for Galaxy to show the user.
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
        sys.exit("One or more tasks failed, e.g. %i from %r gave:\n%s" % (error_level, cmd, output),
                 error_level)
del results
del jobs

out_handle = open(tabular_file, "w")
out_handle.write("#%s\n" % "\t".join(header))
count = 0
for temp in temp_files:
    data_handle = open(temp)
    count += clean_tabular(data_handle, out_handle)
    data_handle.close()
    if not count:
        clean_up(fasta_files + temp_files)
        sys.exit("No output from psortb")
out_handle.close()
print "%i records" % count

clean_up(fasta_files + temp_files)
