#!/usr/bin/env python
"""Wrapper for WoLF PSORT v0.2 for use in Galaxy.

This script takes exactly four command line arguments:
 * the organism type (animal, plant or fungi)
 * number of threads to use (integer)
 * an input protein FASTA filename
 * output tabular filename.

It then calls the standalone WoLF PSORT v0.2 program runWolfPsortSummary
(not the webservice), and coverts the output from something like this:

# k used for kNN is: 27
gi|301087619|ref|XP_002894699.1| extr 12, mito 4, E.R. 3, golg 3, mito_nucl 3
gi|301087623|ref|XP_002894700.1| extr 21, mito 2, cyto 2, cyto_mito 2

In order to make it easier to use in Galaxy, this wrapper script reformats
this to use tab separators, with one line per compartment prediction:

#ID	Compartment	Score	Rank
gi|301087619|ref|XP_002894699.1|	extr	12	1
gi|301087619|ref|XP_002894699.1|	mito	4	2
gi|301087619|ref|XP_002894699.1|	E.R.	3	3
gi|301087619|ref|XP_002894699.1|	golg	3	4
gi|301087619|ref|XP_002894699.1|	mito_nucl	3	5
gi|301087623|ref|XP_002894700.1|	extr	21	1
gi|301087623|ref|XP_002894700.1|	mito	2	2
gi|301087623|ref|XP_002894700.1|	cyto	2	3
gi|301087623|ref|XP_002894700.1|	cyto_mito	2	4

Additionally in order to take full advantage of multiple cores, by subdividing
the input FASTA file multiple copies of WoLF PSORT are run in parallel. I would
normally use Python's multiprocessing library in this situation but it requires
at least Python 2.6 and at the time of writing Galaxy still supports Python 2.4.
"""
import sys
import os
from seq_analysis_utils import sys_exit, split_fasta, run_jobs, thread_count

FASTA_CHUNK = 500
exe = "runWolfPsortSummary"

"""
Note: I had trouble getting runWolfPsortSummary on the path (via a link), other
than by including all of /opt/WoLFPSORT_package_v0.2/bin , so used a wrapper
python script called runWolfPsortSummary as follows:

#!/usr/bin/env python
#Wrapper script to call WoLF PSORT from its own directory.
import os
import sys
import subprocess
saved_dir = os.path.abspath(os.curdir)
os.chdir("/opt/WoLFPSORT_package_v0.2/bin")
args = ["./runWolfPsortSummary"] + sys.argv[1:]
return_code = subprocess.call(args)
os.chdir(saved_dir)
sys.exit(return_code)

For more details on this workaround, see:
https://lists.galaxyproject.org/pipermail/galaxy-dev/2015-December/023386.html
"""

if len(sys.argv) != 5:
    sys_exit("Require four arguments, organism, threads, input protein FASTA file & output tabular file")

organism = sys.argv[1]
if organism not in ["animal", "plant", "fungi"]:
    sys_exit("Organism argument %s is not one of animal, plant, fungi" % organism)

num_threads = thread_count(sys.argv[2], default=4)
fasta_file = sys.argv[3]
tabular_file = sys.argv[4]

def clean_tabular(raw_handle, out_handle):
    """Clean up WoLF PSORT output to make it tabular."""
    for line in raw_handle:
        if not line or line.startswith("#"):
            continue
        name, data = line.rstrip("\r\n").split(None,1)
        for rank, comp_data in enumerate(data.split(",")):
            comp, score = comp_data.split()
            out_handle.write("%s\t%s\t%s\t%i\n" \
                             % (name, comp, score, rank+1))

fasta_files = split_fasta(fasta_file, tabular_file, n=FASTA_CHUNK)
temp_files = [f+".out" for f in fasta_files]
assert len(fasta_files) == len(temp_files)
jobs = ["%s %s < %s > %s" % (exe, organism, fasta, temp)
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
    try:
        output = open(temp).readline()
    except IOError:
        output = ""
    if error_level or output.lower().startswith("error running"):
        clean_up(fasta_files)
        clean_up(temp_files)
        sys_exit("One or more tasks failed, e.g. %i from %r gave:\n%s" % (error_level, cmd, output),
                 error_level)
del results

out_handle = open(tabular_file, "w")
out_handle.write("#ID\tCompartment\tScore\tRank\n")
for temp in temp_files:
    data_handle = open(temp)
    clean_tabular(data_handle, out_handle)
    data_handle.close()
out_handle.close()

clean_up(fasta_files)
clean_up(temp_files)
