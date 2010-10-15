"""A few useful functions for working with FASTA files and running jobs.

This module was originally written to hold common code used in both the TMHMM
and SignalP wrappers in Galaxy.

Given Galaxy currently supports Python 2.4+ this cannot use the Python module
multiprocessing so the function run_jobs instead is a simple pool approach
using just the subprocess library.
"""
import sys
import os
import subprocess
from time import sleep

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

def fasta_iterator(filename):
    """Simple FASTA parser yielding tuples of (title, sequence) strings."""
    handle = open(filename)
    title, seq = "", ""
    for line in handle:
        if line.startswith(">"):
            if title:
                yield title, seq
            title = line[1:].rstrip()
            seq = ""
        elif title:
            seq += line.strip()
        elif not line.strip() or line.startswith("#"):
            #Ignore blank lines, and any comment lines
            #between records (starting with hash).
            pass
        else:
            raise ValueError("Bad FASTA line %r" % line)
    handle.close()
    if title:
        yield title, seq
    raise StopIteration

def split_fasta(filename, n=500, truncate=None):
    """Split FASTA file into sub-files each of at most n sequences.

    Returns a list of the filenames used (based on the input filename).
    Each sequence can also be truncated (since we only need the start for
    SignalP).
    """
    iterator = fasta_iterator(filename)
    files = []
    while True:
        records = []
        for i in range(n):
            try:
                records.append(iterator.next())
            except StopIteration:
                break
        if not records:
            break
        #TODO - Use the tempfile library?
        #That avoids potential problem where we may not have write permission
        #in the folder (can that happen in Galaxy though?)
        new_filename = "%s.%i.tmp" % (filename, len(files))
        handle = open(new_filename, "w")
        if truncate:
            for title, seq in records:
                handle.write(">%s\n%s\n" % (title, seq[:truncate]))
        else:
            for title, seq in records:
                handle.write(">%s\n" % title)
                for i in range(0, len(seq), 60):
                    handle.write(seq[i:i+60] + "\n")
        handle.close()
        files.append(new_filename)
        #print "%i records in %s" % (len(records), new_filename)
    return files

def run_jobs(jobs, threads):
    """Takes list of cmd strings, returns dict with error levels."""
    running = []
    results = {}
    while jobs or running:
        #print "%i jobs pending, %i running, %i completed" \
        #      % (len(jobs), len(running), len(results))
        #See if any have finished
        for (cmd, process) in running:
            return_code = process.wait()
            if return_code is not None:
                results[cmd] = return_code
        running = [(cmd, process) for (cmd, process) in running \
                   if cmd not in results]
        #See if we can start any new threads
        while jobs and len(running) < threads:
            cmd = jobs.pop(0)
            process = subprocess.Popen(cmd, shell=True)
            running.append((cmd, process))
        #Loop...
        sleep(1)
    #print "%i jobs completed" % len(results)
    return results
