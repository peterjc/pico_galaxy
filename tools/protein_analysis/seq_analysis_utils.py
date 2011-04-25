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

__version__ = "0.0.1"

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

def fasta_iterator(filename, max_len=None, truncate=None):
    """Simple FASTA parser yielding tuples of (title, sequence) strings."""
    handle = open(filename)
    title, seq = "", ""
    for line in handle:
        if line.startswith(">"):
            if title:
                if truncate:
                    seq = seq[:truncate]
                if max_len and len(seq) > max_len:
                    raise ValueError("Sequence %s is length %i, max length %i" \
                                     % (title.split()[0], len(seq), max_len))
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
        if truncate:
            seq = seq[:truncate]
        if max_len and len(seq) > max_len:
            raise ValueError("Sequence %s is length %i, max length %i" \
                             % (title.split()[0], len(seq), max_len))
        yield title, seq
    raise StopIteration

def split_fasta(input_filename, output_filename_base, n=500, truncate=None, keep_descr=False, max_len=None):
    """Split FASTA file into sub-files each of at most n sequences.

    Returns a list of the filenames used (based on the input filename).
    Each sequence can also be truncated (since we only need the start for
    SignalP), and have its description discarded (since we don't usually
    care about it and some tools don't like very long title lines).

    If a max_len is given and any sequence exceeds it no temp files are
    created and an exception is raised.
    """
    iterator = fasta_iterator(input_filename, max_len, truncate)
    files = []
    try:
        while True:
            records = []
            for i in range(n):
                try:
                    records.append(iterator.next())
                except StopIteration:
                    break
            if not records:
                break
            new_filename = "%s.%i.tmp" % (output_filename_base, len(files))
            handle = open(new_filename, "w")
            if keep_descr:
                for title, seq in records:
                    handle.write(">%s\n" % title)
                    for i in range(0, len(seq), 60):
                        handle.write(seq[i:i+60] + "\n")
            else:
                for title, seq in records:
                    handle.write(">%s\n" % title.split()[0])
                    for i in range(0, len(seq), 60):
                        handle.write(seq[i:i+60] + "\n")
            handle.close()
            files.append(new_filename)
            #print "%i records in %s" % (len(records), new_filename)
    except ValueError, err:
        #Max length failure from parser - clean up
        try:
            handle.close()
        except:
            pass
        for f in files:
            if os.path.isfile(f):
                os.remove(f)
        raise err
    return files

def run_jobs(jobs, threads, pause=10, verbose=False):
    """Takes list of cmd strings, returns dict with error levels."""
    pending = jobs[:]
    running = []
    results = {}
    while pending or running:
        #See if any have finished
        for (cmd, process) in running:
            return_code = process.poll() #non-blocking
            if return_code is not None:
                results[cmd] = return_code
        running = [(cmd, process) for (cmd, process) in running \
                   if cmd not in results]
        if verbose:
            print "%i jobs pending, %i running, %i completed" \
                  % (len(pending), len(running), len(results))
        #See if we can start any new threads
        while pending and len(running) < threads:
            cmd = pending.pop(0)
            if verbose:
                print cmd
            process = subprocess.Popen(cmd, shell=True)
            running.append((cmd, process))
        #Loop...
        sleep(pause)
    if verbose:
        print "%i jobs completed" % len(results)
    assert set(jobs) == set(results)
    return results
