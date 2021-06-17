"""A few useful functions for working with FASTA files and running jobs.

This module was originally written to hold common code used in both the TMHMM
and SignalP wrappers in Galaxy.

Given Galaxy currently supports Python 2.4+ this cannot use the Python module
multiprocessing so the function run_jobs instead is a simple pool approach
using just the subprocess library.
"""

from __future__ import print_function

import os
import subprocess
import sys

from time import sleep

if sys.version_info[0] < 3:
    range = xrange  # noqa: F821

__version__ = "0.0.5"

try:
    from multiprocessing import cpu_count
except ImportError:
    # Must be under Python 2.5, this is copied from multiprocessing:
    def cpu_count():
        """Return the number of CPUs in the system."""
        if sys.platform == "win32":
            try:
                num = int(os.environ["NUMBER_OF_PROCESSORS"])
            except (ValueError, KeyError):
                num = 0
        elif "bsd" in sys.platform or sys.platform == "darwin":
            comm = "/sbin/sysctl -n hw.ncpu"
            if sys.platform == "darwin":
                comm = "/usr" + comm
                try:
                    with os.popen(comm) as p:
                        num = int(p.read())
                except ValueError:
                    num = 0
        else:
            try:
                num = os.sysconf("SC_NPROCESSORS_ONLN")
            except (ValueError, OSError, AttributeError):
                num = 0

        if num >= 1:
            return num
        else:
            raise NotImplementedError("cannot determine number of cpus")


def thread_count(command_line_arg, default=1):
    """Determine number of threads to use from the command line args."""
    try:
        num = int(command_line_arg)
    except ValueError:
        num = default
    if num < 1:
        sys.exit("Threads argument %r is not a positive integer" % command_line_arg)
    # Cap this with the pysical limit of the machine,
    try:
        num = min(num, cpu_count())
    except NotImplementedError:
        pass
    # For debugging,
    # hostname = os.environ.get("HOSTNAME", "this machine")
    # sys.stderr.write("Using %i cores on %s\n" % (num, hostname))
    return num


def fasta_iterator(filename, max_len=None, truncate=None):
    """Parse FASTA file yielding tuples of (name, sequence)."""
    handle = open(filename)
    title, seq = "", ""
    for line in handle:
        if line.startswith(">"):
            if title:
                if truncate:
                    seq = seq[:truncate]
                if max_len and len(seq) > max_len:
                    raise ValueError(
                        "Sequence %s is length %i, max length %i"
                        % (title.split()[0], len(seq), max_len)
                    )
                yield title, seq
            title = line[1:].rstrip()
            seq = ""
        elif title:
            seq += line.strip()
        elif not line.strip() or line.startswith("#"):
            # Ignore blank lines, and any comment lines
            # between records (starting with hash).
            pass
        else:
            handle.close()
            raise ValueError("Bad FASTA line %r" % line)
    handle.close()
    if title:
        if truncate:
            seq = seq[:truncate]
        if max_len and len(seq) > max_len:
            raise ValueError(
                "Sequence %s is length %i, max length %i"
                % (title.split()[0], len(seq), max_len)
            )
        yield title, seq


def split_fasta(
    input_filename,
    output_filename_base,
    n=500,
    truncate=None,
    keep_descr=False,
    max_len=None,
):
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
                    records.append(next(iterator))
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
                        handle.write(seq[i : i + 60] + "\n")
            else:
                for title, seq in records:
                    handle.write(">%s\n" % title.split()[0])
                    for i in range(0, len(seq), 60):
                        handle.write(seq[i : i + 60] + "\n")
            handle.close()
            files.append(new_filename)
            # print "%i records in %s" % (len(records), new_filename)
    except ValueError as err:
        # Max length failure from parser - clean up
        try:
            handle.close()
        except Exception:
            pass
        for f in files:
            if os.path.isfile(f):
                os.remove(f)
        raise err
    for f in files:
        assert os.path.isfile(f), "Missing split file %r (!??)" % f
    return files


def run_jobs(jobs, threads, pause=10, verbose=False, fast_fail=True):
    """Take list of cmd strings, return dict with error levels."""
    pending = jobs[:]
    running = []
    results = {}
    skipped = []
    if threads == 1:
        # Special case this for speed, don't need the waits
        for cmd in jobs:
            results[cmd] = subprocess.call(cmd, shell=True)
        return results
    failed = False
    while pending or running:
        # See if any have finished
        for (cmd, process) in running:
            return_code = process.poll()  # non-blocking
            if return_code is not None:
                results[cmd] = return_code
                if return_code:
                    failed = True
        running = [(cmd, process) for (cmd, process) in running if cmd not in results]
        if verbose:
            print(
                "%i jobs pending, %i running, %i completed"
                % (len(pending), len(running), len(results))
            )
        # See if we can start any new threads
        if pending and failed and fast_fail:
            # Don't start any more jobs
            if verbose:
                print("Failed, will not start remaining %i jobs" % len(pending))
            skipped = pending
            pending = []
        while pending and len(running) < threads:
            cmd = pending.pop(0)
            if verbose:
                print(cmd)
            process = subprocess.Popen(cmd, shell=True)
            running.append((cmd, process))
        # Loop...
        sleep(pause)
    if verbose:
        print("%i jobs completed" % len(results))
    assert set(jobs) == set(results).union(skipped)
    return results
