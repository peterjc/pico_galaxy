#!/usr/bin/env python
"""BAM coverage statistics using samtools idxstats and depth.

This script expects exactly five command line arguments:
 * Input BAM filename
 * Input BAI filename (via Galaxy metadata)
 * Output tabular filename
 * Max coverage depth (integer)
 * Sliding window size (integer)

This messes about with the filenames to make samtools happy, then
runs "samtools idxstats" and "samtools depth", captures and combines
the output to the desired output tabular file.

Because "samtools depth" treats the max depth a little fuzzily, this
tool tries to account for this and applies a clear max-depth cut off.
"""
import sys
import os
import subprocess
import tempfile

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # Quick and dirty, assumes already a list not an interator
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

if "-v" in sys.argv or "--version" in sys.argv:
    #Galaxy seems to invert the order of the two lines
    print("BAM coverage statistics v0.1.0")
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

def sys_exit(msg, error_level=1):
   """Print error message to stdout and quit with given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)

try:
    import numpy as np
except ImportError:
    sys.exit("Missing Python dependancy numpy")

# TODO - Proper command line API
if len(sys.argv) == 4:
    bam_filename, bai_filename, tabular_filename = sys.argv[1:]
    max_depth = "8000"
elif len(sys.argv) == 5:
    bam_filename, bai_filename, tabular_filename, max_depth = sys.argv[1:]
else:
    sys_exit("Require 3 or 4 arguments: BAM, BAI, tabular filename, [max depth]")

if not os.path.isfile(bam_filename):
    sys_exit("Input BAM file not found: %s" % bam_filename)
if not os.path.isfile(bai_filename):
    if bai_filename == "None":
        sys_exit("Error: Galaxy did not index your BAM file")
    sys_exit("Input BAI file not found: %s" % bai_filename)

try:
    max_depth = int(max_depth)
except ValueError:
    sys_exit("Bad argument for max depth: %r" % max_depth)
if max_depth < 0:
    sys_exit("Bad argument for max depth: %r" % max_depth)

# fuzz factor to ensure can reach max_depth, e.g. region with
# many reads having a deletion present could underestimate the
# coverage by capping the number of reads considered
depth_margin = 100
wanted_percentiles = [1, 5, 25, 50, 75, 95, 99]

#Assign sensible names with real extensions, and setup symlinks:
tmp_dir = tempfile.mkdtemp()
bam_file = os.path.join(tmp_dir, "temp.bam")
bai_file = os.path.join(tmp_dir, "temp.bam.bai")
idxstats_filename = os.path.join(tmp_dir, "idxstats.tsv")
depth_filename = os.path.join(tmp_dir, "depth.tsv")
os.symlink(os.path.abspath(bam_filename), bam_file)
os.symlink(os.path.abspath(bai_filename), bai_file)
assert os.path.isfile(bam_file), bam_file
assert os.path.isfile(bai_file), bai_file
assert os.path.isfile(bam_file + ".bai"), bam_file

class NoCoverage(Exception):
    pass

def clean_up():
    os.remove(idxstats_filename)
    os.remove(depth_filename)
    os.remove(bam_file)
    os.remove(bai_file)
    os.rmdir(tmp_dir)


def samtools_depth_opt_available():
    child = subprocess.Popen(["samtools", "depth"],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Combined stdout/stderr in case samtools is every inconsistent
    output, tmp = child.communicate()
    assert tmp is None
    # Expect to find this line in the help text, exact wording could change:
    #    -d/-m <int>         maximum coverage depth [8000]
    return " -d/-m " in output

if not samtools_depth_opt_available():
    # TODO - gracefull fall back if user wants max depth under 8000?
    sys_exit("The version of samtools installed does not support the -d option.")

# Run samtools idxstats:
cmd = 'samtools idxstats "%s" > "%s"' % (bam_file, idxstats_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    sys_exit("Return code %i from command:\n%s" % (return_code, cmd))

# Run samtools depth:
# TODO - Parse stdout instead?
cmd = 'samtools depth -d %i "%s" > "%s"' % (max_depth + depth_margin, bam_file, depth_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    sys_exit("Return code %i from command:\n%s" % (return_code, cmd))

def get_stats(depth_iterator, length):
    """Calculate min/max/mean coverage.

    Returns tuple (int, int, float).
    """

    # Plus one as counting from 0 to max_depth inclusive
    tallies = np.zeros([max_depth + 1], np.int)

    # First call to iterator may return NoCoverage
    # This also makes initialising the min value easy
    try:
        ref, pos, depth = next(depth_iterator)
    except NoCoverage:
        return 0, 0, 0.0
    except StopIteration:
        raise ValueError("Internal error - was there no coverage?")
    total_cov = min_cov = max_cov = depth
    assert depth <= max_depth, depth
    tallies[depth] += 1

    # Could check pos is strictly increasing and within 1 to length?
    for ref, pos, depth in depth_iterator:
        total_cov += depth
        min_cov = min(min_cov, depth)
        max_cov = max(max_cov, depth)
        tallies[depth] += 1

    mean_cov = total_cov / float(length)

    assert min_cov <= mean_cov <= max_cov
    assert length == sum(tallies)
    assert total_cov == sum(value*count for value, count in enumerate(tallies))

    return min_cov, max_cov, mean_cov, tallies


def depth_iterator(handle, identifier, length):
    """Iterates over ``samtools depth`` output for coverage depths.

    Uses global variables to cache the first line of output from the
    next reference sequence. Will impute missing zero depth lines
    which ``samtools depth`` omits.
    """
    global depth_ref, depth_pos, depth_reads

    assert identifier

    if depth_ref is None:
        # Right at start of file / new contig
        line = depth_handle.readline()
        if not line:
            # Must be at the end of the file.
            # This can happen if the file contig(s) had no reads mapped
            raise NoCoverage("End of file, no reads mapped to %r" % identifier)
        depth_ref, depth_pos, depth_reads = line.rstrip("\n").split()
        depth_pos = int(depth_pos)
        depth_reads = min(max_depth, int(depth_reads))
        # Can now treat as later references where first line cached
    elif identifier != depth_ref:
        # Infer that identifier had coverage zero,
        # and so was not in the 'samtools depth'
        # output.
        raise NoCoverage("%r appears to have no coverage at all" % identifier)

    # Good, at start of expected reference
    if depth_pos > 1:
        for extra_pos in range(1, depth_pos):
            yield identifier, extra_pos, 0
    yield identifier, depth_pos, depth_reads

    last_pos = depth_pos
    depth_ref = None
    depth_pos = 0
    depth_reads = 0
    for line in depth_handle:
        ref, pos, depth = line.rstrip("\n").split()
        pos = int(pos)
        depth = min(max_depth, int(depth))
        if ref != identifier:
            # Reached the end of this identifier's coverage
            # so cache this ready for next identifier
            depth_ref, depth_pos, depth_reads = ref, pos, depth
            break
        if last_pos + 1 < pos:
            for extra_pos in range(last_pos + 1, pos):
                yield identifier, extra_pos, 0
            # print("%s has no coverage between %i and %i" % (identifier, last_pos, pos))
        yield identifier, pos, depth
        last_pos = pos

    #Reached the end of this identifier's coverage or end of file
    if last_pos < length:
        for extra_pos in range(last_pos + 1, length +1):
            yield identifier, extra_pos, 0

    # Function ready to be called again to iterator over next ref


def percentiles_from_counts(tally_table, percentiles):
    """Calculate percentails from a tally table.

    The tally_table should be a list of counts as described below,
    percentiles should be an *increasing* list of percentile points
    desired described as integers between 0 and 100.

    If your data is simply in a list or array, use np.percentile,

    >>> data = [0, 5, 6, 4, 3, 3, 1, 0, 1] # positive or zero integers!
    >>> percentiles = [25, 50, 75]
    >>> import numpy as np
    >>> [np.percentile(data, p) for p in percentiles]
    [1.0, 3.0, 4.0]

    However, especially if dealing with large datasets, you may have
    constructed a simple tally table by iterating over the dataset
    (e.g. while parsing a large file):

    >>> assert 0 <= min(data)
    >>> tally = np.zeros([max(data) + 1], np.int)
    >>> for value in data:
    ...     tally[value] += 1
    ...
    >>> assert sum(tally) == len(data)
    >>> tally
    array([2, 2, 0, 2, 1, 1, 1])

    This kind of tally table array (or a list equivalent) makes sense if
    recording the frequencies of integers in the data set. Now use this
    function:

    >>> percentiles_from_counts(tally, percentiles)
    (1.0, 3.0, 4.0)
    """

    N = sum(list(tally_table))
    assert isinstance(N, int), N

    # TODO: Refactor to avoid resetting for each p
    answer = []
    for p in percentiles:
        target = p / 100.0 * N
        cumulative = 0.0
        # Python 3 TODO: Want to leak count variable from loop...
        for value, count in enumerate(tally_table):
            cumulative += count
            if target < cumulative:
                break
        answer.append(value)
    return answer

tmp = percentiles_from_counts([2, 2, 0, 2, 1, 1, 1], [25, 50, 75])
assert [1.0, 3.0, 4.0] == tmp, tmp
del tmp
        

def load_total_coverage(depth_handle, identifier, length):
    """Parse some of the 'samtools depth' output for coverages.

    Returns min_cov (int), max_cov (int), mean cov (float),
    assorted percentials (float).

    Uses global variables to cache the first line of output from the
    next reference sequence.

    Uses global variables for the max-depth and window size.
    """
    global depth_ref, depth_pos, depth_reads

    assert identifier
    assert 0 < length

    depth_iter = depth_iterator(depth_handle, identifier, length)
    min_cov, max_cov, mean_cov, tallies = get_stats(depth_iter, length)
    percentiles = percentiles_from_counts(tallies, wanted_percentiles)

    return min_cov, max_cov, mean_cov, percentiles

# Parse and combine the output
out_handle = open(tabular_filename, "w")
# TODO - extra headers
p_str = "\t".join("P%i" % p for p in wanted_percentiles)
out_handle.write("#identifer\tlength\tmapped\tplaced\tmin_cov\tmax_cov\tmean_cov\t%s\n" % p_str)

idxstats_handle = open(idxstats_filename)
depth_handle = open(depth_filename)

depth_ref = None
depth_pos = 0
depth_reads = 0
err = None
for line in idxstats_handle:
    identifier, length, mapped, placed = line.rstrip("\n").split()
    length = int(length)
    mapped = int(mapped)
    placed = int(placed)
    if identifier == "*":
        min_cov = 0
        max_cov = 0
        mean_cov = 0.0
        percentiles = [0.0] * len(wanted_percentiles)
    else:
        min_cov, max_cov, mean_cov, percentiles \
            = load_total_coverage(depth_handle, identifier, length)
    p_str = "\t".join("%0.2f" % p for p in percentiles)
    out_handle.write("%s\t%i\t%i\t%i\t%i\t%i\t%0.2f\t%s\n"
                     % (identifier, length, mapped, placed,
                        min_cov, max_cov, mean_cov, p_str))
    if max_cov > max_depth:
        err = ("Using max depth %i yet saw max coverage %i for %s"
              % (max_depth, max_cov, identifier))
    if not (min_cov <= mean_cov <= max_cov):
        err = ("Min/mean/max inconsistent: "
               "expect %r <= %r < %r "
               % (min_cov, mean_cov, max_cov))
    if err:
        out_handle.write("ERROR during %s: %s\n" % (identifier, err))
        idxstats_handle.close()
        depth_handle.close()
        out_handle.close()
        clean_up()
        sys_exit("ERROR during %s: %s\n" % (identifier, err))

idxstats_handle.close()
depth_handle.close()
out_handle.close()

# Remove the temp symlinks and files:
clean_up()

if depth_ref is not None:
    sys_exit("Left over output from 'samtools depth'? %r" % depth_ref)
