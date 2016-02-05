#!/usr/bin/env python
"""BAM coverage statistics using samtools idxstats and depth.

This script takes exactly three command line arguments:
 * Input BAM filename
 * Input BAI filename (via Galaxy metadata)
 * Output tabular filename

Optional fourth argument:
 * Max coverage depth (integer)

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

if "-v" in sys.argv or "--version" in sys.argv:
    # Galaxy seems to invert the order of the two lines
    print("BAM coverage statistics v0.0.5")
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

# TODO - Proper command line API
if len(sys.argv) == 4:
    bam_filename, bai_filename, tabular_filename = sys.argv[1:]
    max_depth = "8000"
elif len(sys.argv) == 5:
    bam_filename, bai_filename, tabular_filename, max_depth = sys.argv[1:]
else:
    sys.exit("Require 3 or 4 arguments: BAM, BAI, tabular filename, [max depth]")

if not os.path.isfile(bam_filename):
    sys.exit("Input BAM file not found: %s" % bam_filename)
if bai_filename == "-":
    # Make this optional for ease of use at the command line by hand:
    if os.path.isfile(bam_filename + ".bai"):
        bai_filename = bam_filename + ".bai"
if not os.path.isfile(bai_filename):
    if bai_filename == "None":
        sys.exit("Error: Galaxy did not index your BAM file")
    sys.exit("Input BAI file not found: %s" % bai_filename)

try:
    max_depth = int(max_depth)
except ValueError:
    sys.exit("Bad argument for max depth: %r" % max_depth)
if max_depth < 0:
    sys.exit("Bad argument for max depth: %r" % max_depth)

# fuzz factor to ensure can reach max_depth, e.g. region with
# many reads having a deletion present could underestimate the
# coverage by capping the number of reads considered
depth_margin = 100

# Assign sensible names with real extensions, and setup symlinks:
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


def clean_up():
    os.remove(idxstats_filename)
    os.remove(depth_filename)
    os.remove(bam_file)
    os.remove(bai_file)
    os.rmdir(tmp_dir)


def samtools_depth_opt_available():
    child = subprocess.Popen(["samtools", "depth"],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Combined stdout/stderr in case samtools is ever inconsistent
    output, tmp = child.communicate()
    assert tmp is None
    # Expect to find this line in the help text, exact wording could change:
    #    -d/-m <int>         maximum coverage depth [8000]
    return " -d/-m " in output

depth_hack = False
if not samtools_depth_opt_available():
    if max_depth + depth_margin <= 8000:
        sys.stderr.write("WARNING: The version of samtools depth installed does not "
                         "support the -d option, however, the requested max-depth "
                         "is safely under the default of 8000.\n")
        depth_hack = True
    else:
        sys.exit("The version of samtools depth installed does not support the -d option.")

# Run samtools idxstats:
cmd = 'samtools idxstats "%s" > "%s"' % (bam_file, idxstats_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    sys.exit("Return code %i from command:\n%s" % (return_code, cmd))

# Run samtools depth:
# TODO - Parse stdout instead?
if depth_hack:
    # Using an old samtools without the -d option, but hard coded default
    # of 8000 should be fine even allowing a margin for fuzzy output
    cmd = 'samtools depth "%s" > "%s"' % (bam_file, depth_filename)
else:
    cmd = 'samtools depth -d %i "%s" > "%s"' % (max_depth + depth_margin, bam_file, depth_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    sys.exit("Return code %i from command:\n%s" % (return_code, cmd))


def load_total_coverage(depth_handle, identifier, length):
    """Parse some of the 'samtools depth' output for coverages.

    Returns min_cov (int), max_cov (int) and mean cov (float).

    Uses global variables to cache the first line of output from the
    next reference sequence.
    """
    global depth_ref, depth_pos, depth_reads

    # print("====")
    # print("%s coverage calculation, length %i, ..." % (identifier, length))

    if depth_ref is None:
        # Right at start of file / new contig
        line = depth_handle.readline()
        # Are we at the end of the file?
        if not line:
            # Must be at the end of the file.
            # This can happen if the file contig(s) had no reads mapped
            return 0, 0, 0.0, 0
        depth_ref, depth_pos, depth_reads = line.rstrip("\n").split()
        depth_pos = int(depth_pos)
        depth_reads = min(max_depth, int(depth_reads))
        # Can now treat as later references where first line cached
    elif identifier != depth_ref:
        # Infer that identifier had coverage zero,
        # and so was not in the 'samtools depth'
        # output.
        # print("%s appears to have no coverage at all" % identifier)
        return 0, 0, 0.0, 0

    # Good, at start of expected reference
    bases = depth_reads
    if depth_pos == 1:
        min_cov = depth_reads
    else:
        # print("%s has no coverage at start" % identifier)
        min_cov = 0
    max_cov = depth_reads

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
        bases += depth
        if last_pos + 1 < pos:
            # print("%s has no coverage between %i and %i" % (identifier, last_pos, pos))
            min_cov = 0
        else:
            min_cov = min(min_cov, depth)
        max_cov = max(max_cov, depth)
        last_pos = pos

    # Reach the end of this identifier's coverage or end of file
    if last_pos < length:
        # print("%s has no coverage at end" % identifier)
        min_cov = 0
    mean_cov = bases / float(length)
    return min_cov, max_cov, mean_cov, bases

# Parse and combine the output
out_handle = open(tabular_filename, "w")
out_handle.write("#identifer\tlength\tmapped\tplaced\tmin_cov\tmax_cov\tmean_cov\n")

idxstats_handle = open(idxstats_filename)
depth_handle = open(depth_filename)

depth_ref = None
depth_pos = 0
depth_reads = 0

global_bases = 0
global_length = 0
for line in idxstats_handle:
    identifier, length, mapped, placed = line.rstrip("\n").split()
    length = int(length)
    mapped = int(mapped)
    placed = int(placed)
    if identifier == "*":
        min_cov = 0
        max_cov = 0
        mean_cov = 0.0
        bases = 0
    else:
        min_cov, max_cov, mean_cov, bases = load_total_coverage(depth_handle, identifier, length)
    if max_cov > max_depth:
        sys.exit("Using max depth %i yet saw max coverage %i for %s"
                 % (max_depth, max_cov, identifier))
    out_handle.write("%s\t%i\t%i\t%i\t%i\t%i\t%0.2f\n"
                     % (identifier, length, mapped, placed,
                        min_cov, max_cov, mean_cov))
    if not (min_cov <= mean_cov <= max_cov):
        out_handle.write("ERROR, expect min_cov <= mean_cov <= max_cov\n")
        idxstats_handle.close()
        depth_handle.close()
        out_handle.close()
        clean_up()
        sys.exit("Problem with coverage for %s, expect min_cov <= mean_cov <= max_cov"
                 % identifier)
    global_length += length
    global_bases += bases

idxstats_handle.close()
depth_handle.close()
out_handle.close()

print("Total reference length %i with overall mean coverage %0.2f" % (global_length, float(global_bases) / global_length))

# Remove the temp symlinks and files:
clean_up()

if depth_ref is not None:
    sys.exit("Left over output from 'samtools depth'? %r" % depth_ref)
