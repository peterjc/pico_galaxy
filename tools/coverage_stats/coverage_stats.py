#!/usr/bin/env python
"""BAM coverage statistics using samtools idxstats and depth.

This script takes exactly three command line arguments:
 * Input BAM filename
 * Input BAI filename (via Galaxy metadata)
 * Output tabular filename

This messes about with the filenames to make samtools happy, then
runs "samtools idxstats" and "samtools depth", captures and combines
the output to the desired output tabular file.
"""
import sys
import os
import subprocess
import tempfile

if "-v" in sys.argv or "--version" in sys.argv:
    #Galaxy seems to invert the order of the two lines
    print("BAM coverage statistics v0.0.1")
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

def stop_err(msg, error_level=1):
   """Print error message to stdout and quit with given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)

if len(sys.argv) != 4:
   stop_err("Require three arguments: BAM, BAI, tabular filenames")

bam_filename, bai_filename, tabular_filename = sys.argv[1:]

if not os.path.isfile(bam_filename):
    stop_err("Input BAM file not found: %s" % bam_filename)
if not os.path.isfile(bai_filename):
    if bai_filename == "None":
        stop_err("Error: Galaxy did not index your BAM file")
    stop_err("Input BAI file not found: %s" % bai_filename)

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

def clean_up():
    os.remove(idxstats_filename)
    os.remove(depth_filename)
    os.remove(bam_file)
    os.remove(bai_file)
    os.rmdir(tmp_dir)

# Run samtools idxstats:
cmd = 'samtools idxstats "%s" > "%s"' % (bam_file, idxstats_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    stop_err("Return code %i from command:\n%s" % (return_code, cmd))

# Run samtools depth:
# TODO - Parse stdout instead?
cmd = 'samtools depth "%s" > "%s"' % (bam_file, depth_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    stop_err("Return code %i from command:\n%s" % (return_code, cmd))

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
        # Right at start of file!
        line = depth_handle.readline()
        depth_ref, depth_pos, depth_reads = line.rstrip("\n").split()
        depth_pos = int(depth_pos)
        depth_reads = int(depth_reads)
        # Can now treat as later references where first line cached
    elif identifier != depth_ref:
        # Infer that identifier had coverage zero,
        # and so was not in the 'samtools depth'
        # output.
        # print("%s appears to have no coverage at all" % identifier)
        return 0, 0, 0.0

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
        depth = int(depth)
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
    return min_cov, max_cov, mean_cov

# Parse and combine the output
out_handle = open(tabular_filename, "w")
out_handle.write("#identifer\tlength\tmapped\tplaced\tmin_cov\tmax_cov\tmean_cov\n")

idxstats_handle = open(idxstats_filename)
depth_handle = open(depth_filename)

depth_ref = None
depth_pos = 0
depth_reads = 0

for line in idxstats_handle:
    identifier, length, mapped, placed = line.rstrip("\n").split()
    length = int(length)
    mapped = int(mapped)
    placed = int(placed)
    if identifier == "*":
        min_cov = 0
        max_cov = 0
        mean_cov = 0.0
    else:
        min_cov, max_cov, mean_cov = load_total_coverage(depth_handle, identifier, length)
    out_handle.write("%s\t%i\t%i\t%i\t%i\t%i\t%0.2f\n"
                     % (identifier, length, mapped, placed,
                        min_cov, max_cov, mean_cov))
    if not (min_cov <= mean_cov <= max_cov):
        out_handle.write("ERROR, expect min_cov <= mean_cov <= max_cov\n")
        idxstats_handle.close()
        depth_handle.close()
        out_handle.close()
        clean_up()
        stop_err("Problem with coverage for %s, expect min_cov <= mean_cov <= max_cov"
                 % identifier)

idxstats_handle.close()
depth_handle.close()
out_handle.close()

# Remove the temp symlinks and files:
clean_up()

if depth_ref is not None:
    stop_err("Left over output from 'samtools depth'? %r" % depth_ref)
