#!/usr/bin/env python
"""Count sequence variants in region of interest in a BAM file.

This script takes exactly four command line arguments:
 * Input BAM filename
 * Input BAI filename (via Galaxy metadata)
 * Output tabular filename
 * Region of interest (reference:start-end as used in samtools)

This messes about with the filenames to make samtools happy, then
runs "samtools view" and parses the reads mapped to the ROI, counts
the observed variants spanning the ROI, and outputs this as a
tabular file.
"""
import sys
import os
import subprocess
import tempfile

if "-v" in sys.argv or "--version" in sys.argv:
    # Galaxy seems to invert the order of the two lines
    print("BAM coverage statistics v0.0.1")
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

# TODO - Proper command line API
if len(sys.argv) == 5:
    bam_filename, bai_filename, tabular_filename, region = sys.argv[1:]
else:
    sys.exit("Require 4 arguments: BAM, BAI, tabular filename, samtools-style region")

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
    #Sanity check we have "ref:start-end" to give clear error message
    #Note can have semi-colon in the reference name
    #Note can have thousand separator commas in the start/end
    ref, start_end = region.rsplit(":", 1)
    start, end = start_end.split("-")
    start = int(start.replace(",", ""))
    end = int(end.replace(",", ""))
except ValueError:
    sys.exit("Bad argument for region: %r" % region)


# Assign sensible names with real extensions, and setup symlinks:
tmp_dir = tempfile.mkdtemp()
bam_file = os.path.join(tmp_dir, "temp.bam")
bai_file = os.path.join(tmp_dir, "temp.bam.bai")
os.symlink(os.path.abspath(bam_filename), bam_file)
os.symlink(os.path.abspath(bai_filename), bai_file)
assert os.path.isfile(bam_file), bam_file
assert os.path.isfile(bai_file), bai_file
assert os.path.isfile(bam_file + ".bai"), bam_file


def clean_up():
    os.remove(bam_file)
    os.remove(bai_file)
    os.rmdir(tmp_dir)


def decode_cigar(cigar):
    """Returns a list of 2-tuples, integer count and operator char."""
    count = ""
    answer = []
    for letter in cigar:
        if letter.isdigit():
            count += letter #string addition
        elif letter in "MIDNSHP=X":
            answer.append((int(count), letter))
            count = ""
        else:
            raise ValueError("Invalid character %s in CIGAR %s" % (letter, cigar))
    return answer

assert decode_cigar("14S15M1P1D3P54M1D34M5S") == [(14,'S'),(15,'M'),(1,'P'),(1,'D'),(3,'P'),(54,'M'),(1,'D'),(34,'M'),(5,'S')]

def align_len(cigar_ops):
    """Sums the CIGAR M/=/X/D/N operators."""
    return sum(count for count, op in cigar_ops if op in "M=XDN")

def expand_cigar(seq, cigar_ops):
    """Yields (ref_offset, seq_base) pairs."""
    ref_offset = 0
    seq_offset = 0
    for count, op in cigar_ops:
        if op in "MX=":
            for (i, base) in enumerate(seq[seq_offset:seq_offset+count]):
                yield ref_offset + i, base
            ref_offset += count
            seq_offset += count
        elif op == "I":
            # Give them all an in-between reference position
            # (Python lets us mix integers and floats, wouldn't work in C)
            for (i, base) in enumerate(seq[seq_offset:seq_offset+count]):
                yield ref_offset - 0.5, base
            # Does not change ref_offset
            seq_offset += count
        elif op in "DN":
            # Deletion/skip, return gap characters (OK?)
            for i in range(count):
                yield ref_offset + i, "-"
            ref_offset += count
        elif op == "S":
            # Soft clipping, silently discard the bases (OK?)
            seq_offset += count
        elif op in "HP":
            # Hard trimming or pad, can ignore
            pass
        else:
            raise NotImplementedError("Unexpected CIGAR operator %s" % op)

assert list(expand_cigar("ACGT", decode_cigar("4M"))) == [(0, "A"), (1, "C"), (2, "G"), (3, "T")]
assert list(expand_cigar("ACGT", decode_cigar("2=1X1="))) == [(0, "A"), (1, "C"), (2, "G"), (3, "T")]
assert list(expand_cigar("ACGT", decode_cigar("2M1D2M"))) == [(0, "A"), (1, "C"), (2, "-"), (3, "G"), (4, "T")]
assert list(expand_cigar("ACtGT", decode_cigar("2M1I2M"))) == [(0, "A"), (1, "C"), (1.5, "t"), (2, "G"), (3, "T")]
assert list(expand_cigar("tACGT", decode_cigar("1I4M"))) == [(-0.5, 't'), (0, 'A'), (1, 'C'), (2, 'G'), (3, 'T')]
assert list(expand_cigar("ACGTt", decode_cigar("4M1I"))) == [(0, 'A'), (1, 'C'), (2, 'G'), (3, 'T'), (3.5, 't')]
assert list(expand_cigar("AAAAGGGGTTTT", decode_cigar("12M"))) == [(0, 'A'), (1, 'A'), (2, 'A'), (3, 'A'), (4, 'G'), (5, 'G\
'), (6, 'G'), (7, 'G'), (8, 'T'), (9, 'T'), (10, 'T'), (11, 'T')]
assert list(expand_cigar("AAAAcGGGGTTTT", decode_cigar("4M1I8M"))) == [(0, 'A'), (1, 'A'), (2, 'A'), (3, 'A'), (3.5, 'c'), (4, 'G'), (5, 'G'), (6, 'G'), (7, 'G'), (8, 'T'), (9, 'T'), (10, 'T'), (11, 'T')]
assert list(expand_cigar("AAAAGGGGcTTTT", decode_cigar("8M1I4M"))) == [(0, 'A'), (1, 'A'), (2, 'A'), (3, 'A'), (4, 'G'), (5, 'G\
'), (6, 'G'), (7, 'G'), (7.5, "c"), (8, 'T'), (9, 'T'), (10, 'T'), (11, 'T')]
assert list(expand_cigar("AAAAcGGGGcTTTT", decode_cigar("4M1I4M1I4M"))) == [(0, 'A'), (1, 'A'), (2, 'A'), (3, 'A'), (3.5, 'c'), (4, 'G'), (5, 'G'), (6, 'G'), (7, 'G'), (7.5, 'c'), (8, 'T'), (9, 'T'), (10, 'T'), (11, 'T')]

def get_roi(seq, cigar_ops, start, end):
    """Extract region of seq mapping to the ROI.

    Expect start and end to be zero based Python style end points.

    i.e. The ROI relative to the mapping start recorded in the POS field.
    Will return part of the SAM/BAM value SEQ based on interpretting the
    passed CIGAR operators.
    """
    if len(cigar_ops) == 1 and cigar_ops[0][1] in "M=X":
        # Easy case, note start/end/pos all one-based
        assert cigar_ops[0][0] == len(seq)
        return seq[start:end]
    # Would use "start <= i < end" if they were all integers, but
    # want to exclude e.g. 3.5 and 7.5 when given start 4 and end 8. 
    return "".join(base for i, base in expand_cigar(seq, cigar_ops) if start <= i <= end-1)

assert "GGGG" == get_roi("AAAAGGGGTTTT", decode_cigar("12M"), 4, 8)
assert "GGGG" == get_roi("AAAAcGGGGTTTT", decode_cigar("4M1I8M"), 4, 8)
assert "GGGG" == get_roi("AAAAGGGGcTTTT", decode_cigar("8M1I4M"), 4, 8)
assert "GGGG" == get_roi("AAAAcGGGGcTTTT", decode_cigar("4M1I4M1I4M"), 4, 8)
assert "GGaGG" == get_roi("AAAAGGaGGTTTT", decode_cigar("6M1I6M"), 4, 8)
assert "GGGgA" == get_roi("AAAAGGGgATTTT", decode_cigar("7M1I5M"), 4, 8)

def count_region():
    # Could recreate the region string (with no commas in start/end)?
    # region = "%s:%i-%i" % (ref, start, end)

    tally = dict()

    # Call samtools view, don't need header so no -h added.
    # Only want mapped reads, thus flag filter -F 4.
    child = subprocess.Popen(["samtools", "view", "-F", "4", bam_file, region],
                             stdout=subprocess.PIPE)
    for line in child.stdout:
        assert line[0] != "@", "Got unexpected SAM header line: %s" % line
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, rest = line.split("\t", 10)
        pos = int(pos)  # one-based
        if start < pos:
            # Does not span the ROI
            continue
        cigar_ops = decode_cigar(cigar)
        if pos + align_len(cigar_ops) - 1 < end:
            # Does not span the ROI
            continue
        # All of start/end/pos are currently one-based, making offsets Python style....
        roi_seq = get_roi(seq, cigar_ops, start - pos, end - pos + 1)
        assert roi_seq, "Error, empty ROI sequence for: %s" % line
        try:
            tally[roi_seq] += 1
        except KeyError:
            tally[roi_seq] = 1

    child.stdout.close()
    return_code = child.wait()
    if return_code:
        sys.exit("Got return code %i from samtools view" % return_code)

    return tally

def record_counts():

    tally = count_region()
    total = sum(tally.values())

    # Using negative count to get sort with highest count first,
    # while tie-breaking by the ROI sequence alphabetically.
    table = sorted((-count, roi_seq) for (roi_seq, count) in tally.items())
    del tally

    with open(tabular_filename, "w") as handle:
        handle.write("Variant\tCount\tPercentage\n")
        for count, roi_seq in table:
            handle.write("%s\t%i\t%0.2f\n" % (roi_seq, -count, -count * 100.0 / total))

    print("Counted %i variants from %i reads spanning %s" % (len(table), total, region))


# Run it!
record_counts()
# Remove the temp symlinks and files:
clean_up()
