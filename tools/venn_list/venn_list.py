#!/usr/bin/env python

"""
Plot up to 3-way Venn Diagram.

It uses the Python libraries matplotlib and matplotlib_venn.

This script is copyright 2010-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI), UK. All rights reserved.

See accompanying text file for licence details (MIT License).
"""

from __future__ import print_function

import sys
from shutil import move

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.1.2")
    sys.exit(0)

try:
    import matplotlib
except ImportError:
    sys.exit("Requires the Python library matplotlib")

matplotlib.use("Agg")
from matplotlib import pyplot as plt

if len(sys.argv) - 1 not in [7, 10, 13]:
    sys.exit(
        "Expected 7, 10 or 13 arguments (for 1, 2 or 3 sets), not %i"
        % (len(sys.argv) - 1)
    )

all_file, all_type, all_label = sys.argv[1:4]
set_data = []
if len(sys.argv) - 1 >= 7:
    set_data.append(tuple(sys.argv[4:7]))
if len(sys.argv) - 1 >= 10:
    set_data.append(tuple(sys.argv[7:10]))
if len(sys.argv) - 1 >= 13:
    set_data.append(tuple(sys.argv[10:13]))
pdf_file = sys.argv[-1]
n = len(set_data)
print("Doing %i-way Venn Diagram" % n)


def load_ids(filename, filetype):
    """Load ids from files."""
    if filetype == "tabular":
        for line in open(filename):
            line = line.rstrip("\n")
            if line and not line.startswith("#"):
                yield line.split("\t", 1)[0]
    elif filetype == "fasta":
        for line in open(filename):
            if line.startswith(">"):
                yield line[1:].rstrip("\n").split(None, 1)[0]
    elif filetype.startswith("fastq"):
        # Use the Galaxy library not Biopython to cope with CS
        from galaxy_utils.sequence.fastq import fastqReader

        handle = open(filename)
        for record in fastqReader(handle):
            # The [1:] is because the fastaReader leaves the @ on the identifer.
            yield record.identifier.split()[0][1:]
        handle.close()
    elif filetype == "sff":
        try:
            from Bio.SeqIO import index
        except ImportError:
            sys.exit("Require Biopython 1.54 or later (to read SFF files)")
        # This will read the SFF index block if present (very fast)
        for this_name in index(filename, "sff"):
            yield this_name
    else:
        sys.exit("Unexpected file type %s" % filetype)


def load_ids_whitelist(filename, filetype, whitelist):
    """Check if ids are in whitelist."""
    for single_name in load_ids(filename, filetype):
        if single_name in whitelist:
            yield single_name
        else:
            sys.exit(
                "Unexpected ID %s in %s file %s" % (single_name, filetype, filename)
            )


if all_file in ["", "-", '""', '"-"']:
    # Load without white list
    sets = [set(load_ids(f, t)) for (f, t, c) in set_data]
    # Take union
    all_ids = set()
    for s in sets:
        all_ids.update(s)
    print("Inferred total of %i IDs" % len(all_ids))
else:
    all_ids = set(load_ids(all_file, all_type))
    print("Total of %i IDs" % len(all_ids))
    sets = [set(load_ids_whitelist(f, t, all_ids)) for (f, t, c) in set_data]

for s, (f, t, c) in zip(sets, set_data):
    print("%i in %s" % (len(s), c))

names = [one_name for one_file, one_type, one_name in set_data]
lengths = [len(one_set) for one_set in sets]

if len(sets) == 3:
    try:
        from matplotlib_venn import venn3_unweighted
    except ImportError:
        sys.exit("Requires the Python library matplotlib_venn")
    venn3_unweighted(
        sets,
        [
            "{} (Total {})".format(name, length)
            for (name, length) in zip(names, lengths)
        ],
    )

if len(sets) == 2:
    try:
        from matplotlib_venn import venn2_unweighted
    except ImportError:
        sys.exit("Requires the Python library matplotlib_venn")
    venn2_unweighted(
        sets,
        [
            "{} (Total {})".format(name, length)
            for (name, length) in zip(names, lengths)
        ],
    )

# not sure what I am doing here,
# matplotlib_venn does not want to create a single Venn circle
# stick to the old behavior (rpy and Limma) as much as possible
if len(sets) == 1:
    try:
        from matplotlib_venn import venn2
    except ImportError:
        sys.exit("Requires the Python library matplotlib_venn")
    venn2((sets[0], set()), [set_data[0][2], ""])

plt.title(all_label)
plt.savefig(pdf_file + ".pdf")

# Format "dat" is not supported.
try:
    move(pdf_file + ".pdf", pdf_file)
except (OSError, IOError) as error:
    sys.exit("Fail to rename file {}".format(str(error)))

plt.close()

print("Done")
