#!/usr/bin/env python
"""Plot up to 3-way Venn Diagram using R limma vennDiagram (via rpy)

This script is copyright 2010-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI), UK. All rights reserved.

See accompanying text file for licence details (MIT License).
"""

from __future__ import print_function

import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.13")
    sys.exit(0)

try:
    import rpy
except ImportError:
    sys.exit("Requires the Python library rpy (to call R)")
except RuntimeError as err:
    sys.exit("The Python library rpy is not availble for the current R version\n\n%s" % err)

try:
    rpy.r.library("limma")
except Exception:
    sys.exit("Requires the R library limma (for vennDiagram function)")


if len(sys.argv) - 1 not in [7, 10, 13]:
    sys.exit("Expected 7, 10 or 13 arguments (for 1, 2 or 3 sets), not %i" % (len(sys.argv) - 1))

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
        handle = open(filename, "rU")
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
        for name in index(filename, "sff"):
            yield name
    else:
        sys.exit("Unexpected file type %s" % filetype)


def load_ids_whitelist(filename, filetype, whitelist):
    for name in load_ids(filename, filetype):
        if name in whitelist:
            yield name
        else:
            sys.exit("Unexpected ID %s in %s file %s" % (name, filetype, filename))


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

# Now call R library to draw simple Venn diagram
try:
    # Create dummy Venn diagram counts object for three groups
    cols = 'c("%s")' % '","'.join("Set%i" % (i + 1) for i in range(n))
    rpy.r('groups <- cbind(%s)' % ','.join(['1'] * n))
    rpy.r('colnames(groups) <- %s' % cols)
    rpy.r('vc <- vennCounts(groups)')
    # Populate the 2^n classes with real counts
    # Don't make any assumptions about the class order
    # print(rpy.r('vc'))
    for index, row in enumerate(rpy.r('vc[,%s]' % cols)):
        if isinstance(row, (int, float)):
            # Hack for rpy being too clever for single element row
            row = [row]
        names = all_ids
        for wanted, s in zip(row, sets):
            if wanted:
                names = names.intersection(s)
            else:
                names = names.difference(s)
        rpy.r('vc[%i,"Counts"] <- %i' % (index + 1, len(names)))
    # print(rpy.r('vc'))
    if n == 1:
        # Single circle, don't need to add (Total XXX) line
        names = [c for (t, f, c) in set_data]
    else:
        names = ["%s\n(Total %i)" % (c, len(s)) for s, (f, t, c) in zip(sets, set_data)]
    rpy.r.assign("names", names)
    rpy.r.assign("colors", ["red", "green", "blue"][:n])
    rpy.r.pdf(pdf_file, 8, 8)
    rpy.r("""vennDiagram(vc, include="both", names=names,
                         main="%s", sub="(Total %i)",
                         circle.col=colors)
                         """ % (all_label, len(all_ids)))
    rpy.r.dev_off()
except Exception as err:
    sys.exit("%s" % str(err))
rpy.r.quit(save="no")
print("Done")
