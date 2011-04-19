#!/usr/bin/env python
"""Plot up to 3-way Venn Diagram using R limma vennDiagram (via rpy)

This script is copyright 2010 by Peter Cock, The James Hutton Institute
(formerly SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.2 of the script.
"""


import sys
import rpy

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

try:
    import rpy
except ImportError:
    stop_err("Requires the Python library rpy (to call R)")

try:
    rpy.r.library("limma")
except:
    stop_err("Requires the R library limma (for vennDiagram function)")


if len(sys.argv)-1 not in [7, 10, 13]:
    stop_err("Expected 7, 10 or 13 arguments (for 1, 2 or 3 sets), not %i" % (len(sys.argv)-1))

all_file, all_type, all_label = sys.argv[1:4]
set_data = []
if len(sys.argv)-1 >= 7:
    set_data.append(tuple(sys.argv[4:7]))
if len(sys.argv)-1 >= 10:
    set_data.append(tuple(sys.argv[7:10]))
if len(sys.argv)-1 >= 13:
    set_data.append(tuple(sys.argv[10:13]))
pdf_file = sys.argv[-1]
print "Doing %i-way Venn Diagram" % len(set_data)

def load_ids(filename, filetype):
    if filetype=="tabular":
        for line in open(filename):
            if not line.startswith("#"):
                yield line.rstrip("\n").split("\t",1)[0]
    elif filetype=="fasta":
        for line in open(filename):
            if line.startswith(">"):
                yield line[1:].rstrip("\n").split(None,1)[0]
    elif filetype.startswith("fastq"):
        #Use the Galaxy library not Biopython to cope with CS
        from galaxy_utils.sequence.fastq import fastqReader
        handle = open(filename, "rU")
        for record in fastqReader(handle):
            #The [1:] is because the fastaReader leaves the @ on the identifer.
            yield record.identifier.split()[0][1:]
        handle.close()
    elif filename=="sff":
        try:
            from Bio.SeqIO.SffIO import SffIterator
        except ImportError:
            stop_err("Require Biopython 1.54 or later (to read SFF files)")
        #TODO - Try and load the index (if present), much faster!
        handle = open(filename, "rb")
        for record in SffIterator(handle):
            yield name
        handle.close()
    else:
        stop_err("Unexpected file type %s" % filetype)

def load_ids_whitelist(filename, filetype, whitelist):
    for name in load_ids(filename, filetype):
        if name in whitelist:
            yield name
        else:
            stop_err("Unexpected ID %s in %s file %s" % (name, filetype, filename))

all = set(load_ids(all_file, all_type))
print "Total of %i IDs" % len(all)
sets = [set(load_ids_whitelist(f,t,all)) for (f,t,c) in set_data]

for s, (f,t,c) in zip(sets, set_data):
    print "%i in %s" % (len(s), c)


A,B,C = sets
a_label, b_label, c_label = "A", "B", "C"

#Now call R library to draw simple Venn diagram
try:
    #Create dummy Venn diagram counts object for three groups
    rpy.r('groups <- cbind(1,1,1)')
    rpy.r('colnames(groups) <- c("A", "B", "C")')
    rpy.r('vc <- vennCounts(groups)')
    #print rpy.r('vc')
    #Populate the 8 classes with real counts
    #Don't make any assumptions about the class order
    for row,(a,b,c) in enumerate(rpy.r('vc[,c("A","B","C")]')):
        names = all
        if a:
            names = names.intersection(A)
        else:
            names = names.difference(A)
        if b:
            names = names.intersection(B)
        else:
            names = names.difference(B)
        if c:
            names = names.intersection(C)
        else:
            names = names.difference(C)
        #print a,b,c,names
        rpy.r('vc[%i,"Counts"] <- %i' % (row+1, len(names)))
    #print rpy.r('vc')
    #rpy.r.assign("title", "%s\n(Total %i)" % (all_label, len(all)))
    rpy.r.assign("names", ["%s\n(Total %i)" % (a_label, len(A)),
                           "%s\n(Total %i)" % (b_label, len(B)),
                           "%s\n(Total %i)" % (c_label, len(C))])
    rpy.r.pdf(pdf_file, 8, 8)
    rpy.r("""vennDiagram(vc, include="both", names=names,
                         main="%s", sub="(Total %i)",
                         circle.col=c("red","green","blue"))
                         """ % (all_label, len(all)))
    rpy.r.dev_off()
except Exception, exc:
    stop_err( "%s" %str( exc ) )
rpy.r.quit( save="no" )
print "Done"
