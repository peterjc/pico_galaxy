#!/usr/bin/env python
"""Plot Chromosome Diagram using Biopython and ReportLab

This script is copyright 2013 by Peter Cock, The James Hutton Institute
(formerly SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script.
"""
import re
import sys

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

try:
    from Bio import SeqIO
    from Bio.Graphics import BasicChromosome
except ImportError:
    stop_err("Requires the Python library Biopython")

try:
    from reportlab.graphics.shapes import Drawing, String, Line, Rect, Wedge
    from reportlab.lib.units import cm
    from reportlab.lib import colors
except:
    stop_err("Requires the Python library ReportLab (for graphical output)")

if len(sys.argv) != 6:
    stop_err("Expected 5 arguments, not %i" % (len(sys.argv)-1))

ref_file, min_gap, tab_file, main_caption, pdf_file = sys.argv[1:]

#Load reference identifiers and their lengths
#TODO - Support reference aliases? e.g. map 'Chr12' -> '12' or 'XII'
try:
    min_gap = int(min_gap)
except ValueError:
    min_gap = None
if min_gap:
    print "Identifying NNNN regions of at least %i" % min_gap
    re_gap = re.compile("N{%i,}" % min_gap)
    refs = []
    n_regions = []
    #TODO - Speed this up (original script would cache this):
    for rec in SeqIO.parse(ref_file, "fasta"):
        refs.append((rec.id, len(rec)))
        #Look for gaps
        seq = str(rec.seq).upper()
        for match in re_gap.finditer(seq):
            n_regions.append((rec.id, match.start, match.end))
else:
    refs = [(rec.id, len(rec)) for rec in SeqIO.parse(ref_file, "fasta")]
    n_regions = []
lengths = dict(refs)
max_length = max(lengths.values())
print "%i chromosomes/references, max length %i" % (len(refs), max_length)

chr_diagram = BasicChromosome.Organism()
#Automate the size and fonts - or make the user pick?
if True:
    chr_diagram.page_size = (29.7*cm, 21*cm)
    chr_percentage = 0.15
    label_percentage = 0.135
    label_size = 5
elif False:
    chr_diagram.page_size = (30*cm, 22*cm)
    chr_percentage = 0.11
    label_percentage = 0.09
    label_size = 6
else:
    chr_diagram.page_size = (33*cm, 24*cm)
    chr_percentage = 0.11
    label_percentage = 0.12
    label_size = 6
chr_diagram._legend_height = 0

telomere_length = 0.01 * max_length
for name, length in refs:
    caption = name
    features = []

    #Add the N-regions
    for n, start, end in n_regions:
        if n==name:
            features.append((start, end, None, "", colors.black, colors.lightgrey))
    
    cur_chromosome = BasicChromosome.Chromosome(caption)
    #Set the length, adding and extra percentage for the tolomeres & spacers:
    cur_chromosome.scale_num = max_length * 1.02
    cur_chromosome.chr_percent = chr_percentage
    cur_chromosome.label_sep_percent = label_percentage
    cur_chromosome.label_size = label_size
    cur_chromosome._color_labels = True

    #Add an opening spacer (to center all chromosomes vertically)
    space = BasicChromosome.SpacerSegment()
    space.scale = (cur_chromosome.scale_num - length) * 0.5 - telomere_length
    space.chr_percent = chr_percentage
    cur_chromosome.add(space)
   
    #Add an opening telomere
    start = BasicChromosome.TelomereSegment()
    start.scale = telomere_length
    start.chr_percent = chr_percentage
    cur_chromosome.add(start)

    #Add a body - using bp as the scale length here.
    body = BasicChromosome.AnnotatedChromosomeSegment(length, features, colors.blue)
    body.scale = length
    body.chr_percent = chr_percentage
    cur_chromosome.add(body)

    #Add a closing telomere
    end = BasicChromosome.TelomereSegment(inverted=True)
    end.scale = telomere_length
    end.chr_percent = chr_percentage
    cur_chromosome.add(end)

    #Add an closing spacer
    space = BasicChromosome.SpacerSegment()
    space.scale = (cur_chromosome.scale_num - length) * 0.5 - telomere_length
    space.chr_percent = chr_percentage
    cur_chromosome.add(space)

    #This chromosome is done
    chr_diagram.add(cur_chromosome)

    print name, caption, length, len(features)

chr_diagram.draw(pdf_file, main_caption)
print "Done"
