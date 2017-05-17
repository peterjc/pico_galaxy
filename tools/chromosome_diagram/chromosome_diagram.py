#!/usr/bin/env python
"""Plot Chromosome Diagram using Biopython and ReportLab

This script is copyright 2013 by Peter Cock, The James Hutton Institute
(formerly SCRI), UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script.
"""

from __future__ import print_function

import re
import sys
import warnings

try:
    from Bio import SeqIO
    from Bio.Graphics import BasicChromosome
    from Bio import BiopythonWarning
except ImportError:
    sys.exit("Requires the Python library Biopython")

try:
    from reportlab.pdfgen import canvas
    from reportlab.graphics import renderPDF
    from reportlab.graphics.shapes import Drawing, String, Line, Rect
    # We don't use these directly, Biopython will though so checking now
    del Drawing, String, Line, Rect
    from reportlab.lib.units import cm, inch
    from reportlab.lib import colors
except Exception:
    sys.exit("Requires the Python library ReportLab (for graphical output)")

if len(sys.argv) - 1 != 13:
    sys.exit("Expected 13 arguments, not %i" % (len(sys.argv) - 1))

ref_file, min_gap, tab_file, chr_col, start_col, end_col, strand_col, \
    caption_col, color_col, fill_col, main_caption, per_page, pdf_file = sys.argv[1:]


def load_column(txt):
    try:
        return int(txt) - 1
    except ValueError:
        if txt.strip():
            sys.exit("Bad column argument %r." % txt)
        return None


chr_col = load_column(chr_col)
start_col = load_column(start_col)
end_col = load_column(end_col)
strand_col = load_column(strand_col)
caption_col = load_column(caption_col)
color_col = load_column(color_col)
fill_col = load_column(fill_col)

try:
    per_page = int(per_page.strip())
except ValueError:
    if per_page.string():
        sys.exit("Bad per-page argument %r" % per_page)
    per_page = 0  # All on one page

# Load reference identifiers and their lengths
# TODO - Support reference aliases? e.g. map 'Chr12' -> '12' or 'XII'
try:
    min_gap = int(min_gap)
except ValueError:
    min_gap = None
if min_gap:
    print("Identifying long NNNN regions of at least %i" % min_gap)
    big_gap = "N" * min_gap
    re_not_gap = re.compile("[^N]")
    refs = []
    n_regions = []
    for rec in SeqIO.parse(ref_file, "fasta"):
        refs.append((rec.id, len(rec)))
        seq = str(rec.seq).upper()
        length = len(seq)
        start = 0
        while True:
            start = seq.find(big_gap, start)
            if start == -1:
                break
            match = re_not_gap.search(seq, start)
            if not match:
                end = length
                # print "%s %i - end (%i)" % (rec.id, start, end-start)
                # assert end-start >= min_gap
                n_regions.append((rec.id, start, end))
                break
            end = match.start()
            # print "%s %i - %i (%i)" % (rec.id, start, end, end-start)
            # assert end-start >=min_gap
            n_regions.append((rec.id, start, end))
            start = end
else:
    refs = [(rec.id, len(rec)) for rec in SeqIO.parse(ref_file, "fasta")]
    n_regions = []
lengths = dict(refs)
max_length = max(lengths.values())
if not per_page:
    per_page = len(refs)
print("%i chromosomes/references, max length %i" % (len(refs), max_length))


def load_color(txt, default=None):
    txt = txt.strip()
    if txt.lower() in ["", ".", "?", "none"]:
        return default
    elif len(txt) == 6 and set("0123456789ABCDEF").issuperset(txt.upper()):
        # Hex color
        return colors.HexColor("#%s" % txt)
    elif len(txt) == 7 and txt[0] == "#" and set("0123456789ABCDEF").issuperset(txt[1:].upper()):
        # Hex color with # prefix
        return colors.HexColor(txt)
    elif len(txt) == 8 and txt[0] == "#" and txt[-1] == ";" and set("0123456789ABCDEF").issuperset(txt[1:-1].upper()):
        # Hex color with # prefix and ; suffix
        return colors.HexColor(txt[:-1])
    else:
        # Let Biopython deal with it...
        return txt


# Load the features
all_features = []
handle = open(tab_file)
for line in handle:
    if line.startswith("#"):
        continue
    parts = line.rstrip("\n").split("\t")
    chr_name = parts[chr_col].strip()
    if chr_name not in lengths:
        sys.exit("Unrecognised reference/chromosome %r in this line:\n%r" % (parts[chr_col], line))
    start = int(parts[start_col])
    if not (0 <= start <= lengths[chr_name]):
        sys.exit("Start %i outside length %i of reference %s in this line:\n%r"
                 % (start, lengths[chr_name], chr_name, line))
    if end_col is None:
        end = start
    else:
        end = int(parts[end_col])
    if not (0 <= end <= lengths[chr_name]):
        sys.exit("End %i outside length %i of reference %s in this line:\n%r"
                 % (end, lengths[chr_name], chr_name, line))
    if strand_col is None:
        strand = None
    else:
        strand = parts[strand_col].strip().lower()
        if strand in ["+1", "1", "+", "f", "forward"]:
            strand = +1
        elif strand in ["-1", "-", "r", "reverse"]:
            strand = -1
        elif strand in ["0", "?", ".", "none", "both", ""]:
            strand = None
        else:
            sys.exit("Bad strand value %r in this line:\n%r" % (parts[strand_col], line))
    caption = parts[caption_col]
    if color_col is None:
        color = colors.black
    else:
        color = load_color(parts[color_col], colors.black)
    if fill_col is None:
        fill_color = color
    else:
        fill_color = load_color(parts[fill_col], color)
    all_features.append((chr_name, start, end, strand, caption, color, fill_color))
handle.close()
telomere_length = 0.01 * max_length
print("%i features loaded" % len(all_features))
if per_page == 1:
    # Use portrait A4
    page_size = (21 * cm, 29.7 * cm)
    chr_percentage = 0.06
    label_percentage = 0.15
    label_size = 12
else:
    # There is a hard-coded half inch left and right margin.
    # These widths get 12 chr on a landscape A4 page nicely.
    page_size = (per_page * 2.2633333 * cm + inch, 21 * cm)
    chr_percentage = 0.15
    label_percentage = 0.135
    label_size = 5


def draw_page(selected_refs):
    chr_diagram = BasicChromosome.Organism()
    chr_diagram.page_size = page_size
    chr_diagram._legend_height = 0

    for name, length in selected_refs:
        features = []

        # Add the N-regions
        for n, start, end in n_regions:
            if n == name:
                # Want to use a border and fill color, needs Biopython 1.62
                features.append((start, end, None, "", colors.black, colors.lightgrey))
        for n, start, end, strand, caption, color, fill_color in all_features:
            if n == name:
                features.append((start, end, strand, caption, color, fill_color))

        cur_chromosome = BasicChromosome.Chromosome(name)
        cur_chromosome.scale_num = max_length + 2 * telomere_length
        cur_chromosome.chr_percent = chr_percentage
        cur_chromosome.label_sep_percent = label_percentage
        cur_chromosome.label_size = label_size
        cur_chromosome._color_labels = True

        # Add an opening spacer (to center all chromosomes vertically)
        space = BasicChromosome.SpacerSegment()
        space.scale = (cur_chromosome.scale_num - length) * 0.5 - telomere_length
        space.chr_percent = chr_percentage
        cur_chromosome.add(space)

        # Add an opening telomere
        start = BasicChromosome.TelomereSegment()
        start.scale = telomere_length
        start.chr_percent = chr_percentage
        cur_chromosome.add(start)

        # Add a body - using bp as the scale length here.
        body = BasicChromosome.AnnotatedChromosomeSegment(length, features, colors.blue)
        body.scale = length
        body.chr_percent = chr_percentage
        cur_chromosome.add(body)

        # Add a closing telomere
        end = BasicChromosome.TelomereSegment(inverted=True)
        end.scale = telomere_length
        end.chr_percent = chr_percentage
        cur_chromosome.add(end)

        # Add an closing spacer
        space = BasicChromosome.SpacerSegment()
        space.scale = (cur_chromosome.scale_num - length) * 0.5 - telomere_length
        space.chr_percent = chr_percentage
        cur_chromosome.add(space)

        # This chromosome is done
        chr_diagram.add(cur_chromosome)

        print("%s %i %i" % (name, length, len(features)))
    return chr_diagram


# Open PDF
c = canvas.Canvas(pdf_file, page_size)
c.setTitle(main_caption)
while refs:
    selected = refs[:per_page]
    refs = refs[per_page:]
    if len(selected) > 1:
        print("New page, %s to %s" % (selected[0][0], selected[-1][0]))
    else:
        print("New page, %s" % (selected[0][0]))
    with warnings.catch_warnings():
        # BiopythonWarning: Too many labels to avoid overlap
        warnings.simplefilter("ignore", BiopythonWarning)
        d = draw_page(selected).draw(None, main_caption)
    renderPDF.draw(d, c, 0, 0)
    c.showPage()
# Close PDF
c.save()
print("Done")
