Galaxy Tool for producing Biopython Chromosome Diagrams
=======================================================

Copyright 2011-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).


Introduction
============

This is a work in progress.

It has not yet been released to the Galaxy Tool Shed, but will eventually be
listed under http://toolshed.g2.bx.psu.edu/view/peterjc with my other Galaxy
tools and wrappers. Further work on font settings etc is probably needed
for broader use.

This tool uses the Biopython and ReportLab libraries to draw the diagram.
If you use this tool in scientific work leading to a publication, please
cite the Biopython application note:

Cock et al (2009) Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

The original Python script (and associated enhancements in Biopython 1.59)
were developed in order to produce Figure 2 of the potato chromosomes and the
detailed supplementary figures in this manuscript:

Jupe, F, Pritchard, L, Etherington, GJ, MacKenzie, K, Cock, PJA, et al. (2012)
Identification and localisation of the NB-LRR gene family within the potato
genome. BMC Genomics 13:75.
https://doi.org/10.1186/1471-2164-13-75

The current general Python script and its preliminary Galaxy wrapper, and
more associated enhancments in Biopython, were developed while producing
Figures 3 (potato chromosomes) and 4 (tomato chromosomes) in:

Jupe, F. et al. (2013) Resistance gene enrichment sequencing (RenSeq) enables
re-annotation of the NB-LRR gene family from sequenced plant genomes and
rapid mapping of resistance loci in segregating populations.
The Plant Journal 76(3):530-544
https://doi.org/10.1111/tpj.12307


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* ``chromosome_diagram.py`` (the Python script)
* ``chromosome_diagram.xml`` (the Galaxy tool definition)

The suggested location is a dedicated ``tools/chromosome_diagram`` folder.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. One suggested location is in the plotting section. Simply add the line::

    <tool file="chromosome_diagram/chromosome_diagram.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id chromosome_diagram

You will also need to install Biopython 1.62 or later, and ReportLab. That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version, not a public release to Tool Shed
v0.0.2  - Reorder XML elements (internal change only).
        - Depends on Biopython 1.67 via legacy Tool Shed package or bioconda.
v0.0.3  - Use ``<command detect_errors="aggressive">`` (internal change only).
        - Single quote command line arguments (internal change only).
        - Python 3 compatible print function.
======= ======================================================================


Developers
==========
Development is on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/chromosome_diagram

The script ``fetch_data.py`` should download *A. thaliana* data in order
to draw a diagram showing the position of the tRNA genes.


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
