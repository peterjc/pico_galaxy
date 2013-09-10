Galaxy tool to draw a Venn Diagram with up to 3 sets
====================================================

This tool is copyright 2011 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool is a short Python script (using both the Galaxy and Biopython library
functions) to extract ID lists from tabular, FASTA, FASTQ or SFF files to build
sets, which are then drawn using the R limma package function vennDiagram
(called from Python using rpy).

There are just two files to install:

* venn_list.py (the Python script)
* venn_list.xml (the Galaxy tool definition)

The suggested location is in the Galaxy folder tools/plotting next to other
graph drawing tools.

You will also need to modify the tools_conf.xml file to tell Galaxy to offer the
tool. The suggested location is in the "Graph/Display Data" section. Simply add
the line::

  <tool file="plotting/venn_list.xml" />

You will also need to install Biopython 1.54 or later, and the R/Bioconductor
pacakge limma. You should already have rpy installed for other Galaxy tools.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.3  - Initial public release.
v0.0.4  - Ignore blank lines when loading IDs from tabular files
v0.0.5  - Explicit Galaxy error handling of return codes
v0.0.6  - Added unit tests.
        - Use reStructuredText for this README file.
======= ======================================================================


Developers
==========

This script and related tools are being developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf venn_list.tar.gz tools/plotting/venn_list.* test-data/venn_list1.pdf test-data/venn_list.tabular test-data/rhodopsin_proteins.fasta

Check this worked::

    $ tar -tzf venn_list.tar.gz
    tools/plotting/venn_list.py
    tools/plotting/venn_list.rst
    tools/plotting/venn_list.xml
    test-data/venn_list1.pdf
    test-data/venn_list.tabular
    test-data/rhodopsin_proteins.fasta


Licence (MIT/BSD style)
=======================

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
