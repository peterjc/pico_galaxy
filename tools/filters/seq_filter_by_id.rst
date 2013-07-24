Galaxy tool to filter FASTA, FASTQ or SFF sequences by ID
=========================================================

This tool is copyright 2010-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool is a short Python script (using both the Galaxy and Biopython library
functions) which divides a FASTA, FASTQ, or SFF file in two, those sequences with
or without an ID present in the specified column(s) of a tabular file. Example uses
include filtering based on search results from a tool like NCBI BLAST before
assembly.

This tool is available from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_filter_by_id

See also sister tools:

* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_select_by_id
* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_rename


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* seq_filter_by_id.py (the Python script)
* seq_filter_by_id.xml (the Galaxy tool definition)

The suggested location is in the Galaxy folder tools/filters next to the tool
for calling sff_extract.py for converting SFF to FASTQ or FASTA + QUAL.

You will also need to modify the tools_conf.xml file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line::

    <tool file="filters/seq_filter_by_id.xml" />

If you wish to run the unit tests, also add this to tools_conf.xml.sample
and move/copy the test-data files under Galaxy's test-data folder. Then::

    $ ./run_functional_tests.sh -id seq_filter_by_id

You will also need to install Biopython 1.54 or later. That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1   - Initial version, combining three separate scripts for each file format.
v0.0.4   - Record script version when run from Galaxy.
         - Faster FASTA code which preserves the original line wrapping.
         - Basic unit test included.
v0.0.5   - Check for errors using Python script's return code.
         - Cope with malformed FASTA entries without an identifier.
v0.0.6   - Link to Tool Shed added to help text and this documentation.
         - Automated installation of the Biopython dependency.
         - Use reStructuredText for this README file.
         - Adopt standard MIT License.
======= ======================================================================



Developers
==========

This script and related tools are being developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

This incorporates the previously used hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/fasta_filter

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf seq_filter_by_id.tar.gz tools/filters/seq_filter_by_id.* test-data/k12_ten_proteins.fasta test-data/k12_hypothetical.fasta test-data/k12_hypothetical.tabular

Check this worked::

    $ tar -tzf seq_filter_by_id.tar.gz
    filter/seq_filter_by_id.py
    filter/seq_filter_by_id.rst
    filter/seq_filter_by_id.xml
    test-data/k12_ten_proteins.fasta
    test-data/k12_hypothetical.fasta
    test-data/k12_hypothetical.tabular


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
