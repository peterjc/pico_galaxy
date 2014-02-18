Galaxy tool to back-translate a protein alignment to nucleotides
================================================================

This tool is copyright 2012-2014 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (using Biopython library functions) to
load a protein alignment, and matching nucleotide FASTA file of unaligned
sequences, which are threaded onto the protein alignment in order to produce
a codon aware nucleotide alignment - which can be viewed as a back translation.

This tool is available from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/align_back_trans


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* ``align_back_trans.py`` (the Python script)
* ``align_back_trans.xml`` (the Galaxy tool definition)

The suggested location is in a dedicated ``tools/align_back_trans`` folder.

You will also need to modify the tools_conf.xml file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line::

    <tool file="align_back_trans/align_back_trans.xml" />

You will also need to install Biopython 1.62 or later. If you want to run
the unit tests, include this line in ``tools_conf.xml.sample`` and the sample
FASTA files under the test-data directory. Then::

    ./run_functional_tests.sh -id align_back_trans

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version, based on a previously written Python script
======= ======================================================================


Developers
==========

This script was initially developed on this repository:
https://github.com/peterjc/picobio/blob/master/align/align_back_trans.py

With the addition of a Galaxy wrapper, developement moved here:
https://github.com/peterjc/pico_galaxy/tree/master/tools/align_back_trans

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf align_back_trans.tar.gz tools/align_back_trans/README.rst tools/align_back_trans/align_back_trans.py tools/align_back_trans/align_back_trans.xml tools/align_back_trans/tool_dependencies.xml test-data/demo_nucs.fasta test-data/demo_prot_align.fasta test-data/demo_nuc_align.fasta

Check this worked::

    $ tar -tzf align_back_trans.tar.gz
    tools/align_back_trans/README.rst
    tools/align_back_trans/align_back_trans.py
    tools/align_back_trans/align_back_trans.xml
    tools/align_back_trans/tool_dependencies.xml
    test-data/demo_nucs.fasta
    test-data/demo_prot_align.fasta
    test-data/demo_nuc_align.fasta


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
