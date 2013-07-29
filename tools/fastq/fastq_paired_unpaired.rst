Galaxy tool to divide FASTQ files into paired and unpaired reads
================================================================

This tool is copyright 2010-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool is a short Python script which divides a FASTQ file into paired
reads, and single or orphan reads. You can have separate files for the
forward/reverse reads, or have them interleaved in a single file.

Note that the FASTQ variant is unimportant (Sanger, Solexa, Illumina, or even
Color Space should all work equally well).

This tool is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/fastq_paired_unpaired


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
the tool from the Galaxy Tool Shed, and run the unit tests


Manual Installation
===================

There are just two files to install:

* fastq_paired_unpaired.py (the Python script)
* fastq_paired_unpaired.xml (the Galaxy tool definition)

The suggested location is in the Galaxy folder tools/fastq next to other FASTQ
tools provided with Galaxy.

You will also need to modify the tools_conf.xml file to tell Galaxy to offer
the tool. One suggested location is next to the fastq_filter.xml entry. Simply
add the line::

    <tool file="fastq/fastq_paired_unpaired.xml" />

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version, using Biopython
v0.0.2  - Help text; cope with multiple pairs per template
v0.0.3  - Galaxy XML wrappers added
v0.0.4  - Use Galaxy library to handle FASTQ files (avoid Biopython dependency)
v0.0.5  - Handle Illumina 1.8 style pair names
v0.0.6  - Record script version when run from Galaxy
        - Added unit test (FASTQ file using Sanger naming)
v0.0.7  - Link to Tool Shed added to help text and this documentation.
v0.0.8  - Use reStructuredText for this README file.
        - Adopt standard MIT License.
======= ======================================================================


Developers
==========

This script and other tools for filtering FASTA, FASTQ and SFF files are
currently being developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/fasta_filter

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf fastq_paired_unpaired.tar.gz tools/fastq/fastq_paired_unpaired.* test-data/sanger-pairs-*.fastq

Check this worked::

    $ tar -tzf fastq_paired_unpaired.tar.gz
    tools/fastq/fastq_paired_unpaired.py
    tools/fastq/fastq_paired_unpaired.rst
    tools/fastq/fastq_paired_unpaired.xml
    test-data/sanger-pairs-forward.fastq
    test-data/sanger-pairs-interleaved.fastq
    test-data/sanger-pairs-mixed.fastq
    test-data/sanger-pairs-reverse.fastq
    test-data/sanger-pairs-singles.fastq


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
