Galaxy tool to primer clip (trim) FASTA, FASTQ or SFF reads
===========================================================

This tool is copyright 2011-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (using the Galaxy library functions and
Biopython). It is available from the Galaxy Tool Shed here:
http://toolshed.g2.bx.psu.edu/view/peterjc/seq_primer_clip


Automated Installation
======================

This should be straightforward, provided you have installed Biopython 1.54 or later.
Galaxy should automatically install the Python script, and run the unit tests.


Manual Installation
===================

There are just two files to install:

* seq_primer_clip.py (the Python script)
* seq_primer_clip.xml (the Galaxy tool definition)

The suggested location is a new tools/primers folder. You will also need to
modify the tools_conf.xml file to tell Galaxy to offer the tool::

  <tool file="primers/seq_primer_clip.xml" />

If you wish to run the unit tests, also add this to tools_conf.xml.sample
and move/copy the test-data files under Galaxy's test-data folder. Then::

    $ ./run_functional_tests.sh -id seq_primer_clip

You will also need to install Biopython 1.54 or later. That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version (not publicly released)
v0.0.2  - Sort primers by length (longest and therefore most specific first)
v0.0.3  - Consider missing bases at start/end of read as mismatches
v0.0.4  - Apply minimum length to sequences with no match too
v0.0.5  - Count clipped & non-matched short reads separately, length bug fixes
v0.0.6  - Added some functional tests
v0.0.7  - Added error check for bad filename arguments
v0.0.8  - Record version of Python script when run from Galaxy.
        - Check for errors using Python script's return code.
v0.0.9  - Moved test data to workaround Galaxy Tool Shed limititation.
v0.0.10 - Include links to Tool Shed in help text and this README file.
        - Use reStructuredText for this README file.
        - Adopted standard MIT licence.
======= ======================================================================


Developers
==========

This script and related tools are being developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

This incorporates the previously used hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/fasta_filter

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf seq_primer_clip.tar.gz tools/primers/seq_primer_clip.* test-data/dop_primers.fasta test-data/MID4_GLZRM4E04_rnd30*

Check this worked::

    $ tar -tzf seq_primer_clip.tar.gz
    tools/primers/seq_primer_clip.py
    tools/primers/seq_primer_clip.rst
    tools/primers/seq_primer_clip.xml
    test-data/dop_primers.fasta
    test-data/MID4_GLZRM4E04_rnd30.fasta
    test-data/MID4_GLZRM4E04_rnd30.fastqsanger
    test-data/MID4_GLZRM4E04_rnd30_fclip.fasta
    test-data/MID4_GLZRM4E04_rnd30_fclip.fastqsanger
    test-data/MID4_GLZRM4E04_rnd30_fclip.sff
    test-data/MID4_GLZRM4E04_rnd30_frclip.fasta
    test-data/MID4_GLZRM4E04_rnd30_frclip.fastqsanger
    test-data/MID4_GLZRM4E04_rnd30_frclip.sff
    test-data/MID4_GLZRM4E04_rnd30.sff


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
