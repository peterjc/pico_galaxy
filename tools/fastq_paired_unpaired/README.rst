Galaxy tool to divide FASTQ files into paired and unpaired reads
================================================================

This tool is copyright 2010-2023 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

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
the tool from the Galaxy Tool Shed, and Biopython, and allow you to run the unit
tests.


Manual Installation
===================

There are just two files to install:

* ``fastq_paired_unpaired.py`` (the Python script)
* ``fastq_paired_unpaired.xml`` (the Galaxy tool definition)

The suggested location is in the Galaxy folder tools/fastq next to other FASTQ
tools provided with Galaxy.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool. One suggested location is next to the ``fastq_filter.xml`` entry, or use
a dedicated folder like ``tools/fastq_paired_unpaired``. Then simply add the line::

    <tool file="fastq_paired_unpaired/fastq_paired_unpaired.xml" />

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
        - Updated citation information (Cock et al. 2013).
        - Development moved to GitHub.
v0.0.9  - Renamed folder and adopted README.rst naming.
        - Removed some unused code in the Python script.
v0.1.0  - Switch to using Biopython (easier to use script outside of Galaxy).
        - Leaves FASTQ plus lines blank (smaller output files).
        - Tool definition now embeds citation information.
v0.1.1  - Reorder XML elements (internal change only).
        - Use ``format_source=...`` tag.
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
v0.1.2  - Belatedly declare Biopython dependency via Tool Shed.
v0.1.3  - Minor internal changes to Python script for error reporting & style.
        - Updated to point at Biopython 1.67 (latest version in Tool Shed).
        - Explicit dependency on ``galaxy_sequence_utils``.
v0.1.4  - Use ``<command detect_errors="aggressive">`` (internal change only).
        - Single quote command line arguments (internal change only).
v0.1.5  - Bump Biopython dependency version for Python 3 fixes.
======= ======================================================================


Developers
==========

This script and other tools for filtering FASTA, FASTQ and SFF files were
initially developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/fasta_filter

Development has now moved to a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/fastq_paired_unpaired

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff tools/fastq_paired_unpaired/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff tools/fastq_paired_unpaired/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only tools/fastq_paired_unpaired/
    ...
    $ tar -tzf shed_upload.tar.gz
    test-data/sanger-pairs-forward.fastq
    test-data/sanger-pairs-interleaved.fastq
    test-data/sanger-pairs-mixed.fastq
    test-data/sanger-pairs-reverse.fastq
    test-data/sanger-pairs-singles.fastq
    tools/fastq_paired_unpaired/README.rst
    tools/fastq_paired_unpaired/fastq_paired_unpaired.py
    tools/fastq_paired_unpaired/fastq_paired_unpaired.xml


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
