Galaxy tool to extract FASTQ paired read names
==============================================

This tool is copyright 2014-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script which handles FASTQ files with paired
reads to identify pairs, and can be used to recover missing partners.

The inputs are FASTQ files, and the output is two tabular files. The first
file has one line per read pair (showing the read pair's name with any
suffix like ``/1`` or ``/2`` removed). The second file has one line per
non-paired read (or paired read with unrecognised naming) giving the read
name.

This tool is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/fastq_pair_names


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
the tool from the Galaxy Tool Shed, and run the unit tests.


Manual Installation
===================

There are just two files to install:

* ``fastq_pair_names.py`` (the Python script)
* ``fastq_pair_names.xml`` (the Galaxy tool definition)

The suggested location is in the Galaxy folder tools/fastq next to other FASTQ
tools provided with Galaxy.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool. One suggested location is next to the ``fastq_filter.xml`` entry, or use
a dedicated folder like ``tools/fastq_pair_names``. Then simply add the line::

    <tool file="fastq_pair_names/fastq_pair_names.xml" />

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version.
v0.0.2  - Tool definition now embeds citation information.
v0.0.3  - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
v0.0.4  - Explicit dependency on ``galaxy_sequence_utils``.
        - Minor internal changes to Python script for error reporting & style.
v0.0.5  - Use ``<command detect_errors="aggressive">`` (internal change only).
        - Single quote command line arguments (internal change only).
        - Python 3 compatible print function.
======= ======================================================================


Developers
==========

Development is on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/fastq_pair_names

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff tools/fastq_pair_names/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff tools/fastq_pair_names/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only tools/fastq_pair_names/
    ...
    $ tar -tzf shed_upload.tar.gz 
    test-data/empty_file.dat
    test-data/sanger-pairs-mixed.fastq
    test-data/sanger-pairs-names.tabular
    tools/fastq_pair_names/README.rst
    tools/fastq_pair_names/fastq_pair_names.py
    tools/fastq_pair_names/fastq_pair_names.xml


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
