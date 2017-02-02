Obsolete
========

This tool is now obsolete, having been replaced	by a more general version
covering the FASTA, FASTQ and SFF sequence formats in a single tool. You
should only install this tool if you need to support existing workflows
which used it.

Galaxy tool to filter FASTQ sequences by ID
===========================================

This tool is copyright 2010-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (using the Galaxy library functions) which
divides a FASTQ file in two, those sequences with or without an ID present in
the specified column(s) of a tabular file. Example uses include filtering based
on search results from a tool like NCBI BLAST before assembly.

There are just two files to install:

* fastq_filter_by_id.py (the Python script)
* fastq_filter_by_id.xml (the Galaxy tool definition)

The suggested location is next to the similarly named fastq_filter.py and
fastq_filter.xml files which are included with Galaxy, i.e. in the Galaxy
folder tools/fastq

You will also need to modify the tools_conf.xml file to tell Galaxy to offer
the tool. The suggested location is next to the fastq_filter.xml entry. Simply
add the line:

<tool file="fastq/fastq_filter_by_id.xml" />

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial verion (not publicly released)
v0.0.2  - Allow both, just pos or just neg output files
        - Preserve the FASTQ variant in the XML wrapper
v0.0.3  - Fixed bug when generating non-matching FASTQ file only
v0.0.4  - Deprecated, marked as hidden in the XML
v0.0.5  - Explicit dependency on ``galaxy_sequence_utils``.
        - Citation information (Cock et al. 2013).
======= ======================================================================


Developers
==========

This script and other tools for filtering FASTA, FASTQ and SFF files were
initially developed on the following hg branches:
http://bitbucket.org/peterjc/galaxy-central/src/tools
http://bitbucket.org/peterjc/galaxy-central/src/fasta_filter

It is now under GitHub https://github.com/peterjc/pico_galaxy/

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff tools/fastq_filter_by_id/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff tools/fastq_filter_by_id/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only tools/fastq_filter_by_id/
    ...
    $ tar -tzf shed_upload.tar.gz
    README.rst
    fastq_filter_by_id.py
    fastq_filter_by_id.xml
    tool_dependencies.xml


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
