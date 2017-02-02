Obsolete
========

This tool is now obsolete, having been replaced	by a more general version
covering the FASTA, FASTQ and SFF sequence formats in a single tool. You
should only install this tool if you need to support existing workflows
which used it.


Galaxy tool to filter FASTA sequences by ID
===========================================

This tool is copyright 2010-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (using the Galaxy library functions) which
divides a FASTA file in two, those sequences with or without an ID present in
the specified column(s) of a tabular file. Example uses include filtering based
on search results from a tool like NCBI BLAST, TMHMM or SignalP.

There are just two files to install:

* fasta_filter_by_id.py (the Python script)
* fasta_filter_by_id.xml (the Galaxy tool definition)

The suggested location is next to the similarly named fasta_filter_by_length.py
and fasta_filter_by_length.xml files which are included with Galaxy, i.e.
in the Galaxy folder tools/fasta_tools

You will also need to modify the tools_conf.xml file to tell Galaxy to offer
the tool. The suggested location is next to the fasta_filter_by_length.xml
entry. Simply add the line:

<tool file="fasta_tools/fasta_filter_by_id.xml" />

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version (not publicly released)
v0.0.2  - Allow both, just pos or just neg output files
v0.0.3  - Include FASTA in tool name
v0.0.4  - Deprecated, marked as hidden in the XML
v0.0.5  - Explicit dependency on ``galaxy_sequence_utils``.
        - Citation information (Cock et al. 2013).
        - Explicitly record version via ``<version_command>``.
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

    $ planemo shed_update -t testtoolshed --check_diff tools/fasta_filter_by_id/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff tools/fasta_filter_by_id/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only tools/fasta_filter_by_id/
    ...
    $ tar -tzf shed_upload.tar.gz
    tools/fasta_filter_by_id/README.rst
    tools/fasta_filter_by_id/fasta_filter_by_id.py
    tools/fasta_filter_by_id/fasta_filter_by_id.xml
    tools/fasta_filter_by_id/tool_dependencies.xml
    test-data/four_human_proteins.fasta
    test-data/blastp_four_human_vs_rhodopsin.tabular
    test-data/four_human_proteins_filter_a.fasta
    test-data/four_human_proteins_filter_b.fasta


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
