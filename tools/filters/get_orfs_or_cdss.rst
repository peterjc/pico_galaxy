Galaxy tool to find ORFs or simple CDSs
=======================================

This tool is copyright 2011-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool is a short Python script (using Biopython library functions)
to search nucleotide sequences for open reading frames (ORFs) or coding
sequences (CDSs) where the first potential start codon is used. See the
help text in the XML file for more information.

This tool is available from the Galaxy Tool Shed at:
* http://toolshed.g2.bx.psu.edu/view/peterjc/get_orfs_or_cdss


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* get_orfs_or_cdss.py (the Python script)
* get_orfs_or_cdss.xml (the Galaxy tool definition)

If you are installing this manually (rather than via the Tool Shed), the
suggested location is in the Galaxy folder tools/filters next to the tool
for calling sff_extract.py for converting SFF to FASTQ or FASTA + QUAL.
You will also need to modify the tools_conf.xml file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line::

    <tool file="filters/get_orfs_or_cdss.xml" />

You will also need to install Biopython 1.54 or later. If you want to run
the unit tests, include this line in tools_conf.xml.sample and the sample
FASTA files under the test-data directory. Then:

    ./run_functional_tests.sh -id get_orfs_or_cdss

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1   - Initial version.
v0.0.2   - Correct labelling issue on reverse strand.
         - Use the new <stdio> settings in the XML wrappers to catch errors
v0.0.3   - Include unit tests.
         - Record Python script version when run from Galaxy.
v0.0.4   - Link to Tool Shed added to help text and this documentation.
v0.0.5   - Automated intallation of the Biopython dependency.
         - Use reStructuredText for this README file.
======= ======================================================================


Developers
==========

This script and related tools are being developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf get_orfs_or_cdss.tar.gz tools/filters/get_orfs_or_cdss.* tools/filters/repository_dependencies.xml test-data/get_orf_input*.fasta test-data/Ssuis.fasta

Check this worked::

    $ tar -tzf get_orfs_or_cdss.tar.gz
    filter/get_orfs_or_cdss.py
    filter/get_orfs_or_cdss.rst
    filter/get_orfs_or_cdss.xml
    tools/filters/repository_dependencies.xml
    test-data/get_orf_input.fasta
    test-data/get_orf_input.Suis_ORF.nuc.fasta
    test-data/get_orf_input.Suis_ORF.prot.fasta
    test-data/get_orf_input.t11_nuc_out.fasta
    test-data/get_orf_input.t11_open_nuc_out.fasta
    test-data/get_orf_input.t11_open_prot_out.fasta
    test-data/get_orf_input.t11_prot_out.fasta
    test-data/get_orf_input.t1_nuc_out.fasta
    test-data/get_orf_input.t1_prot_out.fasta
    test-data/Ssuis.fasta


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
