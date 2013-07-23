Galaxy tool to select FASTA, QUAL, FASTQ or SFF sequences by ID
===============================================================

This tool is copyright 2011-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool is a short Python script (using Biopython library functions) to extract
sequences from a FASTA, QUAL, FASTQ, or SFF file based on the list of IDs given
by a column of a tabular file. The output order follows that of the tabular file,
and if there are duplicates in the tabular file, there will be duplicates in the
output sequence file.

This tool is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/seq_select_by_id

See also the sister tools to filter sequence files according to IDs from column(s)
of a tabular file (where the output order follows the sequence file, and any
duplicate IDs are ignored) and rename sequences:
http://toolshed.g2.bx.psu.edu/view/peterjc/seq_filter_by_id
http://toolshed.g2.bx.psu.edu/view/peterjc/seq_rename


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* seq_select_by_id.py (the Python script)
* seq_select_by_id.xml (the Galaxy tool definition)

The suggested location is in the Galaxy folder tools/filters next to the tool
for calling sff_extract.py for converting SFF to FASTQ or FASTA + QUAL.

You will also need to modify the tools_conf.xml file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line:

<tool file="filters/seq_select_by_id.xml" />

If you wish to run the unit tests, also add this to tools_conf.xml.sample
and move/copy the test-data files under Galaxy's test-data folder. Then:

$ ./run_functional_tests.sh -id seq_select_by_id

You will also need to install Biopython 1.54 or later. That's it.


History
=======

v0.0.1 - Initial version.
v0.0.3 - Ignore blank lines in input.
v0.0.4 - Record script version when run from Galaxy.
       - Basic unit test included.
v0.0.5 - Check for errors using Python script's return code.
v0.0.6 - Link to Tool Shed added to help text and this documentation.
v0.0.7 - Automatic installation of Biopython dependency.


Developers
==========

This script and related tools are being developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder:

$ tar -czf seq_select_by_id.tar.gz tools/filters/seq_select_by_id.* test-data/k12_ten_proteins.fasta test-data/k12_hypothetical.fasta test-data/k12_hypothetical.tabular

Check this worked:

$ tar -tzf seq_select_by_id.tar.gz
filter/seq_select_by_id.py
filter/seq_select_by_id.rst
filter/seq_select_by_id.xml
test-data/k12_ten_proteins.fasta
test-data/k12_hypothetical.fasta
test-data/k12_hypothetical.tabular


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
