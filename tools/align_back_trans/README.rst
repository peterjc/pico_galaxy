Galaxy tool to back-translate a protein alignment to nucleotides
================================================================

This tool is copyright 2012-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (using Biopython library functions) to
load a protein alignment, and matching nucleotide FASTA file of unaligned
sequences, which are threaded onto the protein alignment in order to produce
a codon aware nucleotide alignment - which can be viewed as a back translation.

This tool is available from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/align_back_trans

The underlying Python script can also be used outside of Galaxy, for
details run::

    $ python align_back_trans.py

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

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool. One suggested location is in the multiple alignments section. Simply
add the line::

    <tool file="align_back_trans/align_back_trans.xml" />

You will also need to install Biopython 1.62 or later.

If you wish to run the unit tests, also	move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    ./run_tests.sh -id align_back_trans

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version, based on a previously written Python script
v0.0.2  - Optionally check the translation is consistent
v0.0.3  - First official release
v0.0.4  - Simplified XML to apply input format to output data.
        - Fixed error message when sequence length not a multiple of three.
v0.0.5  - More explicit error messages when seqences lengths do not match.
        - Tool definition now embeds citation information.
v0.0.6  - Reorder XML elements (internal change only).
        - Use ``format_source=...`` tag.
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
v0.0.7  - Minor Python code style improvements (internal change only).
======= ======================================================================


Developers
==========

This script was initially developed on this repository:
https://github.com/peterjc/picobio/blob/master/align/align_back_trans.py

With the addition of a Galaxy wrapper, developement moved here:
https://github.com/peterjc/pico_galaxy/tree/master/tools/align_back_trans

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/align_back_trans/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/align_back_trans/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/align_back_trans/
    ...
    $ tar -tzf shed_upload.tar.gz 
    test-data/demo_nucs.fasta
    test-data/demo_nucs_trailing_stop.fasta
    test-data/demo_prot_align.fasta
    test-data/demo_nuc_align.fasta
    tools/align_back_trans/README.rst
    tools/align_back_trans/align_back_trans.py
    tools/align_back_trans/align_back_trans.xml
    tools/align_back_trans/tool_dependencies.xml


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
