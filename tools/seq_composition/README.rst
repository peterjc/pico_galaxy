Galaxy tool reporting sequence composition
==========================================

This tool is copyright 2014-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (using Biopython library functions) to
loop over given sequence files (in a range of formats including FASTA, FASTQ,
and SFF), and report the count of each letter (i.e. amino acids or bases).

This can be useful for sanity checking assemblies (e.g. proportion of N
bases) or looking at differences in base composition.

This tool is available from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_composition


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* ``seq_composition.py`` (the Python script)
* ``seq_composition.xml`` (the Galaxy tool definition)

The suggested location is in a dedicated ``tools/seq_composition`` folder.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line::

    <tool file="seq_composition/seq_composition.xml" />

You will also need to install Biopython 1.62 or later.

If you wish to run the unit tests, also	move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    ./run_tests.sh -id seq_composition

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version.
        - Tool definition now embeds citation information.
v0.0.2  - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
======= ======================================================================


Developers
==========

This script and related tools are being developed on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/seq_composition


For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/seq_composition/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/seq_composition/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/seq_composition/
    ...
    $ tar -tzf shed_upload.tar.gz 
    test-data/MID4_GLZRM4E04_rnd30_frclip.sff
    test-data/MID4_GLZRM4E04_rnd30_frclip.seq_composition.tabular
    test-data/ecoli.fastq
    test-data/ecoli.seq_composition.tabular
    test-data/four_human_proteins.fasta
    test-data/four_human_proteins.seq_composition.tabular
    tools/seq_composition/README.rst
    tools/seq_composition/seq_composition.py
    tools/seq_composition/seq_composition.xml
    tools/seq_composition/tool_dependencies.xml


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
