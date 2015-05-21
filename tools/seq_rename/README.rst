Galaxy tool to rename FASTA, QUAL, FASTQ or SFF sequences
=========================================================

This tool is copyright 2011-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool is a short Python script (using Biopython library functions) to rename
sequences from a FASTA, QUAL, FASTQ, or SFF file based on an ID mapping gives as
two columns of a tabular file. The output order follows that of the sequence file,
and if there are duplicates in the input sequence file, there will be duplicates
in the output sequence file.

This tool is available from the Galaxy Tool Shed,

* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_rename

See also the sister tools to filter or select sequence files according to IDs
from column(s) of tabular file:

* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_filter_by_id
* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_select_by_id


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* ``seq_rename.py`` (the Python script)
* ``seq_rename.xml`` (the Galaxy tool definition)

The suggested location is in a dedicated ``tools/seq_rename`` folder.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line::

    <tool file="seq_rename/seq_rename.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id seq_rename

You will also need to install Biopython 1.54 or later. That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version.
v0.0.2  - Record script version when run from Galaxy.
        - Add unit test.
        - Check for errors using Python script's return code.
v0.0.3  - Link to Tool Shed added to help text and this documentation.
v0.0.4  - Automated installation of Biopython dependency.
        - Use reStructuredText for this README file.
        - Adopt standard MIT License.
        - Updated citation information (Cock et al. 2013).
        - Development moved to GitHub, https://github.com/peterjc/pico_galaxy
        - Renamed folder and adopted README.rst naming.
v0.0.5  - Correct automated dependency definition.
v0.0.6  - Simplified XML to apply input format to output data.
        - Tool definition now embeds citation information.
        - If white space is found in the requested tabular field then only
          the first word is used as the identifier (with a warning to stderr).
v0.0.7  - Use the ``format_source=...`` tag.
        - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
======= ======================================================================


Developers
==========

This script and related tools were initially developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

Development has now moved to a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update --shed_target testtoolshed --check_diff ~/repositories/pico_galaxy/tools/seq_rename/
    ...

or::

    $ planemo shed_update --shed_target toolshed --check_diff ~/repositories/pico_galaxy/tools/seq_rename/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/seq_rename/
    ...
    $ tar -tzf shed_upload.tar.gz 
    test-data/four_human_proteins.fasta
    test-data/four_human_proteins.rename.fasta
    test-data/four_human_proteins.rename.tabular
    tools/seq_rename/README.rst
    tools/seq_rename/seq_rename.py
    tools/seq_rename/seq_rename.xml
    tools/seq_rename/tool_dependencies.xml


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
