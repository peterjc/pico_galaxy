Galaxy tool to select FASTA, QUAL, FASTQ or SFF sequences by ID
===============================================================

This tool is copyright 2011-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool is a short Python script (using Biopython library functions) to extract
sequences from a FASTA, QUAL, FASTQ, or SFF file based on the list of IDs given
by a column of a tabular file. The output order follows that of the tabular file,
and if there are duplicates in the tabular file, there will be duplicates in the
output sequence file.

This tool is available from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_select_by_id

See also the sister tools to filter sequence files according to IDs from column(s)
of a tabular file (where the output order follows the sequence file, and any
duplicate IDs are ignored) and rename sequences:

* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_filter_by_id
* http://toolshed.g2.bx.psu.edu/view/peterjc/seq_rename


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* ``seq_select_by_id.py`` (the Python script)
* ``seq_select_by_id.xml`` (the Galaxy tool definition)

The suggested location is a dedicated ``tools/seq_select_by_id`` folder.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line::

    <tool file="seq_select_by_id/seq_select_by_id.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id seq_select_by_id

You will also need to install Biopython 1.54 or later. That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version.
v0.0.3  - Ignore blank lines in input.
v0.0.4  - Record script version when run from Galaxy.
        - Basic unit test included.
v0.0.5  - Check for errors using Python script's return code.
v0.0.6  - Link to Tool Shed added to help text and this documentation.
        - Automatic installation of Biopython dependency.
        - Use reStructuredText for this README file.
        - Adopt standard MIT License.
v0.0.7  - Updated citation information (Cock et al. 2013).
        - Fixed Biopython dependency setup.
        - Development moved to GitHub, https://github.com/peterjc/pico_galaxy
        - Renamed folder and adopted README.rst naming.
v0.0.8  - Corrected automated dependency definition.
v0.0.9  - Simplified XML to apply input format to output data.
        - Tool definition now embeds citation information.
        - Include input dataset name in output dataset names.
        - If white space is found in the requested tabular field then only
          the first word is used as the identifier (with a warning to stderr).
v0.0.10 - Includes testing of stdout messages.
        - Includes testing of failure modes.
v0.0.11 - Use the ``format_source=...`` tag.
        - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
        - Quote filenames in case of spaces (internal change only).
v0.0.12 - Python style changes (internal change only).
        - Use ``<command detect_errors="aggressive">`` (internal change only).
        - Depends on Biopython 1.67 via legacy Tool Shed package or bioconda.
        - Python 3 compatible print function.
v0.0.13 - Python 3 compatible exception handling.
v0.0.14 - Script works on Python 2 and 3 (fixed output file mode).
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

    $ planemo shed_update -t testtoolshed --check_diff tools/seq_select_by_id/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff tools/seq_select_by_id/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only tools/seq_select_by_id/
    ...
    $ tar -tzf shed_upload.tar.gz
    test-data/k12_hypothetical.fasta
    test-data/k12_hypothetical.tabular
    test-data/k12_hypothetical_alt.tabular
    test-data/k12_ten_proteins.fasta
    tools/seq_select_by_id/README.rst
    tools/seq_select_by_id/seq_select_by_id.py
    tools/seq_select_by_id/seq_select_by_id.xml
    tools/seq_select_by_id/tool_dependencies.xml


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
