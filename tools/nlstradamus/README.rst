Galaxy wrapper for NLStradamus v1.7 or v1.8 (C++ version)
=========================================================

This wrapper is copyright 2011-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

NLStradamus is a command line tool for predicting nuclear localization
signals (NLSs) in a FASTA file of proteins using a Hidden Markov Model (HMM).

This wrapper is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/nlstradamus

A. N. Nguyen Ba, A. Pogoutse, N. Provart, A. M. Moses.
NLStradamus: a simple Hidden Markov Model for nuclear localization signal prediction.
BMC Bioinformatics. 2009 Jun 29;10(1):202.
http://dx.doi.org/10.1186/1471-2105-10-202

http://www.moseslab.csb.utoronto.ca/NLStradamus

Early versions of NLStradamus did not have a native tabular output format, this
was added in version 1.7. Additionally a fast C++ implementation was added at
this point (early versions of NLStradamus came as a perl script only).

Version 1.8 fixed a C++ compilation issue on modern compilers, but is otherwise
unchanged.


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
the C++ implementation of NLStradamus v1.8, and run the unit tests.


Manual Installation
===================
This wrapper expects the compiled C++ binary "NLStradamus" to be on the ``$PATH``.

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a tools/protein_analysis folder:

* ``nlstradamus.xml`` (the Galaxy tool definition)
* ``nlstradamus.txt`` (this README file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool. If you are using other protein analysis tools like TMHMM or SignalP,
put it next to them. Just add the line (matching the chosen install path)::

  <tool file="protein_analysis/nlstradamus.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id nlstradamus

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.3  - Initial public release
v0.0.4  - Adding DOI link to reference
          (Documentation change only)
v0.0.5  - Assume non-zero return codes are errors
v0.0.6  - Show output help text using a table
        - Added unit tests
v0.0.7  - Automatic installation of the NLStradamus binary when installed
          via the Galaxy Tool Shed
v0.0.8  - Link to Tool Shed added to help text and this documentation.
        - Use reStructuredText for this README file.
        - Adopted standard MIT licence.
        - Updated citation information (Cock et al. 2013).
        - Development moved to GitHub, https://github.com/peterjc/pico_galaxy
v0.0.9  - Tool definition now embeds citation information.
v0.0.10 - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
======= ======================================================================


Developers
==========

This script and related tools were initially developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

Development has now moved to a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/nlstradamus


For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/nlstradamus/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/nlstradamus/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/nlstradamus/
    ...
    $ tar -tzf shed_upload.tar.gz
    test-data/empty.fasta
    test-data/empty_nlstradamus.tabular
    test-data/four_human_proteins.fasta
    test-data/four_human_proteins.nlstradamus.tabular
    tools/nlstradamus/README.rst
    tools/nlstradamus/nlstradamus.xml
    tools/nlstradamus/tool_dependencies.xml
    test-data/four_human_proteins.fasta


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
