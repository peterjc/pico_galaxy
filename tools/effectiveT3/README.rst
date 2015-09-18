Galaxy wrapper for EffectiveT3 v1.0.1
=====================================

This wrapper is copyright 2011-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This is a wrapper for the command line Java tool EffectiveT3, v1.0.1,

Jehl, Arnold and Rattei.
Effective - a database of predicted secreted bacterial proteins
Nucleic Acids Research, 39(Database issue), D591-5, 2011.
http://dx.doi.org/10.1093/nar/gkq1154

Arnold, Brandmaier, Kleine, Tischler, Heinz, Behrens, Niinikoski, Mewes, Horn and Rattei.
Sequence-based prediction of type III secreted proteins.
PLoS Pathog. 5(4):e1000376, 2009.
http://dx.doi.org/10.1371/journal.ppat.1000376

http://effectors.org/

This wrapper is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/effectivet3


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
the Jar files for effectiveT3 v1.0.1 and the three models (animal, plant and std).


Manual Installation
===================

You can change the path by setting the environment variable ``$EFFECTIVET3`` to the
relevant folder, but by default it expects the following files to be installed
at these locations::

    /opt/EffectiveT3/TTSS_GUI-1.0.1.jar
    /opt/EffectiveT3/module/TTSS_ANIMAL-1.0.1.jar
    /opt/EffectiveT3/module/TTSS_PLANT-1.0.1.jar
    /opt/EffectiveT3/module/TTSS_STD-1.0.1.jar

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a tools/effectiveT3 folder:

* ``effectiveT3.xml`` (the Galaxy tool definition)
* ``effectiveT3.py`` (the Python wrapper script)
* ``README.rst`` (this file)

Also copy ``effectiveT3.loc.sample`` to ``effectiveT3.loc`` in the ``tool-data``
folder (and edit if appropriate, e.g. to add or remove a model).

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. If you are using other protein analysis tools like TMHMM or SignalP, put
it next to them. Just add the line::

  <tool file="effectiveT3/effectiveT3.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id effectiveT3

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.7  - Initial public release
v0.0.8  - Include effectiveT3.loc.sample in Tool Shed
v0.0.9  - Check the return code for errors in the XML
v0.0.10 - Added unit test
v0.0.11 - Automated installation
        - Record version of Python script when called from Galaxy
        - Link to Tool Shed added to help text and this documentation.
v0.0.12 - More explicit naming of the output dataset.
        - Adopt standard MIT licence.
        - Use reStructuredText for this README file.
        - Updated citation information (Cock et al. 2013).
        - Development moved to GitHub, https://github.com/peterjc/pico_galaxy
v0.0.13 - Relax unit test to allow for small floating point score difference.
        - Tool definition now embeds citation information.
v0.0.14 - Fixed error handling in ``effectiveT3.py``.
v0.0.15 - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
v0.0.16 - Updated URLs to download the tool and models.
        - Includes new standard classification model v2.0.1 (Sep 2015)
          as the default entry in ``tool-data/effectiveT3.loc``
======= ======================================================================


Developers
==========

This script and related tools were initially developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

Development has now moved to a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/effectiveT3

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/effectiveT3/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/effectiveT3/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/effectiveT3/
    ...
    $ tar -tzf shed_upload.tar.gz
    tool-data/effectiveT3.loc.sample
    test-data/empty.fasta
    test-data/empty_effectiveT3.tabular
    test-data/four_human_proteins.fasta
    test-data/four_human_proteins.effectiveT3.tabular
    tool-data/effectiveT3.loc.sample
    tools/effectiveT3/README.rst
    tools/effectiveT3/effectiveT3.py
    tools/effectiveT3/effectiveT3.xml
    tools/effectiveT3/tool_dependencies.xml


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

NOTE: This is the licence for the Galaxy Wrapper only.
EffectiveT3 is available and licenced separately.
