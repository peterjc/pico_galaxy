Galaxy wrapper for drawing dotplots using MUMmer 3
==================================================

This wrapper is copyright 2014-2016 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This is a dotplot tool wrapping functionality from the MUMmer v3 suite.

S. Kurtz et al. (2004).
Versatile and open software for comparing large genomes.
Genome Biology (2004), 5:R12.
http://dx.doi.org/10.1186/gb-2004-5-2-r12

This wrapper is available to install into other Galaxy Instances via the Galaxy
Tool Shed at http://toolshed.g2.bx.psu.edu/view/peterjc/mummer

Automated Installation
======================

Automated installation from the Galaxy Tool Shed should be straightforward.
In addition to the wrapper itself, Galaxy should automatically download and
install the MUMmer files, gnuplot, and GhostScript which provides the
``ps2pdf`` binary.


Manual Installation
===================

This expects MUMmer binaries (at least ``mummer``, ``nucmer``, ``promer``, and
``mummerplot``) and the tools ``gnuplot`` and ``ps2pdf`` (from GhostScript) to
be on the system ``$PATH``.

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/mummer`` folder:

* ``mummer.xml`` (the Galaxy tool definition)
* ``mummer.py`` (the Python wrapper script)
* ``README.rst`` (this file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. Just add the line::

  <tool file="mummer/mummer.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id mummerplot_wrapper

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial public release.
v0.0.2  - Correct typo in Tool identifier.
        - Include MUMmer citation in Galaxy XML markup.
v0.0.3  - Install ``ps2pdf`` using Tool Shed's GhostScript package.
v0.0.4  - Install ``gnuplot`` using Tool Shed's gnuplot package.
        - Test case added.
v0.0.5  - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
v0.0.6  - PEP8 style changes to Python script (internal change only).
        - Fixed display of input parameter help text.
======= ======================================================================


Developers
==========

Development is on GitHub at:
https://github.com/peterjc/pico_galaxy/tree/master/tools/mummer


For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/mummer/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/mummer/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/mummer/
    ...
    $ tar -tzf shed_upload.tar.gz 
    test-data/magic.pdf
    test-data/magic.png
    test-data/rhodopsin_nucs.fasta
    test-data/three_human_mRNA.fasta
    tools/mummer/README.rst
    tools/mummer/mummer.py
    tools/mummer/mummer.xml
    tools/mummer/tool_dependencies.xml


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
MUMmer is available and licenced separately.
