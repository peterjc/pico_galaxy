Galaxy wrapper for pyani's average_nucleotide_identity.py
=========================================================

This wrapper is copyright 2017 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

pyani is a Python3 module that provides support for calculating average
nucleotide identity (ANI) and related measures for whole genome comparisons,
and rendering relevant graphical summary output. 


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
pyani and its dependencies (NCBI BLAST+ etc) using bioconda.


Manual Installation
===================
This wrapper expects ``average_nucleotide_identity.py`` script and underlying
binaries like ``blastn`` all to be on the ``$PATH``.

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/pyani`` folder:

* ``average_nucleotide_identity.xml`` (the Galaxy tool definition)
* ``README.rst`` (this README file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool by inserting a line like this (matching the chosen install path)::

  <tool file="pyani/average_nucleotide_identity.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id average_nucleotide_identity

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.0  - Started work
======= ======================================================================


Developers
==========

Development is on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/pyani

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff tools/pyani/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff tools/pyani/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only tools/pyani/
    ...
    $ tar -tzf shed_upload.tar.gz
    tools/pyani/README.rst
    tools/pyani/average_nucleotide_identity.xml


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
