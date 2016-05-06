Count sequence variants in region of interest in BAM file
=========================================================

This tool is copyright 2016 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This tool runs the command ``samtools view`` (taking advantage of an
indexed BAM file) to access only those reads mapped to the region of
interest (ROI), and then counts the different sequence variants found.

Internally this tool uses the command-line samtools suite.

This tool is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/count_roi_variants


Use outside of Galaxy
=====================

You just need the ``count_roi_variants.py`` script and to have samtools
on the ``$PATH``.  If you move/copy the script somewhere on your ``$PATH``
and then you can run it like this::

    $ count_roi_variants.py --help

Or, call the script at an explicit path::

    $ /path/to/my/stuff/count_roi_variants.py --help

Run like this it will use the current default Python. This was written and
tested under Python 2.7, but should also work under Python 2.6 and Python 3.
e.g.::

    $ python3 /path/to/my/stuff/count_roi_variants.py --help

The sample data and tests are designed to be run via Galaxy.


Automated Galaxy Installation
=============================

This should be straightforward, Galaxy should automatically download and install
samtools if required.


Manual Galaxy Installation
==========================

This expects samtools to be on the ``$PATH``, and was tested using v0.1.3

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/count_roi_variants`` folder:

* ``count_roi_variants.xml`` (the Galaxy tool definition)
* ``count_roi_variants.py`` (the Python wrapper script)
* ``README.rst`` (this file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool. Just add the line, perhaps under the NGS tools section::

  <tool file="count_roi_variants/count_roi_variants.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id count_roi_variants

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial public release
v0.0.2  - Cope with pipes in reference name (e.g. NCBI style FASTA naming)
v0.0.3  - Include a functional test for using an unrecognised reference.
======= ======================================================================


Developers
==========

Development is on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/count_roi_variants

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/count_roi_variants/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/count_roi_variants/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/count_roi_variants/
    ...
    $ tar -tzf shed_upload.tar.gz
    tools/count_roi_variants/README.rst
    tools/count_roi_variants/count_roi_variants.xml
    tools/count_roi_variants/count_roi_variants.py
    tools/count_roi_variants/tool_dependencies.xml
    ...


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
samtools is available and licenced separately.
