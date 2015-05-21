Galaxy wrapper for samtools depth
=================================

This wrapper is copyright 2013-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This is a wrapper for part of the command line samtools suite.

This wrapper is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/samtools_idxstats


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
samtools 1.1 if required.


Manual Installation
===================

This expects samtools to be on the $PATH, and was tested using v1.1

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/samtools_depth`` folder:

* ``samtools_depth.xml`` (the Galaxy tool definition)
* ``README.rst`` (this file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool. Just add the line, perhaps under the NGS tools section::

  <tool file="samtools_depth/samtools_depth.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id samtools_depth

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial public release
v0.0.2  - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
======= ======================================================================


Developers
==========

Development is on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/samtools_depth


For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update --shed_target testtoolshed --check_diff ~/repositories/pico_galaxy/tools/samtools_depth/
    ...

or::

    $ planemo shed_update --shed_target toolshed --check_diff ~/repositories/pico_galaxy/tools/samtools_depth/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/samtools_depth/
    ...
    $ tar -tzf shed_upload.tar.gz 
    test-data/ex1.bam
    test-data/ex1.depth.tabular
    test-data/sam_spec_padded.sam
    test-data/sam_spec_padded.bam
    test-data/sam_spec_padded.depth.tabular
    tools/samtools_depth/README.rst
    tools/samtools_depth/samtools_depth.xml
    tools/samtools_depth/tool_dependencies.xml


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
