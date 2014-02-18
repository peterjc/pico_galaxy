Galaxy wrapper for samtools idxstats
====================================

This wrapper is copyright 2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This is a wrapper for part of the command line samtools suite, v0.1.19

This wrapper is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/samtools_idxstats


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
samtools 0.1.19 if required.


Manual Installation
===================

This expects samtools to be on the $PATH, and was tested using v0.1.19.

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/samtools_idxstats`` folder:

* ``samtools_idxstats.xml`` (the Galaxy tool definition)
* ``samtools_idxstats.py`` (the Python wrapper script)
* ``README.rst`` (this file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool. Just add the line, perhaps under the NGS tools section::

  <tool file="samtools_idxstats/samtools_idxstats.xml" />

If you wish to run the unit tests, also add this to ``tools_conf.xml.sample``
and move/copy the ``test-data`` files under Galaxy's ``test-data`` folder. Then::

    $ ./run_functional_tests.sh -id samtools_idxstats

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial public release
======= ======================================================================


Developers
==========

Development is one this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/samtools_idxstats

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf samtools_idxstats.tar.gz tools/samtools_idxstats/README.rst tools/samtools_idxstats/samtools_idxstats.xml tools/samtools_idxstats/samtools_idxstats.py tools/samtools_idxstats/tool_dependencies.xml test-data/ex1.bam test-data/ex1.idxstats.tabular

Check this worked::

    $ tar -tzf samtools_idxstats.tar.gz
    tools/samtools_idxstats/README.rst
    tools/samtools_idxstats/samtools_idxstats.xml
    tools/samtools_idxstats/samtools_idxstats.py
    tools/samtools_idxstats/tool_dependencies.xml
    test-data/ex1.bam
    test-data/ex1.idxstats.tabular


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
