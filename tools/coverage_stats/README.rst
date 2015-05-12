BAM coverage statistics using samtools idxstats and depth
=========================================================

This tool is copyright 2014-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

Internally this tool uses the command-line samtools suite.

This tool is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/coverage_stats


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
samtools 0.1.19 if required.


Manual Installation
===================

This expects samtools to be on the ``$PATH``, and was tested using v0.1.19.

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/coverage_stats`` folder:

* ``coverage_stats.xml`` (the Galaxy tool definition)
* ``coverage_stats.py`` (the Python wrapper script)
* ``README.rst`` (this file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tool. Just add the line, perhaps under the NGS tools section::

  <tool file="coverage_stats/coverage_stats.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id coverage_stats

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial public release
v0.0.2  - Cope with samtools' default depth limit using modified samtools,
          see https://github.com/samtools/samtools/pull/322
v0.0.3  - Cope with no coverage in final contigs.
v0.0.4  - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
======= ======================================================================


Developers
==========

Development is on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/coverage_stats

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_upload --shed_target testtoolshed --check_diff ~/repositories/pico_galaxy/tools/coverage_stats/
    ...

or::

    $ planemo shed_upload --shed_target toolshed --check_diff ~/repositories/pico_galaxy/tools/coverage_stats/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/coverage_stats/
    ...
    $ tar -tzf shed_upload.tar.gz
    tools/coverage_stats/README.rst
    tools/coverage_stats/coverage_stats.xml
    tools/coverage_stats/coverage_stats.py
    tools/coverage_stats/tool_dependencies.xml
    test-data/coverage_test.bam
    test-data/coverage_test.coverage_stats.tabular
    test-data/ex1.bam
    test-data/ex1.coverage_stats.tabular
    tools/coverage_stats/README.rst
    tools/coverage_stats/coverage_stats.xml
    tools/coverage_stats/coverage_stats.py
    tools/coverage_stats/tool_dependencies.xml


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
