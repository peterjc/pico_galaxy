Galaxy tool to wrap the MIRA sequence assembly program (v3.4)
=============================================================

This tool is copyright 2011-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (to collect the MIRA output and move it
to where Galaxy expects the files) and associated Galaxy wrapper XML file.

It is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/mira4_assembler 


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and
install the precompiled binary for MIRA v4.0 for the Galaxy wrapper, and
run any tests.


Manual Installation
===================

There are just two Galaxy files to install:

* mira4.py (the Python script)
* mira4.xml (the Galaxy tool definition)

The suggested location is a new tools/mira_4_0 folder. You will also need to
modify the tools_conf.xml file to tell Galaxy to offer the tool, and also do
this to tools_conf.xml.sample in order to run any tests::

  <tool file="mira_4_0/mira4.xml" />

You will also need to install MIRA, we used version 4.0 RC2. See:

* http://chevreux.org/projects_mira.html
* http://sourceforge.net/projects/mira-assembler/

WARNING: This tool was developed to construct viral genome assembly and
mapping pipelines, for which the run time and memory requirements are
negligible. For larger tasks, be aware that MIRA can require vast amounts
of RAM and run-times of over a week are possible. This tool wrapper makes
no attempt to spot and reject such large jobs.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version (prototype using MIRA 4.0 RC2, and wrapper for v3.4)
======= ======================================================================


Developers
==========

Development is on a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/mira_4_0

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf mira4_wrapper.tar.gz tools/mira_4_0/README.rst tools/mira_4_0/mira4.xml tools/mira_4_0/mira4.py tools/mira_4_0/tool_dependencies.xml test-data/tvc_mini.fastq test-data/tvc_contigs_mira4.fasta

Check this worked::

    $ tar -tzf mira_wrapper.tar.gz
    tools/mira_4_0/README.rst
    tools/mira_4_0/mira4.xml
    tools/mira_4_0/mira4.py
    tools/mira_4_0/tool_dependencies.xml
    test-data/tvc_mini.fastq
    test-data/tvc_contigs_mira4.fasta


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
