Galaxy tool to wrap the MIRA sequence assembly program (v3.4)
=============================================================

This tool is copyright 2011-2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (to collect the MIRA output and move it
to where Galaxy expects the files, and convert MIRA's TCS file into a tab
separated file for use in Galaxy).

It is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/mira_assembler 


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and
install the precompiled binary for MIRA v3.4.0 for the Galaxy wrapper,
and run any tests.


Manual Installation
===================

There are just two Galaxy files to install:

* mira.py (the Python script)
* mira.xml (the Galaxy tool definition)

The suggested location is a new tools/mira_3_4 folder. You will also need to
modify the tools_conf.xml file to tell Galaxy to offer the tool, and also do
this to tools_conf.xml.sample in order to run any tests::

  <tool file="mira_3_4/mira.xml" />

You will also need to install MIRA, we used version 3.4.1.1. See:

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
v0.0.1  - Initial version (working prototype, using MIRA 3.2.1)
v0.0.2  - Improve capture of stdout/stderr (should see it as it runs)
v0.0.3  - Support Ion Torrent reads, now requires MIRA 3.4.0 or later
          (some other switches changed, e.g. -OUT rrol to rrot, which
          means the wrapper no longer works with MIRA 3.2.x)
        - The contig summary file (TCS file) was removed in MIRA 3.4
        - Report all missing output files (not just first missing one)
v0.0.4  - Fix problem with backbone arguments inroduced in v0.0.3
v0.0.5  - Implement the <version_command> tag to record the wrapper
          version and the MIRA version being used.
        - Check using MIRA 3.4 (later versions have a different API)
v0.0.6  - Tell MIRA to use /tmp for temporary files
        - Tell MIRA to ignore long read names (otherwise it aborts)
v0.0.7  - Automated installation of the 64 bit Linux MIRA binary.
v0.0.8  - Basic unit test added.
        - Link to Tool Shed added to help text and this documentation.
        - Use reStructuredText for this README file.
        - Adopted standard MIT licence.
        - Updated citation information (Cock et al. 2013).
======= ======================================================================


Developers
==========

This script and related tools are being developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf mira_wrapper.tar.gz tools/mira_3_4/mira.* tools/mira_3_4/tool_dependencies.xml test-data/tvc_mini.fastq test-data/tvc_contigs.fasta

Check this worked:

    $ tar -tzf mira_wrapper.tar.gz
    tools/mira_3_4/mira.py
    tools/mira_3_4/mira.rst
    tools/mira_3_4/mira.xml
    tools/mira_3_4/tool_dependencies.xml
    test-data/tvc_mini.fastq
    test-data/tvc_contigs.fasta


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
