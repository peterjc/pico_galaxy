Galaxy wrapper for the MIRA assembly program (v4.9)
===================================================

This tool is copyright 2011-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (to collect the MIRA output and move it
to where Galaxy expects the files) and associated Galaxy wrapper XML file.

It is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/mira4_9

It uses a Galaxy datatype definition 'mira' for the MIRA Assembly Format,
http://toolshed.g2.bx.psu.edu/view/peterjc/mira_datatypes

A separate wrapper for MIRA v4.0 is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/mira4_assembler

A separate wrapper for MIRA v3.4 is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/mira_assembler

Automated Installation
======================

This should be straightforward. Via the Tool Shed, Galaxy should automatically
install the `mira` datatype, samtools, and download and install the precompiled
binary for MIRA v4.9.5 for the Galaxy wrapper, and offer to run any tests.

For MIRA 4, the Galaxy wrapper has been split in two, allowing separate
cluster settings for de novo usage (high RAM) and mapping (lower RAM).
Consult the Galaxy adminstration documentation for your cluster setup.

WARNING: For larger tasks, be aware that MIRA can require vast amounts
of RAM and run-times of over a week are possible. This tool wrapper makes
no attempt to spot and reject such large jobs.


Manual Installation
===================

First install the 'mira' datatype for Galaxy, available here:

* http://toolshed.g2.bx.psu.edu/view/peterjc/mira_datatypes 

There are various Python and XML files to install into Galaxy:

* ``mira_bait_4_9.xml`` (the Galaxy tool definition for mirabait)
* ``mira_bait_4_9.py`` (the Python wrapper script for mirabait)

The suggested location is a new ``tools/mira4_9`` folder. You will
also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool::

  <tool file="mira4_9/mira_bait_4_9.xml" />
  ...

You will also need to install MIRA 4.9, we used version 4.9.5, and define the
environment variable ``$MIRA4_9`` pointing at the folder containing the
binaries. See:

* http://chevreux.org/projects_mira.html
* http://sourceforge.net/projects/mira-assembler/

You may wish to use different cluster setups for the de novo and mapping
tools, see above.

You will also need to install samtools (for generating a BAM file from MIRA's
SAM output).

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id mira_bait_4_9
    ...


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version for MIRA 4.9.5, based on wrapper for v4.0.2
======= ======================================================================


Developers
==========

Development is on a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/mira_assembler_4_9

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/mira4_9/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/mira4_9/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/mira4_9/
    ...
    $ tar -tzf shed_upload.tar.gz 
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
