Galaxy wrapper for the CLC Assembly Cell suite from CLCbio
==========================================================

This wrapper is copyright 2013-2017 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

CLC Assembly Cell is the commercial command line assembly suite from CLCbio.
It uses SIMD instructions to parallelize and accelerate their assembly
algorithms, and is also very memory efficient making it an appealing choice
for complex genomes where the RAM requirements exclude other popular tools.

For more information:
http://www.clcbio.com/products/clc-assembly-cell/

You can download the CLC Assembly Cell User Manual here, currently v4.2
http://www.clcbio.com/files/usermanuals/CLC_Assembly_Cell_User_Manual.pdf

There is also an online manual here:
http://clcsupport.com/clcassemblycell/current/index.php?manual=Introduction.html

There is currently a free trial download here:
http://www.clcbio.com/?action=transfer_user&productVersion=4.2&productID=6982&productName=CLC+Assembly+Cell&nonce=db842e3f95

This wrapper is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/clc_assembly_cell

This Galaxy wrapper was written and tested using CLC Assembly Cell v4.1.0,
specifically ``clc_assember`` and ``clc_mapper`` binaries for 64 bit Linux
which report version 4.10.86742 at the command line.


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and
install the wrapper from the Galaxy Tool Shed. However, you will need to
manually install the CLC Assembly Cell software, and setup the environment
variable ``$CLC_ASSEMBLY_CELL`` to the directory containing the binaries
(in particular, binaries ``clc_assembler``, ``clc_mapper`` and
``clc_cas_to_sam``). For example::

    $ export CLC_ASSEMBLY_CELL=/opt/clcbio/clc-assembly-cell-4.1.0-linux_64/

If your CLC Bio licence is restricted to specific machines on your cluster,
use Galaxy's job configuration settings to ensure CLC jobs are only sent
to those licenced computers. For SGE, we use the ``-l hostname="..."``
option to do this. Alternatively your cluster administrator might setup
a dedicated job queue.


Manual Installation
===================

First install the CLC Assembly Cell sortware as described above.

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/clcbio/`` folder:

* ``clc_assembler.xml`` (Galaxy tool definition)
* ``clc_mapper.xml`` (Galaxy tool definition)
* ``README.rst`` (this file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer
the tools. Just all these line, for example next to other assembly tools::

  <tool file="clc_assembly_cell/clc_assembler.xml" />
  <tool file="clc_assembly_cell/clc_mapper.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then run::

    $ ./run_tests.sh -id clc_assembler
    $ ./run_tests.sh -id clc_mapper

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial public release.
v0.0.2  - Actually use the ``$CLC_ASSEMBLY_CELL`` environment variable.
        - Enable and fixed the tests.
v0.0.3  - Reorder XML elements (internal change only).
        - Added citation tags.
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
v0.0.4  - Bug fix for ``<version_command>`` to capture tool version.
v0.0.5  - Support the ``-u`` or ``--discardunmapped`` option to discard
          unmapped reads in the CLC Mapper wrapper.
        - Bug fix to use the ``$CLC_ASSEMBLY_CELL`` environment variable when
          calling ``clc_cas_to_sam`` in the CLC MApper wrapper.
v0.0.6  - Update tool XML for optional vs required parameters.
        - Update URLs in documentation.
======= ======================================================================


Developers
==========

Development is on this itHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/clc_assembly_cell

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/clc_assembly_cell/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/clc_assembly_cell/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/clc_assembly_cell/
    ...
    $ tar -tzf shed_upload.tar.gz
    test-data/NC_010642.fna
    tools/clc_assembly_cell/README.rst
    tools/clc_assembly_cell/clc_assembler.xml
    tools/clc_assembly_cell/clc_mapper.xml
    tools/clc_assembly_cell/tool_dependencies.xml


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

NOTE: This is the licence for the Galaxy Wrapper only. The CLCbio tools are
commercial, and are available and licenced separately.
