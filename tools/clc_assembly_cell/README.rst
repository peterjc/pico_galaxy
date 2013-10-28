Galaxy wrapper for the CLC Assembly Cell suite from CLCBio
==========================================================

This wrapper is copyright 2013 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

CLC Assembly Cell is the commercial command line assembly suite from CLCBio.
It uses SIMD instructions to parallelize and accelerate their assembly
algorithms, and is also very memory efficient making it an appealing choice
for complex genomes where the RAM requirements exclude other popular tools.

For more information:
http://www.clcbio.com/products/clc-assembly-cell/

You can download the CLC Assembly Cell User Manual here, currently v4.2
http://www.clcbio.com/files/usermanuals/CLC_Assembly_Cell_User_Manual.pdf

There is currently a free trial download here:
http://www.clcbio.com/?action=transfer_user&productVersion=4.2&productID=6982&productName=CLC+Assembly+Cell&nonce=db842e3f95

This wrapper is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/clc_assembly_cell

This Galaxy wrapper was written and tested using CLC Assembly Cell
version 4.10.86742


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and
install the wrapper from the Galaxy Tool Shed. However, you will need to
manually install the CLC Assembly Cell software, and setup the environment
variable ``$CLC_ASSRMBLY_CELL`` to the directory containing the binaries
(and in particular, the ``clc_assembler`` binary). For example:

$ export CLC_ASSRMBLY_CELL=/opt/clcBio/clc-assembly-cell-4.1.0-linux_64/


Manual Installation
===================

First install the CLC Assembly Cell sortware as described above.

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a tools/clcbio folder:

* clc_assembler.xml (the Galaxy tool definition)
* README.rst (this file)

You will also need to modify the tools_conf.xml file to tell Galaxy to offer the
tools. Just all these line, for example next to other assembly tools::

  <tool file="clc_assembly_cell/clc_assembler.xml" />
  <tool file="clc_assembly_cell/clc_mapper.xml" />

If you wish to run the unit tests, also add this to tools_conf.xml.sample
and move/copy the test-data files under Galaxy's test-data folder. Then::

    $ ./run_functional_tests.sh -id clc_assembler

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

Development is on this itHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/clc_assembly_cell

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf clcbio.tar.gz tools/clc_assembly_cell/README.rst tools/clc_assembly_cell/clc_assembler.xml tools/clc_assembly_cell/clc_mapper.xml

Check this worked::

    $ tar -tzf clcbio.tar.gz
    tools/clc_assembly_cell/README.rst
    tools/clc_assembly_cell/clc_assembler.xml
    tools/clc_assembly_cell/clc_mapper.xml


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

NOTE: This is the licence for the Galaxy Wrapper only. The CLCBio tools are
commercial, and are available and licenced separately.
