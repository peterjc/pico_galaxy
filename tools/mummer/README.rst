Galaxy wrapper for EffectiveT3 v1.0.1
=====================================

This wrapper is copyright 2014 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

This is a dotplot tool wrapping functionality from the MUMmer v3 suite.

S. Kurtz et al. (2004).
Versatile and open software for comparing large genomes.
Genome Biology (2004), 5:R12.
http://dx.doi.org/10.1186/gb-2004-5-2-r12

This wrapper is available to install into other Galaxy Instances via the Galaxy
Tool Shed at http://toolshed.g2.bx.psu.edu/view/peterjc/mummer

Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
the MUMmer files.

It also needs the tools ``gnuplot`` and ``ps2pdf`` to be installed and on the
system ``$PATH``.


Manual Installation
===================

This expects MUMmer binaries (at least ``mummer``, ``nucmer``, ``promer``, and
``mummerplot``) and the tools ``gnuplot`` and ``ps2pdf`` to be on the system
``$PATH``.

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/mummer`` folder:

* ``mummer.xml`` (the Galaxy tool definition)
* ``mummer.py`` (the Python wrapper script)
* ``README.rst`` (this file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. Just add the line::

  <tool file="mummer/mummer.xml" />

If you wish to run the unit tests, also add this to ``tools_conf.xml.sample``
and move/copy the ``test-data`` files under Galaxy's ``test-data`` folder. Then::

    $ ./run_tests.sh -id mumerplot_wrapper

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial public release.
v0.0.2  - Include MUMmer citation in Galaxy XML markup.
======= ======================================================================


Developers
==========

Development is on GitHub at:
https://github.com/peterjc/pico_galaxy/tree/master/tools/mummer

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf mummer.tar.gz tools/mummer/README.rst tools/mummer/mummer.xml tools/mummer/mummer.py tools/mummer/tool_dependencies.xml

Check this worked::

    $ tar -tzf mummer.tar.gz
    tools/mummer/README.rst
    tools/mummer/mummer.xml
    tools/mummer/mummer.py
    tools/mummer/tool_dependencies.xml


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
MUMmer is available and licenced separately.
