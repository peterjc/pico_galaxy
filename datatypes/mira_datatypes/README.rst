Galaxy datatypes for MIRA Assembly Format (MAF)
===============================================

These Galaxy datatypes are copyright 2013 by Peter Cock, The James Hutton
Institute (formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.

See the licence text below (MIT licence).

This code defines a 'mira' datatype within Galaxy for the MIRA Assembly Format
(specifically v2 of the format introduced by MIRA 4.0 and the development
preview releases of MIRA v3.9). See http://chevreux.org/projects_mira.html
and https://sourceforge.net/projects/mira-assembler/ for background.

This format is not to be confused with the existing 'maf' datatype within Galaxy
for the unrelated Multiple (sequence) Alignment Format (MAF).

It is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/mira_datatypes


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - First public release
======= ======================================================================


Installation
============

Doing this automatically via the Galaxy Tool Shed is probably simplest, and will
happen automatically as a dependency of the Galaxy wrapper for MIRA v4.0.


Manual Installation
===================

Normally you would install this via the Galaxy ToolShed, which would move
the provided mira.py file into a suitable location and process the
datatypes_conf.xml entry to be combined with your local configuration.

However, if you really want to this should work for a manual install. Add
the following lines to the datatypes_conf.xml file in the Galaxy main folder::

    <datatype extension="mira" type="galaxy.datatypes.mira:MiraAssemblyFormat" mimetype="text/plain" display_in_upload="true"/>

and later in the sniffer section::

    <sniffer type="galaxy.datatypes.mira:MiraAssemblyFormat"/>

Also create the file lib/galaxy/datatypes/mira.py by moving, copying or linking
the mira.py file provided in this tar-ball.  Finally add 'import mira' near
the start of file lib/galaxy/datatypes/registry.py (after the other import
lines).


Bug Reports
===========

You can file an issue here https://github.com/peterjc/pico_galaxy/issues or ask
us on the Galaxy development list http://lists.bx.psu.edu/listinfo/galaxy-dev


Developers
==========

Development is done on this GitHub repository:
https://github.com/peterjc/pico_galaxy

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball I use
the following command from the mira_datatypes  folder::

    $ tar -czf mira_datatypes.tar.gz README.rst datatypes_conf.xml mira.py

Check this worked::

    $ tar -tzf mira_datatypes.tar.gz
    README.rst
    datatypes_conf.xml
    mira.py

For development, rather than having a local ToolShed running, I currently
use a symlink from lib/galaxy/datatypes/mira.py to the actual file as
described above.


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

NOTE: This is the licence for the Galaxy MIRA datatypes **only**. MIRA itself
is available and licenced separately.
