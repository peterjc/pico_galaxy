Galaxy wrapper for Command line NoD predictor (v1.3)
====================================================

This wrapper is copyright 2011-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

Command line NoD predictor is a tool for predicting nucleolar localization
sequences (NoLSs) in a FASTA file of proteins using a neural network. There
is also a webtool version at http://www.compbio.dundee.ac.uk/www-nod/

This NoD wrapper is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/peterjc/clinod


Automated Installation
======================

This should be straightforward, Galaxy should automatically download and install
the JAR file for clinod v1.3 and the batchman executable from the Stuttgart
Neural Network Simulator (SNNS), and set the ``$CLINOD`` to their folder.


Manual Installation
===================
This wrapper expects the java binary to be on the system PATH, and to be able
to access command line NoD as ``$CLINOD/clinod-1.3.jar`` which means if you
used ``/opt/clinod/clinod-1.3.jar`` set the environment variable ``$CLINOD``
to ``/opt/clinod``

Internally NoD calls the binary batchman v1.0 from the Stuttgart Neural Network
Simulator (SNNS) v4.2 or 4.3 software suite. This binary can either be on the
system path or located next to the JAR file, i.e. ``/opt/clinod/batchman``

To install the wrapper copy or move the following files under the Galaxy tools
folder, e.g. in a ``tools/clinod`` folder:

* ``clinod.xml`` (the Galaxy tool definition)
* ``README.rst`` (this file)

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. If you are using other protein analysis tools like TMHMM or SignalP, put
it next to them. Just add the line::

  <tool file="clinod/clinod.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    $ ./run_tests.sh -id clinod

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial public release
v0.0.2  - Treat non-zero return codes as errors
v0.0.3  - Describe output table in help
v0.0.4  - Added unit test
v0.0.5  - Link to Tool Shed added to help text and this documentation.
        - Automated tool installation.
v0.0.6  - Adopted standard MIT licence.
        - Use reStructuredText for this README file.
        - Updated citation information (Cock et al. 2013).
        - Development moved to GitHub, https://github.com/peterjc/pico_galaxy
v0.0.7  - Use the new ``$GALAXY_SLOTS`` environment variable for thread count.
        - Tool definition now embeds citation information.
v0.0.8  - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
======= ======================================================================


Developers
==========

This script and related tools were initially developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

Development has now moved to a dedicated GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/clinod

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update --shed_target testtoolshed --check_diff ~/repositories/pico_galaxy/tools/clinod/
    ...

or::

    $ planemo shed_update --shed_target toolshed --check_diff ~/repositories/pico_galaxy/tools/clinod/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/clinod/
    ...
    $ tar -tzf shed_upload.tar.gz 
    test-data/four_human_proteins.fasta
    test-data/four_human_proteins.clinod-1.3.tabular
    tools/clinod/README.rst
    tools/clinod/clinod.xml
    tools/clinod/tool_dependencies.xml


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

NOTE: This is the licence for the Galaxy Wrapper only. Command line
NoD is available and licenced separately.
