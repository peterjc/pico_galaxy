Python re-implementation of predictNLS with Galaxy wrapper
==========================================================

This Galaxy tool is copyright 2011-2015 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below.

The tool consists of a Galaxy interface definition (``predictnls.xml``), and a Python
script (``predictnls.py``) which re-implements the command line tool predictNLS. This
should match the behaviour of predictNLS v1.0.20 (July 2011), the current latest
release from the Rost Lab, see http://rostlab.org and their paper:

Murat Cokol, Rajesh Nair, and Burkhard Rost.
Finding nuclear localization signals.
EMBO reports 1(5), 411–415, 2000
http://dx.doi.org/10.1093/embo-reports/kvd092

This wrapper is available from the Galaxy Tool Shed at
http://toolshed.g2.bx.psu.edu/view/peterjc/predictnls


Automatic Installation
======================

This Galaxy tool is self contained, and so should install automatically via the
Galaxy Tool Shed.


Manual Installation
===================

There are just four files which should be moved under the Galaxy tools folder,
e.g. in a ``tools/protein_analysis`` filter:

* ``predictlns.xml`` (the Galaxy tool definition)
* ``predictlns.py`` (the Python script)
* ``predictlns.txt`` (this README file)
* ``My_NLS_list`` (the default set of NLS motifs from the Rost Lab)

You will also need to modify the tools_conf.xml file to tell Galaxy to offer the
tool. If you are using other protein analysis tools like TMHMM or SignalP, put
it next to them. Just add the line::

  <tool file="protein_analysis/predictnls.xml" />

If you wish to run the unit tests, also move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    ./run_tests.sh -id predictnls

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.4  - Initial public release.
v0.0.5  - Treat non-zero return codes as errors.
v0.0.6  - Link to Tool Shed added to help text and this documentation.
        - Use reStructuredText for this README file.
        - Updated citation information (Cock et al. 2013).
        - Development moved to GitHub, https://github.com/peterjc/pico_galaxy
v0.0.7  - Tool definition now embeds citation information.
v0.0.8  - Reorder XML elements (internal change only).
        - Planemo for Tool Shed upload (``.shed.yml``, internal change only).
======= ======================================================================


Developers
==========

This script and related tools were initially developed on the following hg branch:
http://bitbucket.org/peterjc/galaxy-central/src/tools

For pushing a release to the test or main "Galaxy Tool Shed", use the following
Planemo commands (which requires you have set your Tool Shed access details in
``~/.planemo.yml`` and that you have access rights on the Tool Shed)::

    $ planemo shed_update -t testtoolshed --check_diff ~/repositories/pico_galaxy/tools/predictnls/
    ...

or::

    $ planemo shed_update -t toolshed --check_diff ~/repositories/pico_galaxy/tools/predictnls/
    ...

To just build and check the tar ball, use::

    $ planemo shed_upload --tar_only  ~/repositories/pico_galaxy/tools/predictnls/
    ...
    $ tar -tzf shed_upload.tar.gz 
    test-data/four_human_proteins.fasta
    test-data/four_human_proteins.predictnls.tabular
    tools/predictnls/My_NLS_list
    tools/predictnls/README.rst
    tools/predictnls/predictnls.py
    tools/predictnls/predictnls.xml


Licence (GPL)
=============

This tool is open source, licensed under the GNU GENERAL PUBLIC LICENSE
version 3 (GNU v3), see http://www.gnu.org/licenses/gpl.html

The Python script is my reimplementation of the original Perl program from
the Rost Lab, which was released under the GPL v3. Therefore, as I consider
this to be a derivative work, this too is released under the GPL v3.

Please note that the ``My_NLS_list`` file should be an exact copy of the file of the
same name included with ``predictnls-1.0.7.tar.gz`` to ``predictnls-1.0.20.tar.gz``
inclusive (the list was extended in v1.0.7 in August 2010, see the change log
included in those tar-balls), available from ftp://rostlab.org/predictnls/
