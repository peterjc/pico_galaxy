.. image:: https://img.shields.io/travis/peterjc/pico_galaxy/master.svg
   :alt: Linux testing with TravisCI
   :target: https://travis-ci.org/peterjc/pico_galaxy/branches
.. image:: https://landscape.io/github/peterjc/pico_galaxy/master/landscape.svg?style=flat
   :alt: Landscape Code Metrics
   :target: https://landscape.io/github/peterjc/pico_galaxy

Introduction
============

Galaxy is a web-based platform for biological data analysis, supporting
extension with additional tools (often wrappers for existing command line
tools) and datatypes. See http://www.galaxyproject.org/ and the public
server at http://usegalaxy.org for an example.

This repository is for the development of several of my Galaxy tools and
wrappers for use within Galaxy, with a focus on sequence analysis. These
Galaxy tools are released on the Galaxy Tool Shed here:

* http://toolshed.g2.bx.psu.edu/view/peterjc

Previews of works in progress are often on the Galaxy Test Tool Shed,
but these should not be used in a production Galaxy server:

* http://testtoolshed.g2.bx.psu.edu/view/peterjc

Most of these Galaxy tools and wrappers were originally developed on
the 'tools' branch of my BitBucket hosted mercurial (hg) repository,
https://bitbucket.org/peterjc/galaxy-central - the old commit history
has been preserved were possible, and/or the releases made to the
Galaxy Tool Shed were used as snapshots. Some of the content was first
developed on https://github.com/peterjc/picobio and also moved here.

See also https://github.com/peterjc/galaxy_blast for my repository for
the Galaxy wrappers for the NCBI BLAST+ suite and other related tools
like Blast2GO.

The repository name (pico_galaxy) is a play on pico meaning small (10^-12),
and the Japanense phonetics of my name (romaji).


Folder Structure
================

Within the ``tools`` folder is one folder for each Tool or Tool Suite released
on the Galaxy Tool Shed, similarly for the ``workflows`` folder.

Additionally there is a shared ``test-data`` folder used for functional test
sample data, and a shared ``tool-data`` folder used for configuration files.


Testing
=======

Most of these Galaxy tools include a ``<tests>`` section in the tool XML files,
which defines one or more functional tests - listing sample input files and
user parameters, along with the expected output. If you install the tools,
you can run these tests via Galaxy's ``run_functional_tests.sh`` script -
and/or do this automatically if installing the tools via the Tool Shed.

The Galaxy team run nightly tests on all the tools which have been uploaded
to the main Tool Shed and the Test Tool Shed, simulating how they would
behave in a local Galaxy instance once installed via the Tool Shed.

In addition we are running the same functional tests via TravisCI whenever
this GitHub repository is updated:

.. image:: https://travis-ci.org/peterjc/pico_galaxy.png?branch=master
   :alt: Current status of TravisCI build for master branch
   :target: https://travis-ci.org/peterjc/pico_galaxy/builds

This TravisCI integration is still somewhat experimental, but simulates a
manual install of these Galaxy Tools and their dependencies. See the
special ``.travis.yml`` file for more technical details.


Bug Reports
===========

You can file an issue here https://github.com/peterjc/pico_galaxy/issues or
for more general Galaxy Tool Shed problems please ask on the Galaxy development
list http://lists.bx.psu.edu/listinfo/galaxy-dev


License
=======

Please see the README file in each folder, but by default the MIT license is
being used.
