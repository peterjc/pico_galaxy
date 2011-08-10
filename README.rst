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


Bug Reports
===========

You can file an issue here https://github.com/peterjc/pico_galaxy/issues or
for more general Galaxy Tool Shed problems please ask on the Galaxy development
list http://lists.bx.psu.edu/listinfo/galaxy-dev


License
=======

Please see the README file in each folder, but by default the MIT license is
being used.
