Galaxy tool to sub-sample sequence files
========================================

This tool is copyright 2014-2014 by Peter Cock, The James Hutton Institute
(formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
See the licence text below (MIT licence).

This tool is a short Python script (using Biopython library functions)
to sub-sample sequence files (in a range of formats including FASTA, FASTQ,
and SFF). This can be useful for preparing a small sample of data to test
or time a new pipeline, or for reducing the read coverage in a de novo
assembly.

This tool is available from the Galaxy Tool Shed at:

* http://toolshed.g2.bx.psu.edu/view/peterjc/sample_seqs


Automated Installation
======================

This should be straightforward using the Galaxy Tool Shed, which should be
able to automatically install the dependency on Biopython, and then install
this tool and run its unit tests.


Manual Installation
===================

There are just two files to install to use this tool from within Galaxy:

* ``sample_seqs.py`` (the Python script)
* ``sample_seqs.xml`` (the Galaxy tool definition)

The suggested location is in a dedicated ``tools/sample_seqs`` folder.

You will also need to modify the ``tools_conf.xml`` file to tell Galaxy to offer the
tool. One suggested location is in the filters section. Simply add the line::

    <tool file="sample_seqs/sample_seqs.xml" />

You will also need to install Biopython 1.62 or later.

If you wish to run the unit tests, also	move/copy the ``test-data/`` files
under Galaxy's ``test-data/`` folder. Then::

    ./run_tests.sh -id sample_seqs

That's it.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - Initial version.
v0.1.1  - Using optparse to provide a proper command line API.
v0.1.2  - Interleaved mode for working with paired records.
        - Tool definition now embeds citation information.
v0.2.0  - Option to give number of sequences (or pairs) desired.
          This works by first counting all your sequences, then calculates
          the percentage required in order to sample them uniformly (evenly).
          This makes two passes through the input and is therefore slower. 
v0.2.1  - Was missing a file for the functional tests.
        - Included testing of stdout messages.
        - Includes testing of failure modes.
======= ======================================================================


Developers
==========

This script and related tools are being developed on this GitHub repository:
https://github.com/peterjc/pico_galaxy/tree/master/tools/sample_seqs

For making the "Galaxy Tool Shed" http://toolshed.g2.bx.psu.edu/ tarball use
the following command from the Galaxy root folder::

    $ tar -czf sample_seqs.tar.gz tools/sample_seqs/README.rst tools/sample_seqs/sample_seqs.py tools/sample_seqs/sample_seqs.xml tools/sample_seqs/tool_dependencies.xml test-data/ecoli.fastq test-data/ecoli.sample_N100.fastq test-data/ecoli.pair_sample_N100.fastq test-data/ecoli.sample_C10.fastq test-data/get_orf_input.Suis_ORF.prot.fasta test-data/get_orf_input.Suis_ORF.prot.sample_N100.fasta test-data/get_orf_input.Suis_ORF.prot.pair_sample_N100.fasta test-data/get_orf_input.Suis_ORF.prot.sample_C10.fasta test-data/get_orf_input.Suis_ORF.prot.pair_sample_C10.fasta test-data/MID4_GLZRM4E04_rnd30_frclip.sff test-data/MID4_GLZRM4E04_rnd30_frclip.sample_N5.sff test-data/MID4_GLZRM4E04_rnd30_frclip.pair_sample_N5.sff test-data/MID4_GLZRM4E04_rnd30_frclip.sample_C1.sff

Check this worked::

    $ tar -tzf sample_seqs.tar.gz
    tools/sample_seqs/README.rst
    tools/sample_seqs/sample_seqs.py
    tools/sample_seqs/sample_seqs.xml
    tools/sample_seqs/tool_dependencies.xml
    test-data/ecoli.fastq
    test-data/ecoli.sample_N100.fastq
    test-data/ecoli.pair_sample_N100.fastq
    test-data/ecoli.sample_C10.fastq
    test-data/get_orf_input.Suis_ORF.prot.fasta
    test-data/get_orf_input.Suis_ORF.prot.sample_N100.fasta
    test-data/get_orf_input.Suis_ORF.prot.pair_sample_N100.fasta
    test-data/get_orf_input.Suis_ORF.prot.sample_C10.fasta
    test-data/get_orf_input.Suis_ORF.prot.pair_sample_C10.fasta
    test-data/MID4_GLZRM4E04_rnd30_frclip.sff
    test-data/MID4_GLZRM4E04_rnd30_frclip.sample_N5.sff
    test-data/MID4_GLZRM4E04_rnd30_pair_sample.sff
    test-data/MID4_GLZRM4E04_rnd30_frclip.pair_sample_N5.sff
    test-data/MID4_GLZRM4E04_rnd30_frclip.sample_C1.sff


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
