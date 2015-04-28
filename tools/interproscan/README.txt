Created 07/01/2011 - Konrad Paszkiewicz, Exeter Sequencing Service, University of Exeter, UK
Revisions 2013-5015 by Peter Cock, The James Hutton Institute, UK

The attached is a crude wrapper script for Interproscan. Typically this is useful when one wants to produce an annotation which is not based on sequence 
similarity. E.g after a denovo transcriptome assembly, each transcript could be translated and run through this tool.

Prerequisites:

1. A working installation of Interproscan on your Galaxy server/cluster.

Limitations:

Currently it is setup to work with PFAM only due to the heavy computational demands Interproscan makes. 

Input formats:

The standard interproscan input is either genomic or protein sequences. In the case of genomic sequences Interproscan will of run an ORF 
prediction tool. However this tends to lose the ORF information (e.g. start/end co-ordinates) from the header. As such the requirement here is to input ORF 
sequences (e.g. from EMBOSS getorf) and to then replace any spaces in the FASTA header with underscores. This workaround generally preserves the relevant 
positional information. 



