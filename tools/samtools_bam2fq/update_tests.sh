#!/bin/sh
set -e
echo "This will update test files using the current version of BLAST+"

if [ -f "tools/samtools_bam2fq/update_tests.sh" ]
then
echo "Good, in the expected directory"
else
echo "ERROR. Run this from the GitHub repository root directory."
exit 1
fi

cd test-data

echo "sam_spec_padded.bam2fq.fastq (defaults)"
samtools bam2fq sam_spec_padded.bam > sam_spec_padded.bam2fq.fastq
# Should give idential output to:
# $ samtools bam2fq sam_spec_padded.sam > sam_spec_padded.sam2fq.fastq
# $ samtools bam2fq sam_spec_padded.depad.bam > sam_spec_padded.depad.bam2fq.fastq

echo "sam_spec_padded.bam2fq_no_suf.fastq (without suffices)"
samtools bam2fq -n sam_spec_padded.bam > sam_spec_padded.bam2fq_no_suf.fastq
