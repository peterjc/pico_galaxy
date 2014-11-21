#!/bin/sh
set -e
echo "This will update test files using the current version of samtools and the script"

if [ -f "tools/coverage_stats/update_tests.sh" ]
then
echo "Good, in the expected directory"
else
echo "ERROR. Run this from the GitHub repository root directory."
exit 1
fi

cd test-data

echo ex1.coverage_stats.tabular
../tools/coverage_stats/coverage_stats.py ex1.bam ex1.bam.bai ex1.coverage_stats.tabular 

echo coverage_test.bam
samtools view -b coverage_test.sam > coverage_test.bam
samtools index coverage_test.bam

echo coverage_test.coverage_stats.tabular
../tools/coverage_stats/coverage_stats.py coverage_test.bam coverage_test.bam.bai coverage_test.coverage_stats.tabular
