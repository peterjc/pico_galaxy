#!/usr/bin/env python
"""Wrapper script using miraconvert & samtools to get BAM from MIRA.
"""
import os
import sys
import shutil
import tempfile

def stop_err(msg, err=1):
    sys.stderr.write(msg+"\n")
    sys.exit(err)

def make_bam(mira_convert, maf_file, fasta_file, bam_file):
    if not os.path.isfile(mira_convert):
        return "Missing binary %r" % mira_convert
    if not os.path.isfile(maf_file):
        return "Missing input MIRA file: %r" % maf_file
    if not os.path.isfile(fasta_file):
        return "Missing padded FASTA file: %r" % fasta_file
    #Debug:
    print("Converting %s (MIRA) to %s (BAM) using reference %s (FASTA)" % (maf_file, bam_file, fasta_file))

    tmp_dir = tempfile.mkdtemp()
    sam_file = os.path.join(tmp_dir, "x.sam")

    # Note add nbb to the template name, possible MIRA 4.0 RC4 bug
    cmd = '"%s" -f maf -t samnbb "%s" "%snbb"' % (mira_convert, maf_file, sam_file)
    return_code = os.system(cmd)
    if return_code:
        return "Error %i from command:\n%s" % (return_code, cmd)
    if not os.path.isfile(sam_file):
        return "Conversion from MIRA to SAM failed"

    #Also doing SAM to (uncompressed) BAM during depad
    bam_stem = bam_file + ".tmp" # Have write permissions and want final file in this folder
    cmd = 'samtools depad -S -u -T "%s" "%s" | samtools sort - "%s"' % (fasta_file, sam_file, bam_stem)
    return_code = os.system(cmd)
    if return_code:
        return "Error %i from command:\n%s" % (return_code, cmd)
    if not os.path.isfile(bam_stem + ".bam"):
        return "samtools depad or sort failed to produce BAM file"

    os.remove(sam_file)
    os.rmdir(tmp_dir)

    cmd = 'samtools index "%s.bam"' % bam_stem
    return_code = os.system(cmd)
    if return_code:
        return "Error %i from command:\n%s" % (return_code, cmd)
    if not os.path.isfile(bam_stem + ".bam.bai"):
        return "samtools indexing of BAM file failed to produce BAI file"

    shutil.move(bam_stem + ".bam", bam_file)
    os.remove(bam_stem + ".bam.bai") #Let Galaxy handle that...

    return None #Good :)

if __name__ == "__main__":
    mira_convert, maf_file, fasta_file, bam_file = sys.argv[1:]
    msg = make_bam(mira_convert, maf_file, fasta_file, bam_file)
    if msg:
        stop_err(msg)
