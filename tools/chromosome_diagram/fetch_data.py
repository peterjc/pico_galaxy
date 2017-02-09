#!/usr/bin/env python
import os
import sys

data = [("I", "NC_003070"),
        ("II", "NC_003071"),
        ("III", "NC_003074"),
        ("IV", "NC_003075"),
        ("V", "NC_003076")]

merged_fasta = "A_thaliana.fasta"
merged_trna = "A_thaliana_tRNA.tsv"

BLUE = ("blue", "none")
RED = ("red", "none")
GREEN = ("green", "none")
YELLOW = ("black", "yellow")

# hydrophobicity scheme
colors = {
    # AGILPV --> Blue
    "Ala": BLUE, "Gly": BLUE, "Ile": BLUE, "Leu": BLUE, "Pro": BLUE, "Val": BLUE,
    # FYW --> Red
    "Phe": RED, "Tyr": RED, "Trp": RED,
    # DENQRHSTK --> Green
    "Asp": GREEN, "Glu": GREEN, "Asn": GREEN, "Gln": GREEN, "Arg": GREEN,
    "His": GREEN, "Ser": GREEN, "Thr": GREEN, "Lys": GREEN,
    # CM --> Yellow
    "Cys": YELLOW, "Met": YELLOW,
}

protein_letters_1to3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp',
    'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met',
    'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp',
    'Y': 'Tyr',
}


def run(cmd):
    print cmd
    err = os.system(cmd)
    if err:
        sys.stderr.write("Error %i from command:\n%s\n" % (err, cmd))
        sys.exit(err)


for chr, acc in data:
    for ext in ["fna", "rnt"]:
        filename = "%s.%s" % (acc, ext)
        url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/Arabidopsis_thaliana/CHR_%s/%s" % (chr, filename)
        if not os.path.isfile(filename):
            run("curl -O %s" % url)

if not os.path.isfile(merged_fasta):
    print("Merging FASTA...")
    handle = open(merged_fasta, "w")
    for chr, acc in data:
        filename = "%s.fna" % acc
        with open(filename) as h:
            line = h.readline()
            assert line.startswith(">"), line
            # Rename to chrI etc
            handle.write(">chr%s %s" % (chr, line[1:]))
            print(line.rstrip())
            for line in h:
                assert not line.startswith(">"), line
                handle.write(line)
    handle.close()
    print("Merged FASTA")

if not os.path.isfile(merged_trna):
    print("Merging tRNA...")
    handle = open(merged_trna, "w")
    handle.write("#Name\tChr\tStart\tEnd\tStrand\tProduct\tForeground\tBackground\n")
    for chr, acc in data:
        filename = "%s.rnt" % acc
        with open(filename) as h:
            for line in h:
                if not line.endswith(" tRNA\n"):
                    continue
                fields = line.rstrip("\n").split("\t")
                # Location\tStrand\tLength\tPID\tGene\tSynonymCode\tCOG\tProduct\n
                start, end = fields[0].split("..")
                strand = fields[1]
                name = fields[5]
                product = fields[8]
                assert product.endswith(" tRNA")
                product = product[:-5].strip()
                fore, back = colors.get(product, (".", "."))
                handle.write("\t".join([name, "chr" + chr, start, end, strand, product, fore, back]) + "\n")
    handle.close()
    print("Merged tRNA")

print("Done, see %s and %s" % (merged_fasta, merged_trna))
print("")
print("Try producing a PDF by running something like this at the command line:")
print("")
print('$ ./chromosome_diagram.py A_thaliana.fasta 5000 A_thaliana_tRNA.tsv 2 3 4 5 1 7 8 "A. thaliana tRNA" 0 A_thaliana.pdf')
