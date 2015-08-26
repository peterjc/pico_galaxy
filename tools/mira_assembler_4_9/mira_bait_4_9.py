#!/usr/bin/env python
"""A simple wrapper script to call MIRA4's mirabait and collect its output.
"""
import os
import sys
import subprocess
import shutil
import time

WRAPPER_VER = "0.0.1" #Keep in sync with the XML file

def sys_exit(msg, err=1):
    sys.stderr.write(msg+"\n")
    sys.exit(err)


def get_version(mira_binary):
    """Run MIRA to find its version number"""
    # At the commend line I would use: mira -v | head -n 1
    # however there is some pipe error when doing that here.
    cmd = [mira_binary, "-v"]
    try:
        child = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
    except Exception, err:
        sys.stderr.write("Error invoking command:\n%s\n\n%s\n" % (" ".join(cmd), err))
        sys.exit(1)
    ver, tmp = child.communicate()
    del child
    #Workaround for -v not working in mirabait 4.0RC4
    if "invalid option" in ver.split("\n", 1)[0]:
        for line in ver.split("\n", 1):
            if " version " in line:
                line = line.split()
                return line[line.index("version")+1].rstrip(")")
        sys_exit("Could not determine MIRA version:\n%s" % ver)
    return ver.split("\n", 1)[0]

try:
    mira_path = os.environ["MIRA4"]
except KeyError:
    sys_exit("Environment variable $MIRA4 not set")
mira_binary = os.path.join(mira_path, "mirabait")
if not os.path.isfile(mira_binary):
    sys_exit("Missing mirabait under $MIRA4, %r\nFolder contained: %s"
             % (mira_binary, ", ".join(os.listdir(mira_path))))
mira_ver = get_version(mira_binary)
if not mira_ver.strip().startswith("4.9"):
    sys_exit("This wrapper is for MIRA V4.9, not:\n%s" % mira_ver)
if "-v" in sys.argv or "--version" in sys.argv:
    print "%s, MIRA wrapper version %s" % (mira_ver, WRAPPER_VER)
    sys.exit(0)


format, output_choice, strand_choice, kmer_length, min_occurance, bait_file, in_file, out_file = sys.argv[1:]

if format.startswith("fastq"):
    format = "fastq"
elif format == "mira":
    format = "maf"
elif format != "fasta":
    sys_exit("Was not expected format %r" % format)

assert out_file.endswith(".dat")
out_file_stem = out_file[:-4]

cmd_list = [mira_binary, "-f", format, "-t", format,
            "-k", kmer_length, "-n", min_occurance,
            bait_file, in_file, out_file_stem]
if output_choice == "pos":
    pass
elif output_choice == "neg":
    #Invert the selection...
    cmd_list.insert(1, "-i")
else:
    sys_exit("Output choice should be 'pos' or 'neg', not %r" % output_choice)
if strand_choice == "both":
    pass
elif strand_choice == "fwd":
    #Ingore reverse strand...
    cmd_list.insert(1, "-r")
else:
    sys_exit("Strand choice should be 'both' or 'fwd', not %r" % strand_choice)

cmd = " ".join(cmd_list)
#print cmd
start_time = time.time()
try:
    #Run MIRA
    child = subprocess.Popen(cmd_list,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
except Exception, err:
    sys.stderr.write("Error invoking command:\n%s\n\n%s\n" % (cmd, err))
    sys.exit(1)
#Use .communicate as can get deadlocks with .wait(),
stdout, stderr = child.communicate()
assert stderr is None # Due to way we ran with subprocess
run_time = time.time() - start_time
return_code = child.returncode
print "mirabait took %0.2f minutes" % (run_time / 60.0)

if return_code:
    sys.stderr.write(stdout)
    sys_exit("Return error code %i from command:\n%s" % (return_code, cmd),
             return_code)

#Capture output
out_tmp = out_file_stem + "." + format
if not os.path.isfile(out_tmp):
    sys.stderr.write(stdout)
    sys_exit("Missing output file from mirabait: %s" % out_tmp)
shutil.move(out_tmp, out_file)
