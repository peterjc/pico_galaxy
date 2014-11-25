#!/usr/bin/env python
"""Wrapper for EffectiveT3 v1.0.1 for use in Galaxy.

This script takes exactly five command line arguments:
 * model name (e.g. TTSS_STD-1.0.1.jar)
 * threshold (selective or sensitive)
 * an input protein FASTA filename
 * output tabular filename

It then calls the standalone Effective T3 v1.0.1 program (not the
webservice), and reformats the semi-colon separated output into
tab separated output for use in Galaxy.
"""
import sys
import os
import subprocess

#The Galaxy auto-install via tool_dependencies.xml will set this environment variable
effectiveT3_dir = os.environ.get("EFFECTIVET3", "/opt/EffectiveT3/")
effectiveT3_jar = os.path.join(effectiveT3_dir, "TTSS_GUI-1.0.1.jar")

if "-v" in sys.argv or "--version" in sys.argv:
    #TODO - Get version of the JAR file dynamically?
    print("Wrapper v0.0.14, TTSS_GUI-1.0.1.jar")
    sys.exit(0)

def stop_err(msg, error_level=1):
   """Print error message to stdout and quit with given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)

if len(sys.argv) != 5:
   stop_err("Require four arguments: model, threshold, input protein FASTA file & output tabular file")

model, threshold, fasta_file, tabular_file = sys.argv[1:]

if not os.path.isfile(fasta_file):
   stop_err("Input FASTA file not found: %s" % fasta_file)

if threshold not in ["selective", "sensitive"] \
and not threshold.startswith("cutoff="):
   stop_err("Threshold should be selective, sensitive, or cutoff=..., not %r" % threshold)

def clean_tabular(raw_handle, out_handle):
    """Clean up Effective T3 output to make it tabular."""
    count = 0
    positive = 0
    errors = 0
    for line in raw_handle:
        if not line or line.startswith("#") \
        or line.startswith("Id; Description; Score;"):
            continue
        assert line.count(";") >= 3, repr(line)
        #Normally there will just be three semi-colons, however the
        #original FASTA file's ID or description might have had
        #semi-colons in it as well, hence the following hackery:
        try:
            id_descr, score, effective = line.rstrip("\r\n").rsplit(";",2)
            #Cope when there was no FASTA description
            if "; " not in id_descr and id_descr.endswith(";"):
                id = id_descr[:-1]
                descr = ""
            else:
                id, descr = id_descr.split("; ",1)
        except ValueError:
            stop_err("Problem parsing line:\n%s\n" % line)
        parts = [s.strip() for s in [id, descr, score, effective]]
        out_handle.write("\t".join(parts) + "\n")
        count += 1
        if float(score) < 0:
           errors += 1
        if effective.lower() == "true":
           positive += 1
    return count, positive, errors

def run(cmd):
    #Avoid using shell=True when we call subprocess to ensure if the Python
    #script is killed, so too is the child process.
    try:
        child = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception, err:
        stop_err("Error invoking command:\n%s\n\n%s\n" % (" ".join(cmd), err))
    #Use .communicate as can get deadlocks with .wait(),
    stdout, stderr = child.communicate()
    return_code = child.returncode
    if return_code:
        cmd_str= " ".join(cmd) # doesn't quote spaces etc
        if stderr and stdout:
            stop_err("Return code %i from command:\n%s\n\n%s\n\n%s" % (return_code, cmd_str, stdout, stderr))
        else:
            stop_err("Return code %i from command:\n%s\n%s" % (return_code, cmd_str, stderr))

if not os.path.isdir(effectiveT3_dir):
    stop_err("Effective T3 folder not found: %r" % effectiveT3_dir)

if not os.path.isfile(effectiveT3_jar):
    stop_err("Effective T3 JAR file not found: %r" % effectiveT3_jar)

if not os.path.isdir(os.path.join(effectiveT3_dir, "module")):
    stop_err("Effective T3 module folder not found: %r" % os.path.join(effectiveT3_dir, "module"))

effectiveT3_model = os.path.join(effectiveT3_dir, "module", model)
if not os.path.isfile(effectiveT3_model):
    sys.stderr.write("Contents of %r is %s\n"
                     % (os.path.join(effectiveT3_dir, "module"),
                        ", ".join(repr(p) for p in os.listdir(os.path.join(effectiveT3_dir, "module")))))
    sys.stderr.write("Main JAR was found: %r\n" % effectiveT3_jar)
    stop_err("Effective T3 model JAR file not found: %r" % effectiveT3_model)

#We will have write access whereever the output should be,
temp_file = os.path.abspath(tabular_file + ".tmp")

#Use absolute paths since will change current directory...
tabular_file = os.path.abspath(tabular_file)
fasta_file = os.path.abspath(fasta_file)

cmd = ["java", "-jar", effectiveT3_jar,
       "-f", fasta_file,
       "-m", model,
       "-t", threshold,
       "-o", temp_file,
       "-q"]

try:
    #Must run from directory above the module subfolder:
    os.chdir(effectiveT3_dir)
except:
    stop_err("Could not change to Effective T3 folder: %s" % effectiveT3_dir)

run(cmd)

if not os.path.isfile(temp_file):
   stop_err("ERROR - No output file from Effective T3")

out_handle = open(tabular_file, "w")
out_handle.write("#ID\tDescription\tScore\tEffective\n")
data_handle = open(temp_file)
count, positive, errors = clean_tabular(data_handle, out_handle)
data_handle.close()
out_handle.close()

os.remove(temp_file)

if errors:
   print("%i sequences, %i positive, %i errors"
         % (count, positive, errors))
else:
   print("%i/%i sequences positive" % (positive, count))

if count and count==errors:
   #Galaxy will still  allow them to see the output file
   stop_err("All your sequences gave an error code")
