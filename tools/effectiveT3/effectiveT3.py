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
import os

# We want to be able to use shutil.which, but need Python 3.3+
# import shutil
import subprocess
import sys

# The Galaxy auto-install via tool_dependencies.xml will set the
# environment variable $EFFECTIVET3 pointing at the folder with
# the JAR file.
#
# The BioConda recipe will put a wrapper script on the $PATH,
# which we can use to find the JAR file.
#
# We fall back on /opt/EffectiveT3/
#
effective_t3_jarname = "TTSS_GUI-1.0.1.jar"

if "-v" in sys.argv or "--version" in sys.argv:
    # TODO - Get version of the JAR file dynamically?
    print("Wrapper v0.0.20, for %s" % effective_t3_jarname)
    sys.exit(0)

if len(sys.argv) != 5:
    sys.exit(
        "Require four arguments: model, threshold, input protein FASTA file & output tabular file"
    )

model, threshold, fasta_file, tabular_file = sys.argv[1:]

if not os.path.isfile(fasta_file):
    sys.exit("Input FASTA file not found: %s" % fasta_file)

if threshold not in ["selective", "sensitive"] and not threshold.startswith("cutoff="):
    sys.exit(
        "Threshold should be selective, sensitive, or cutoff=..., not %r" % threshold
    )


def clean_tabular(raw_handle, out_handle):
    """Clean up Effective T3 output to make it tabular."""
    count = 0
    positive = 0
    errors = 0
    for line in raw_handle:
        if (
            not line
            or line.startswith("#")
            or line.startswith("Id; Description; Score;")
        ):
            continue
        assert line.count(";") >= 3, repr(line)
        # Normally there will just be three semi-colons, however the
        # original FASTA file's ID or description might have had
        # semi-colons in it as well, hence the following hackery:
        try:
            id_descr, score, effective = line.rstrip("\r\n").rsplit(";", 2)
            # Cope when there was no FASTA description
            if "; " not in id_descr and id_descr.endswith(";"):
                id = id_descr[:-1]
                descr = ""
            else:
                id, descr = id_descr.split("; ", 1)
        except ValueError:
            sys.exit("Problem parsing line:\n%s\n" % line)
        parts = [s.strip() for s in [id, descr, score, effective]]
        out_handle.write("\t".join(parts) + "\n")
        count += 1
        if float(score) < 0:
            errors += 1
        if effective.lower() == "true":
            positive += 1
    return count, positive, errors


def run(cmd):
    """Run the command line string via subprocess."""
    # Avoid using shell=True when we call subprocess to ensure if the Python
    # script is killed, so too is the child process.
    try:
        child = subprocess.Popen(
            cmd, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except Exception as err:
        sys.exit("Error invoking command:\n%s\n\n%s\n" % (" ".join(cmd), err))
    # Use .communicate as can get deadlocks with .wait(),
    stdout, stderr = child.communicate()
    return_code = child.returncode
    if return_code or stderr.startswith("Exception in thread"):
        cmd_str = " ".join(cmd)  # doesn't quote spaces etc
        if stderr and stdout:
            sys.exit(
                "Return code %i from command:\n%s\n\n%s\n\n%s"
                % (return_code, cmd_str, stdout, stderr)
            )
        else:
            sys.exit(
                "Return code %i from command:\n%s\n%s" % (return_code, cmd_str, stderr)
            )


try:
    from shutil import which
except ImportError:
    # Likely running on Python 2, use backport:
    def which(cmd, mode=os.F_OK | os.X_OK, path=None):
        """Python implementation of command line tool which.

        Given a command, mode, and a PATH string, return the path which
        conforms to the given mode on the PATH, or None if there is no such
        file.

        `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
        of os.environ.get("PATH"), or can be overridden with a custom search
        path.
        """
        # Check that a given file can be accessed with the correct mode.
        # Additionally check that `file` is not a directory, as on Windows
        # directories pass the os.access check.
        def _access_check(fn, mode):
            return os.path.exists(fn) and os.access(fn, mode) and not os.path.isdir(fn)

        # Short circuit. If we're given a full path which matches the mode
        # and it exists, we're done here.
        if _access_check(cmd, mode):
            return cmd

        path = (path or os.environ.get("PATH", os.defpath)).split(os.pathsep)

        if sys.platform == "win32":
            # The current directory takes precedence on Windows.
            if os.curdir not in path:
                path.insert(0, os.curdir)

            # PATHEXT is necessary to check on Windows.
            pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
            # See if the given file matches any of the expected path extensions.
            # This will allow us to short circuit when given "python.exe".
            matches = [cmd for ext in pathext if cmd.lower().endswith(ext.lower())]
            # If it does match, only test that one, otherwise we have to try
            # others.
            files = [cmd] if matches else [cmd + ext.lower() for ext in pathext]
        else:
            # On other platforms you don't have things like PATHEXT to tell you
            # what file suffixes are executable, so just pass on cmd as-is.
            files = [cmd]

        seen = set()
        for dir in path:
            dir = os.path.normcase(dir)
            if dir not in seen:
                seen.add(dir)
                for thefile in files:
                    name = os.path.join(dir, thefile)
                    if _access_check(name, mode):
                        return name
        return None


# Try in order the following to find the JAR file:
# - Location of any wrapper script, e.g. from BioConda installation
# - The $EFFECTIVET3 env var, e.g. old-style Galaxy tool installation
# - The /opt/EffectiveT3/ folder.
effective_t3_jar = None
effective_t3_dir = None
dirs = ["/opt/EffectiveT3/"]
if "EFFECTIVET3" in os.environ:
    dirs.insert(0, os.environ.get("EFFECTIVET3"))
if which("effectivet3"):
    # Assuming this is a BioConda installed wrapper for effective T3,
    # this will get the directory of the wrapper script which is where
    # the JAR file will be:
    dirs.insert(0, os.path.split(os.path.realpath(which("effectivet3")))[0])
for effective_t3_dir in dirs:
    effective_t3_jar = os.path.join(effective_t3_dir, effective_t3_jarname)
    if os.path.isfile(effective_t3_jar):
        # Good
        break
    effective_t3_jar = None
if not effective_t3_dir or not effective_t3_jar:
    sys.exit("Effective T3 JAR file %r not found in %r" % (effective_t3_jarname, dirs))

if not os.path.isdir(os.path.join(effective_t3_dir, "module")):
    sys.exit(
        "Effective T3 module folder not found: %r"
        % os.path.join(effective_t3_dir, "module")
    )

effective_t3_model = os.path.join(effective_t3_dir, "module", model)
if not os.path.isfile(effective_t3_model):
    sys.stderr.write(
        "Contents of %r is %s\n"
        % (
            os.path.join(effective_t3_dir, "module"),
            ", ".join(
                repr(p) for p in os.listdir(os.path.join(effective_t3_dir, "module"))
            ),
        )
    )
    sys.stderr.write("Main JAR was found: %r\n" % effective_t3_jar)
    sys.exit("Effective T3 model JAR file not found: %r" % effective_t3_model)

# We will have write access wherever the output should be,
if tabular_file == "/dev/stdout":
    temp_file = os.path.abspath("effectivet3_tabular_output.tmp")
else:
    temp_file = os.path.abspath(tabular_file + ".tmp")

# Use absolute paths since will change current directory...
tabular_file = os.path.abspath(tabular_file)
fasta_file = os.path.abspath(fasta_file)

cmd = [
    "java",
    "-jar",
    effective_t3_jar,
    "-f",
    fasta_file,
    "-m",
    model,
    "-t",
    threshold,
    "-o",
    temp_file,
    "-q",
]

try:
    # Must run from directory above the module subfolder:
    os.chdir(effective_t3_dir)
except Exception:
    sys.exit("Could not change to Effective T3 folder: %s" % effective_t3_dir)

run(cmd)

if not os.path.isfile(temp_file):
    sys.exit("ERROR - No output file from Effective T3")

out_handle = open(tabular_file, "w")
out_handle.write("#ID\tDescription\tScore\tEffective\n")
data_handle = open(temp_file)
count, positive, errors = clean_tabular(data_handle, out_handle)
data_handle.close()
out_handle.close()

os.remove(temp_file)

if errors:
    print("%i sequences, %i positive, %i errors" % (count, positive, errors))
else:
    print("%i/%i sequences positive" % (positive, count))

if count and count == errors:
    # Galaxy will still  allow them to see the output file
    sys.exit("All your sequences gave an error code")
