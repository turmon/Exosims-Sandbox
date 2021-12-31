#!/usr/bin/env bash
#
# sandbox-test-run.sh: perform a test run to exercise sandbox functions
#
# Usage:
#   sandbox-test-run.sh [ -L ] [ -R ] SCRIPT
#
# where SCRIPT is a script name (e.g., Scripts/HabEx-foobar.json).
#
# By default, the SCRIPT is copied to a temporary file (in the above
# example, Scripts/HabEx-foobar-testrun.json) before running tests,
# so that the test runs land in sims/HabEx-foobar-testrun/... rather than a
# possibly-useful existing directory.
#
# This means that SCRIPT can be an existing script name, and a test run
# will not pollute the corresponding sims/... files.
# 
# Options are:
#   -L: use the SCRIPT literally. This does *not* use the copied script-file
#       scheme above, so the new sims and plots (if they succeed) 
#       will land in sims/SCRIPT/... among whatever might be already there.
#   -R: remove the simulation directory corresponding to the given SCRIPT, before
#       the test run, so the run is clean (except for EXOSIMS file caches).
#       This code prompts for confirmation of -R.
#
# Note that using -R and -L *together* will delete the possibly-existing sims/SCRIPT
# directory, and all its contents.
#
# Description:
#
# The test sequence exercises several sandbox functions within the processing
# pipeline. This is approximately:
#   - single test run  (add-sims.sh ... =777)
#   - ensemble of runs (add-sims.sh ... 100)
# and then various invocations of "make":
#   - reduce data (make reduce)
#   - basic plots (make html)
#   - movies and final frames
#   - webpage update (make html-only)
#   - other ancillary graphics [*]
#   - final webpage update (make html-only)
# The intent is to shake out any errors in the infrastructure by exercising most
# of the mainline functions, all at once. Note that the ancillary graphics [*] above
# can be fragile for some scripts, and they need a 3-year mission.
# 
# See also: python-version-check.py 
#
# turmon dec 2021
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# attempt to give group-write to created files
umask 002

literal=0
remove=0
while getopts "hLR" opt; do
    case $opt in
	h)
	    # help text
	    sed 's/^#//' $(which $0) | awk '/^#/{exit};NR>1{print}'
	    exit 2
	    ;;
	L)
	    # use script literally, do not copy
	    literal=1
	    ;;
	R)
	    # remove results dir in sims/...
	    remove=1
	    ;;
	\?)
	    echo "${PROGNAME}: Invalid option, exiting.  Try -h." >&2
	    exit 2
	    ;;
    esac
done
shift $((OPTIND-1))

# enforce argument count
if [ $# -ne 1 ]; then
   echo "${PROGNAME}: Error: Script filename required." >&2
   exit 1
fi

# original, as-supplied, script name
script_orig="$1"

# check for existence of the given script_orig
if [ ! -r "$script_orig" ]; then
    echo "${PROGNAME}: Error: Given script \`$script_orig' not readable, exiting." >&2
    exit 1
fi    

################################################################################

## Generate the script to use

# should we make a copy or use script_orig literally
if [ $literal -eq 1 ]; then
    echo "${PROGNAME}: Using the supplied script \`$script_orig' for this test."
    script="$script_orig"
else
    # generate new script name
    script=Scripts/$(basename "$script_orig" .json)-testrun.json
    cp -f "$script_orig" "$script"
    echo "${PROGNAME}: Copying the supplied script \`$script_orig' for the test."
fi

echo ""
echo "${PROGNAME}: Test will use script: \`$script'"
if [ $literal -eq 0 ]; then
    echo "${PROGNAME}: You may delete this file any time after the test."
fi

# check for existence of script
if [ ! -r "$script" ]; then
    echo "${PROGNAME}: Error: Script \`$script' not readable, exiting (should not happen)." >&2
    exit 2
fi    

## Process the simulation directory

# the name
sim_dir=sims/$(basename $script .json)

# guard against some fail of basename
if [ $sim_dir == sims/ ]; then
    echo "${PROGNAME}: Name mangling failed, should not happen." >&2
    exit 2
fi

## Empty the sim directory if needed
if [ $remove -eq 1 ]; then
    # in this case, work backward to ensure the file is there, so we know name-mangling worked
    test_script=Scripts/$(basename $sim_dir).json
    if [ ! -r $test_script ]; then
	echo "${PROGNAME}: Name mangling failed, should not happen." >&2
	echo "${PROGNAME}: We think the simulation directory is \`$sim_dir' but cannot find \`$test_script'"
	echo "${PROGNAME}: You must clear the sim dir yourself." >&2
	exit 2
    fi
    # ok to proceed
    echo ""
    echo "${PROGNAME}: Removing the simulation directory \`$sim_dir', as directed."
    echo "${PROGNAME}: Directory contains $(ls $sim_dir | wc -l) files/directories."
    read -p "${PROGNAME}: Answer yes to permanently remove \`$sim_dir' "
    if [ "$REPLY" != "yes" ]; then
	echo "${PROGNAME}: Aborted."
	exit 2
    fi
    # a dangerous command
    rm -rf "$sim_dir"
fi
mkdir -p "$sim_dir"

# check for existence of sim_dir
if [ ! -d "$sim_dir" ]; then
    echo "${PROGNAME}: Error: Simulation directory \`$sim_dir' cannot be created. Exiting." >&2
    exit 1
fi    

echo ""
echo "${PROGNAME}: Test results will go into: \`$sim_dir'"
if [ $literal -eq 0 ]; then
    echo "${PROGNAME}: You may delete this directory any time after the test."
fi
echo ""

################################################################################
# announce progress of this script

SEP="####"
PREV_STEP=
PREV_TIME=
N_STEPS=10

function announce_progress() {
    if [ -n "$PREV_STEP" ]; then
	DELTA=$(echo $(date +%s) - $PREV_TIME | bc)
	echo "$SEP"
	echo "$SEP Finished $PREV_STEP (took $DELTA seconds)"
	echo ""
    fi
    echo "$SEP"
    echo "$SEP Phase: $2"
    echo "$SEP Step $1 of $N_STEPS"
    echo "$SEP"
    PREV_TIME=$(date "+%s")
    PREV_STEP=$2
}

################################################################################
# check for existence of some programs

function die() {
    echo "${PROGNAME}: Abort: $1" 1>&2
    exit 2
}

function check_programs() {
    which bash || die 'bash not found'
    /usr/bin/env bash --version | head -1
    which python || die 'python not found'
    python --version || die "No python"
    which parallel || die 'parallel not found'
    parallel --version | head -1
    which matlab || die 'matlab not found'
    echo "exit" | matlab -nodisplay -nosplash -nodesktop |& grep -v '^$' | head -4
    }

################################################################################
# check that we can ssh without password, for parallel runs

REMOTES="aftac1 aftac2 aftac3 mustang2 mustang3 mustang4"

function check_ssh() {
    ssh_fails=
    for r in $REMOTES; do
	if ! ssh -o PasswordAuthentication=no -o BatchMode=yes $r exit &> /dev/null; then
	    echo "passwordless ssh to $r failed"
	    ssh_fails="$ssh_fails $r"
	else
	    echo "passwordless ssh to $r OK"
	fi
    done
    }

################################################################################
# main sequence of checks

announce_progress 1 "Auxiliary Programs"
check_programs

announce_progress 2 "Passwordless SSH"
check_ssh $script
if [ -n "$ssh_fails" ]; then
    echo "${PROGNAME}: ssh failed to $ssh_fails -- parallel add-sims will fail"
fi

announce_progress 3 "Initial add-sims"
add-sims.sh $script =777
announce_progress 4 "Parallel add-sims"
add-sims.sh -Z $script 100

# the target for make
script_base=$(basename $script .json)

announce_progress 4 "Reductions"
make S=$script_base reduce
announce_progress 5 "Plots (make html)"
make S=$script_base html

announce_progress 6 "Path Ensemble"
make -j 10 S=$script_base path-ensemble
announce_progress 7a "Path Movies"
make -j 5 S=$script_base path-movie-5
announce_progress 7b "Path Final Frames"
make -j 10 S=$script_base path-final-10

announce_progress 8 "HTML page re-make"
make S=$script_base html-only

announce_progress 9a "Observing-target timelines"
make -j 5 S=$script_base obs-timeline-5
announce_progress 9b "Keepout vs. time"
make -j 5 S=$script_base keepout-5

announce_progress 10 "HTML page re-make (final)"
make S=$script_base html-only

echo ""
echo "${PROGNAME}: Finished running ($(date))."
echo "${PROGNAME}: All checks finished OK."


